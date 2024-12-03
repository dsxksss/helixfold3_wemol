import hashlib
import json
import logging
import os
from datetime import datetime
import sys
from typing import Dict, List, Tuple, Union


def mkdir_if_not_exist(path: str) -> None:
    if not os.path.exists(path):
        os.makedirs(path)


def parse_timestamp(timestamp, custom_strfmt="%Y-%m-%d %H:%M:%S"):
    # 将时间戳转换为datetime对象
    dt = datetime.fromtimestamp(timestamp)

    # 格式化日期和时间
    formatted = dt.strftime(custom_strfmt)

    return formatted


def set_logging_default_config(
    log_level: int = 15,
    file_log_level: int = 15,
    console_log_level: int = 15,
    log_file_save_div: str = "./logs",
) -> None:
    logging.addLevelName(15, "TEMPLATE_DEBUG")

    console_handler = logging.StreamHandler()
    mkdir_if_not_exist(log_file_save_div)

    file_handler = logging.FileHandler(
        f"./{log_file_save_div}/Helixfold_WeMol.log", encoding="utf-8"
    )
    console_handler.setLevel(console_log_level)
    file_handler.setLevel(file_log_level)

    console_format = logging.Formatter(
        "[%(asctime)s] %(funcName)s - %(levelname)s | %(message)s"
    )
    file_format = logging.Formatter(
        "[%(asctime)s] %(funcName)s - %(levelname)s | %(message)s"
    )

    console_handler.setFormatter(console_format)
    file_handler.setFormatter(file_format)

    logging.basicConfig(level=log_level, handlers=[console_handler, file_handler])


def get_sha256_hash_of_file(file_path):
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        # 读取文件直到结束
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()


def load_fasta(fst: str) -> list:
    """读取FASTA格式的序列文件

    Args:
        fst (str): FASTA文件路径

    Returns:
        list: 包含(序列ID, 序列)元组的列表
    """
    fasta, sid, seq = [], None, ""
    for line in open(fst, "r").readlines():
        line = line.strip()
        if line.startswith(">"):
            if len(seq) > 0:
                fasta.append((sid, seq))
                seq = ""
            sid = line[1:]
        else:
            seq += line.upper()
            continue
    if len(seq) > 0:
        fasta.append((sid, seq))
    return fasta


def load_fasta_batch(fst: str) -> list:
    fasta, sid, seq = [], None, ""
    unique_id = []
    for line in open(fst, "r").readlines():
        line = line.strip()
        if line.startswith(">"):
            if len(seq) > 0:
                fasta.append((sid, seq))
                seq = ""
            sid = line[1:]
            if sid in unique_id:
                sys.stderr.write(
                    f"The seq ID in '{fst}' should be unique in batch format: f{sid}"
                )
            else:
                unique_id.append(sid)
        else:
            seq += line.upper()
            continue
    if len(seq) > 0:
        fasta.append((sid, seq))
    return fasta


def load_ligands(lig_file: str, batch: bool = False) -> list | dict:
    """读取配体信息文件

    Args:
        lig_file: 配体文件路径
        batch: 是否使用批处理模式

    Returns:
        list | dict: 非批处理模式返回配体列表，批处理模式返回字典
    """
    if batch:
        ligands = {}
        with open(lig_file) as fin:
            for line in fin.readlines():
                ligs = []
                tmp = line.strip().split(":")
                uid = tmp[0].strip()
                for lig in tmp[1:]:
                    lig = lig.strip()
                    if lig.startswith("CCD,") or lig.startswith("ccd,"):
                        lig = lig.split(",")
                        ligs.append([t.strip() for t in lig])
                    else:
                        ligs.append(lig)
                ligands[uid] = ligs
    else:
        ligands = []
        with open(lig_file) as fin:
            for line in fin.readlines():
                tmp = line.strip()
                if tmp.startswith("CCD,") or tmp.startswith("ccd,"):
                    tmp = tmp.split(",")
                    ligands.append([t.strip() for t in tmp])
                else:
                    ligands.append(tmp)

    return ligands


def load_modification(mdf_file: str, batch: bool = False) -> list | dict:
    """读取修饰信息文件

    Args:
        mdf_file: 修饰文件路径
        batch: 是否使用批处理模式

    Returns:
        list | dict: 非批处理模式返回修饰列表，批处理模式返回字典
    """
    if batch:
        mdfs_batch = {}
        with open(mdf_file) as fin:
            for line in fin.readlines():
                tmp = line.strip().split(":")
                uid = tmp[0].strip()
                mdfs = []
                for mdf in tmp[1:]:
                    mdf = mdf.strip().split(",")
                    assert len(mdf) == 3, f"Error format in modifications: {line}"
                    seq_index = int(mdf[0])
                    type_mdf = mdf[1].strip()
                    pos_aa = int(mdf[2])
                    assert (
                        seq_index >= 1 and pos_aa >= 1 and len(type_mdf) == 3
                    ), f"Error format in modifications: {line}"
                    mdfs.append((seq_index, type_mdf, pos_aa))
                mdfs_batch[uid] = mdfs
        return mdfs_batch
    else:
        mdfs = []
        with open(mdf_file) as fin:
            for line in fin.readlines():
                tmp = line.strip().split(",")
                assert len(tmp) == 3, f"Error format in modifications: {line}"
                seq_index = int(tmp[0])
                type_mdf = tmp[1].strip()
                pos_aa = int(tmp[2])
                assert (
                    seq_index >= 1 and pos_aa >= 1 and len(type_mdf) == 3
                ), f"Error format in modifications: {line}"
                mdfs.append((seq_index, type_mdf, pos_aa))
        return mdfs


def load_bonds(bond_file: str, batch: bool = False) -> list | dict:
    """读取共价键信息文件

    Args:
        bond_file: 共价键文件路径
        batch: 是否使用批处理模式

    Returns:
        list | dict: 非批处理模式返回键列表，批处理模式返回字典
    """

    def parse_bond(bond_str: str) -> tuple:
        atoms = bond_str.strip().split(";")
        assert len(atoms) == 2, f"Error format in covalent bonds: {bond_str}"
        left = atoms[0].strip().split(",")
        right = atoms[1].strip().split(",")
        assert (
            len(left) == 3 and len(right) == 3
        ), f"Error format in covalent bonds: {bond_str}"

        seq_index_l = int(left[0])
        pos_l = int(left[1])
        atom_name_l = left[2].strip()
        seq_index_r = int(right[0])
        pos_r = int(right[1])
        atom_name_r = right[2].strip()

        assert (
            seq_index_l >= 1 and pos_l >= 1 and seq_index_r >= 1 and pos_r >= 1
        ), f"Error format in covalent bonds: {bond_str}"

        atom_l = (seq_index_l, pos_l, atom_name_l)
        atom_r = (seq_index_r, pos_r, atom_name_r)
        return (atom_l, atom_r)

    if batch:
        bonds_batch = {}
        with open(bond_file) as fin:
            for line in fin.readlines():
                tmp = line.strip().split(":")
                uid = tmp[0].strip()
                bonds = []
                for bond in tmp[1:]:
                    bonds.append(parse_bond(bond))
                bonds_batch[uid] = bonds
        return bonds_batch
    else:
        bonds = []
        with open(bond_file) as fin:
            for line in fin.readlines():
                bonds.append(parse_bond(line))
        return bonds


chain_MARK = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def set_progress_value(value: int):
    if value < 0 or value > 100:
        raise ValueError(
            "ProgressValue must be between 0 and 100, but got [{value}].".format(value)
        )

    with open("./state.json", "w", encoding="utf-8") as f:
        json.dump({"ProgressValue": value}, f)
