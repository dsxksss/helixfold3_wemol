import glob
import logging
import os
import json
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from subprocess import Popen
from typing import Optional
import shutil
import zipfile

from dotenv import load_dotenv

from py_script_template.cli import get_cli_argument
from .utils import (
    set_logging_default_config,
    set_progress_value,
)


def run_ext_cmder(cmds: list, query: Optional[str] = None) -> None:
    """执行外部命令行工具

    Args:
        cmds (list): 要执行的命令行参数列表
        query (Optional[str]): 可选的查询字符串

    Raises:
        RuntimeError: 当外部命令执行失败时抛出
    """
    try:
        with open("out.log", "w") as fout, open("err.log", "w") as ferr:
            process = Popen(cmds, stdout=fout, stderr=ferr)
            retcode = process.wait()
            if retcode:
                raise RuntimeError("Run Ext Cmd Failed")
    except Exception as e:
        logging.error(f"Error running command: {e}")
        raise


def cif_to_pdb(cif, pdb):
    # clean LIG_* in cif file genarated by AF3 on SMILES
    os.system(f"sed -i 's/LIG_[A-Z]/LIG  /g' {cif}")

    parser = MMCIFParser()
    structure = parser.get_structure("tmp", cif)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb)


def data_to_json(
    protein_str,
    dna_str,
    rna_str,
    ligand_str,
    recycle,
    ensemble,
):
    """将输入数据转换为JSON格式列表

    Args:
        protein_str: 蛋白质序列FASTA字符串
        dna_str: DNA序列FASTA字符串
        rna_str: RNA序列FASTA字符串
        ligand_str: 配体信息字符串
        recycle (int): 循环次数，默认10
        ensemble (int): 集成预测次数，默认1

    Returns:
        list: 包含任务数据的JSON格式列表
    """
    json_data = {}
    json_data["job_name"] = "complex"
    json_data["recycle"] = recycle
    json_data["ensemble"] = ensemble
    entities = []

    def parse_fasta_str(fasta_str):
        if not fasta_str:
            return []
        seqs = []
        sid, seq = None, ""
        for line in fasta_str.split("\n"):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if len(seq) > 0:
                    seqs.append((sid, seq))
                    seq = ""
                sid = line[1:]
            else:
                seq += line.upper()
        if len(seq) > 0:
            seqs.append((sid, seq))
        return seqs

    # 处理蛋白质序列
    if protein_str:
        pro_seqs_batch = parse_fasta_str(protein_str)
        for i, pro_seq in enumerate(pro_seqs_batch):
            entity = {"type": "protein", "sequence": pro_seq[1], "count": 1}
            entities.append(entity)

    # 处理DNA序列
    if dna_str:
        dna_seqs_batch = parse_fasta_str(dna_str)
        for i, dna_seq in enumerate(dna_seqs_batch):
            entity = {"type": "dna", "sequence": dna_seq[1], "count": 1}
            entities.append(entity)

    # 处理RNA序列
    if rna_str:
        rna_seqs_batch = parse_fasta_str(rna_str)
        for i, rna_seq in enumerate(rna_seqs_batch):
            entity = {"type": "rna", "sequence": rna_seq[1], "count": 1}
            entities.append(entity)

    # 处理配体和离子
    if ligand_str:
        for line in ligand_str.strip().split("\n"):
            line = line.strip()
            if not line:
                continue
            if line.startswith("CCD,") or line.startswith("ccd,"):  # CCD codes
                parts = line.split(",")
                for part in parts[1:]:
                    part = part.strip()
                    if len(part) == 2:  # 离子通常是2字符代码
                        entity = {"type": "ion", "ccd": part, "count": 1}
                    else:
                        entity = {"type": "ligand", "ccd": part, "count": 1}
                    entities.append(entity)
            else:  # SMILES
                entity = {"type": "ligand", "smiles": line, "count": 1}
                entities.append(entity)

    json_data["entities"] = entities
    return [json_data]  # 返回列表以符合多任务格式


def run_hf3(json_data):
    """运行HelixFold3预测
    Args:
        json_data: 输入的JSON数据列表或单个任务数据
    Returns:
        bool: 是否运行成功
    """
    try:
        from paddlehelix.task import helixfold3

        # 确保json_data是列表格式
        if not isinstance(json_data, list):
            json_data = [json_data]

        for task in json_data:
            helixfold3.execute(data=task, output_dir="./")
        return True
    except Exception as e:
        logging.error(f"Failed to run HelixFold3: {str(e)}")
        return False


def get_results_single_run(topN: int):
    """获取并保存预测结果中置信度最高的前N个结构

    Args:
        topN (int): 需要保存的最佳结构数量
    """
    try:
        results_dir = "./"
        all_ensemble_results = []

        # 从input.json读取ensemble值和job_name
        with open("input.json") as f:
            input_data = json.load(f)
            job_name = input_data[0].get("job_name", "complex")
            ensemble = input_data[0].get("ensemble", 1)

        # 检查结果目录
        result_base_dir = f"{results_dir}/helixfold3_result_to_download_{job_name}"
        if not os.path.exists(result_base_dir):
            logging.error(f"Result directory does not exist: {result_base_dir}")
            raise RuntimeError(f"预测结果目录不存在: {result_base_dir}")

        # 查找所有ensemble结果目录
        ensemble_dirs = glob.glob(os.path.join(result_base_dir, f"*-*-rank1"))
        if not ensemble_dirs:
            logging.error(f"No ensemble result directories found in: {result_base_dir}")
            return

        # 收集所有ensemble结果
        for ensemble_dir in ensemble_dirs:
            all_results_file = os.path.join(ensemble_dir, "all_results.json")
            if not os.path.exists(all_results_file):
                logging.error(f"Results file does not exist: {all_results_file}")
                continue

            with open(all_results_file) as f:
                result = json.load(f)
                result["source_dir"] = ensemble_dir
                all_ensemble_results.append(result)

        # 根据ranking_confidence排序并获取topN结果
        sorted_results = sorted(
            all_ensemble_results, key=lambda x: x["ranking_confidence"], reverse=True
        )[:topN]

        # 保存topN结果
        for i, result in enumerate(sorted_results):
            output_dir = f"./result_rank{i+1}"
            os.makedirs(output_dir, exist_ok=True)

            source_dir = result["source_dir"]

            # 复制结构文件
            shutil.copy2(
                os.path.join(source_dir, "predicted_structure.cif"),
                os.path.join(output_dir, f"{job_name}.cif"),
            )

            # 复制all_results.json
            shutil.copy2(
                os.path.join(source_dir, "all_results.json"),
                os.path.join(output_dir, "all_results.json"),
            )

            # 复制chain_id_to_input_mapping.json
            mapping_file = os.path.join(source_dir, "chain_id_to_input_mapping.json")
            if os.path.exists(mapping_file):
                shutil.copy2(
                    mapping_file,
                    os.path.join(output_dir, "chain_id_to_input_mapping.json"),
                )

        # 复制user_input.json到结果目录
        if os.path.exists("input.json"):
            shutil.copy2("input.json", os.path.join("./", "user_input.json"))

        logging.info(f"Successfully saved top {topN} results")

    except Exception as e:
        logging.error(f"Failed to process results: {str(e)}")
        raise


ERR_CODE = 100


def main() -> int:
    """主函数

    Returns:
        int: 0表示成功，100表示失败
    """
    try:
        load_dotenv()
        set_logging_default_config()
        logging.info("Starting HelixFold3 prediction pipeline...")

        script_path = __file__
        script_dir = os.path.dirname(script_path)

        try:
            set_progress_value(5)
        except Exception as e:
            logging.error(f"Failed to set progress value: {e}")
            return ERR_CODE

        # 获取参数
        try:
            config_file = os.path.join(script_dir, "..", "..", "cli_config.toml")
            arguments = get_cli_argument(config_file)
            logging.debug(f"Input Arguments: {arguments}")
        except Exception as e:
            logging.error(f"Failed to get CLI arguments: {e}")
            return ERR_CODE

        # 生成JSON
        try:
            json_data = data_to_json(
                arguments.get("protein"),
                arguments.get("dna"),
                arguments.get("rna"),
                arguments.get("ligand"),
                arguments.get("recycle", 10),
                arguments.get("ensemble", 1),
            )
            if not json_data:
                logging.error("Failed to generate input JSON: empty data")
                return ERR_CODE
        except Exception as e:
            logging.error(f"Failed to generate input JSON: {e}")
            return ERR_CODE

        try:
            set_progress_value(20)
        except Exception as e:
            logging.error(f"Failed to set progress value: {e}")
            return ERR_CODE

        # 保存JSON
        try:
            with open("input.json", "w") as fout:
                json.dump(json_data, fout, indent=2)
            set_progress_value(46)
        except Exception as e:
            logging.error(f"Failed to save input JSON: {e}")
            return ERR_CODE

        # 运行预测
        try:
            if not run_hf3(json_data):
                logging.error("Failed to run HelixFold3")
                return ERR_CODE
            set_progress_value(60)
        except Exception as e:
            logging.error(f"Failed to run HelixFold3: {e}")
            return ERR_CODE

        # 处理结果
        try:
            get_results_single_run(arguments.get("topN", 5))
            set_progress_value(89)
        except RuntimeError as e:
            logging.error(f"Failed to process results: {e}")
            return ERR_CODE
        except Exception as e:
            logging.error(f"Unexpected error while processing results: {e}")
            return ERR_CODE

        try:
            set_progress_value(100)
        except Exception as e:
            logging.error(f"Failed to set final progress value: {e}")
            return ERR_CODE

        logging.info("Prediction completed successfully!")
        return 0

    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        return ERR_CODE
