import glob
import logging
import os
import json
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
import zipfile
import time

from dotenv import load_dotenv

from py_script_template.cli import get_cli_argument
from .utils import (
    set_logging_default_config,
    set_progress_value,
)


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


def run_hf3(json_data) -> bool:
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


def get_results_single_run():
    try:
        all_ensemble_results = []

        # 从input.json读取job_name
        with open("input.json") as f:
            input_data = json.load(f)
            job_name = input_data[0].get("job_name", "complex")

        # 查找所有ensemble结果文件夹
        ensemble_dirs = glob.glob(f"{job_name}-*-rank1")
        if not ensemble_dirs:
            raise RuntimeError("未找到预测结果文件夹")

        output_zip = f"helixfold3_result_to_download_{job_name}_{int(time.time())}.zip"

        with zipfile.ZipFile(output_zip, "w") as zf:
            # 这两个文件的处理可以合并到一个循环中
            for input_file in ["input.json", "terms_of_use.md"]:
                if os.path.exists(input_file):
                    dest_name = (
                        "user_input.json" if input_file == "input.json" else input_file
                    )
                    zf.write(input_file, dest_name)

            required_files = [
                "all_results.json",
                "chain_id_to_input_mapping.json",
                "predicted_structure.cif",
            ]

            # 处理每个ensemble结果文件夹
            for ensemble_dir in ensemble_dirs:
                # 检查必需文件并添加到zip
                for file in required_files:
                    file_path = os.path.join(ensemble_dir, file)
                    if not os.path.exists(file_path):
                        raise RuntimeError(f"在{ensemble_dir}中未找到必需文件: {file}")
                    zf.write(file_path, os.path.join(ensemble_dir, file))

                # 读取评分信息
                with open(os.path.join(ensemble_dir, "all_results.json")) as f:
                    results = json.load(f)
                    all_ensemble_results.append(
                        {
                            "ensemble_dir": ensemble_dir,
                            "ptm": results.get("ptm", 0),
                            "iptm": results.get("iptm", 0),
                            "mean_plddt": results.get("mean_plddt", 0),
                            "ranking_confidence": results.get("ranking_confidence", 0),
                        }
                    )

        # 找出最佳结果并输出信息
        best_result = max(all_ensemble_results, key=lambda x: x["ranking_confidence"])

        # 使用f-string合并多行日志输出
        logging.info(
            f"最佳预测结果在 {best_result['ensemble_dir']}\n"
            f"PTM: {best_result['ptm']:.3f}\n"
            f"iPTM: {best_result['iptm']:.3f}\n"
            f"Mean pLDDT: {best_result['mean_plddt']:.3f}\n"
            f"Ranking confidence: {best_result['ranking_confidence']:.3f}"
        )

        # 转换最佳结构格式
        best_cif = os.path.join(best_result["ensemble_dir"], "predicted_structure.cif")
        best_pdb = os.path.join(best_result["ensemble_dir"], "predicted_structure.pdb")
        cif_to_pdb(best_cif, best_pdb)

        logging.info(f"结果已保存到: {output_zip}")
        return True

    except Exception as e:
        logging.error(f"处理结果时发生错误: {str(e)}")
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
            get_results_single_run()
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
