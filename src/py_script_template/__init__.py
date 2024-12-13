import logging
import os
import json
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
import zipfile
import shutil
import csv

from dotenv import load_dotenv

from py_script_template.cli import get_cli_argument
from .utils import (
    set_logging_default_config,
    set_progress_value,
    load_fasta,
    load_ligands,
)


def cif_to_pdb(cif, pdb):
    # 使用Python读取和处理文件，替代sed命令
    with open(cif, "r") as f:
        content = f.read()

    # 替换LIG_[A-Z]为LIG
    import re

    content = re.sub(r"LIG_[A-Z]", "LIG  ", content)

    with open(cif, "w") as f:
        f.write(content)

    parser = MMCIFParser()
    structure = parser.get_structure("tmp", cif)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb)


def data_to_json(
    protein_file,
    dna_file,
    rna_file,
    ligand_file,
    ion,
    recycle,
    ensemble,
    job_name="complex",
):
    """Convert input data to JSON format list

    Args:
        protein_file: Path to protein sequence FASTA file
        dna_file: Path to DNA sequence FASTA file
        rna_file: Path to RNA sequence FASTA file
        ligand_file: Path to ligand entries file
        ion: Ion CCD code string
        recycle (int): Number of recycles, default 10
        ensemble (int): Number of ensemble predictions, default 1
        job_name (str): Name of the job, default "complex"

    Returns:
        list: List of task data in JSON format
    """
    # 验证基本参数
    valid_recycle_values = {10, 20, 50, 100}
    if recycle not in valid_recycle_values:
        raise ValueError(
            f"recycle must be one of {valid_recycle_values}, got {recycle}"
        )

    valid_ensemble_values = {1, 5, 10, 100}
    if ensemble not in valid_ensemble_values:
        raise ValueError(
            f"ensemble must be one of {valid_ensemble_values}, got {ensemble}"
        )

    if not isinstance(job_name, str) or not job_name:
        raise ValueError("job_name must be a non-empty string")

    json_data = {
        "job_name": job_name,
        "recycle": recycle,
        "ensemble": ensemble,
        "entities": [],
    }

    # 处理蛋白质序列
    if protein_file:
        seqs = load_fasta(protein_file)
        for sid, sequence in seqs:
            if sid.endswith(".A"):  # 只处理蛋白质序列
                if len(sequence) > 2000:
                    raise ValueError(f"Protein sequence length exceeds 2000: {sid}")

                # 检查是否只包含20种标准氨基酸
                valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
                invalid_aa = set(sequence) - valid_aa
                if invalid_aa:
                    raise ValueError(
                        f"Invalid amino acids found in protein sequence {sid}: {invalid_aa}"
                    )

                entity = {"type": "protein", "sequence": sequence, "count": 1}
                json_data["entities"].append(entity)

    # 处理DNA序列
    if dna_file:
        seqs = load_fasta(dna_file)
        for sid, sequence in seqs:
            if sid.endswith((".B", ".C")):  # 只处理DNA序列
                if len(sequence) > 2000:
                    raise ValueError(f"DNA sequence length exceeds 2000: {sid}")

                # 将序列转换为大写进行检查
                sequence = sequence.upper()
                # 检查是否只包含ATCG
                valid_bases = set("ATCG")
                invalid_bases = set(sequence) - valid_bases
                if invalid_bases:
                    raise ValueError(
                        f"Invalid bases found in DNA sequence {sid}: {invalid_bases}"
                    )

                entity = {
                    "type": "dna",
                    "sequence": sequence.upper(),  # DNA序列使用小写
                    "count": 1,
                }
                json_data["entities"].append(entity)

    # 处理RNA序列
    if rna_file:
        seqs = load_fasta(rna_file)
        for sid, sequence in seqs:
            if sid.endswith(".R"):  # 只处理RNA序列
                if len(sequence) > 2000:
                    raise ValueError(f"RNA sequence length exceeds 2000: {sid}")

                # 将序列转换为大写并替换T为U
                sequence = sequence.upper().replace("T", "U")

                # 检查是否只包含AUCG
                valid_bases = set("AUCG")
                invalid_bases = set(sequence) - valid_bases
                if invalid_bases:
                    raise ValueError(
                        f"Invalid bases found in RNA sequence {sid}: {invalid_bases}"
                    )

                entity = {
                    "type": "rna",
                    "sequence": sequence.upper(),  # RNA序列使用大写
                    "count": 1,
                }
                json_data["entities"].append(entity)

    # 处理配体
    if ligand_file:
        ligands = load_ligands(ligand_file)
        for ligand in ligands:
            if isinstance(ligand, list):  # CCD codes
                for part in ligand[1:]:  # 跳过 "CCD," 前缀
                    entity = {"type": "ligand", "ccd": part, "count": 1}
                    json_data["entities"].append(entity)
            else:  # SMILES
                entity = {"type": "ligand", "smiles": ligand, "count": 1}
                json_data["entities"].append(entity)

    # 处理离子
    VALID_IONS = {
        "MG",  # MG2+
        "ZN",  # ZN2+
        "CL",  # CL-
        "CA",  # CA2+
        "NA",  # NA+
        "MN",  # MN2+
        "MN3",  # MN3+
        "K",  # K+
        "FE",  # FE3+
        "FE2",  # FE2+
        "CU",  # CU2+
        "CU1",  # CU1+
        "CU3",  # CU3+
        "CO",  # CO2+
    }

    if ion:  # 如果提供了ion参数
        if ion in VALID_IONS:
            entity = {"type": "ion", "ccd": ion, "count": 1}
            json_data["entities"].append(entity)
        else:
            raise ValueError(
                f"Invalid ion CCD code: {ion}. Must be one of: {', '.join(sorted(VALID_IONS))}"
            )

    # 验证实体数量和总长度
    if not json_data["entities"]:
        raise ValueError("No valid entities found in input files")

    # 计算总token数量
    total_tokens = 0
    for entity in json_data["entities"]:
        if entity["type"] in ["protein", "dna", "rna"]:
            total_tokens += len(entity["sequence"])
        elif entity["type"] == "ligand" and "smiles" in entity:
            # 配体中的一个原子算做一个token
            # TODO: 实现更准确的SMILES token计算
            total_tokens += len(entity["smiles"])  # 这是一个简化的计算方式

    if total_tokens > 2000:
        raise ValueError(f"Total sequence length exceeds 2000: {total_tokens}")

    # 在返回之前打印生成的JSON，方便调试
    logging.debug(f"Generated JSON data: {json.dumps(json_data, indent=2)}")

    return [json_data]  # 始终返回列表格式


def run_hf3(json_data) -> bool:
    """Run HelixFold3 prediction
    Args:
        json_data: Input JSON data list or single task data
        mode: TODO Execution mode, "single" or "batch"
    Returns:
        bool: Whether the run was successful
    """
    try:
        from paddlehelix.task import helixfold3

        # Ensure json_data is in list format
        if not isinstance(json_data, list):
            json_data = [json_data]

        # In single mode, only process the first task
        helixfold3.execute(data=json_data[0], output_dir="./", overwrite=True)
        return True
    except Exception as e:
        logging.error(f"Failed to run HelixFold3: {str(e)}")
        return False


def get_results_single_run():
    try:
        # 从input.json读取job_name和ensemble数量
        with open("input.json") as f:
            input_data = json.load(f)
            job_name = input_data[0].get("job_name", "complex")
            ensemble_size = input_data[0].get("ensemble", 1)

        # 检查数据目录下的结果
        data_dir = os.path.join("data", f"{job_name}_0")
        if not os.path.exists(data_dir):
            raise RuntimeError(f"未找到预测结果文件夹: {data_dir}")

        # 处理和解压结果文件
        result_zip = None
        for file in os.listdir(data_dir):
            if file.startswith("helixfold3_result_to_download_") and file.endswith(
                ".zip"
            ):
                result_zip = os.path.join(data_dir, file)
                break

        if not result_zip:
            raise RuntimeError("未找到结果zip文件")

        output_dir = os.path.join(data_dir, "extracted_results")
        os.makedirs(output_dir, exist_ok=True)

        with zipfile.ZipFile(result_zip, "r") as zf:
            zf.extractall(output_dir)

        # 查找所有ensemble结果文件夹
        ensemble_dirs = []
        for root, dirs, files in os.walk(output_dir):
            for dir_name in dirs:
                if dir_name.endswith("-rank1"):
                    ensemble_dirs.append(os.path.join(root, dir_name))

        if not ensemble_dirs:
            raise RuntimeError(f"在解压后的目录中未找到预测结果文件夹: {output_dir}")

        # 处理每个ensemble结果
        all_ensemble_results = []
        for ensemble_dir in ensemble_dirs:
            # 检查必需文件
            required_files = ["all_results.json", "predicted_structure.cif"]
            for file in required_files:
                file_path = os.path.join(ensemble_dir, file)
                if not os.path.exists(file_path):
                    raise RuntimeError(f"在{ensemble_dir}中未找到必需文件: {file}")

            # 获取ensemble编号
            ensemble_id = os.path.basename(ensemble_dir).split("-")[1]

            with open(os.path.join(ensemble_dir, "all_results.json")) as f:
                results = json.load(f)
                all_ensemble_results.append(
                    {
                        "ensemble_dir": ensemble_dir,
                        "ensemble_id": ensemble_id,
                        "ranking_confidence": results.get("ranking_confidence", 0),
                        "ptm": results.get("ptm", 0),
                        "iptm": results.get("iptm", 0),
                        "mean_plddt": results.get("mean_plddt", 0),
                    }
                )

        # 按ranking_confidence排序结果并分配rank
        all_ensemble_results.sort(key=lambda x: x["ranking_confidence"], reverse=True)
        for rank, result in enumerate(all_ensemble_results, 1):
            result["rank"] = rank

        best_result = all_ensemble_results[0]  # 排序后的第一个结果就是最佳结果

        # 生成简化的CSV文件，只包含rank和ranking_score
        with open("./ranking_scores.csv", "w", newline="") as csvfile:
            fieldnames = ["Rank", "Ranking_Score"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for result in all_ensemble_results:
                writer.writerow(
                    {
                        "Rank": f"rank{result['rank']}",
                        "Ranking_Score": f"{result['ranking_confidence']:.3f}",
                    }
                )

        # 复制并重命名所有结构文件，按rank排序
        for result in all_ensemble_results:
            src_cif = os.path.join(result["ensemble_dir"], "predicted_structure.cif")
            rank_num = result["rank"]

            # 复制并重命名CIF文件
            dst_cif = f"./rank_{rank_num}.cif"
            shutil.copy2(src_cif, dst_cif)

        # 出结果信息
        logging.info(
            f"\n最佳结果 (rank{best_result['rank']}):\n"
            f"PTM: {best_result['ptm']:.3f}\n"
            f"iPTM: {best_result['iptm']:.3f}\n"
            f"Mean pLDDT: {best_result['mean_plddt']:.3f}\n"
            f"Ranking confidence: {best_result['ranking_confidence']:.3f}"
        )

        # 输出所有ensemble结果的要信息
        logging.info("\n所有预测结果排名:")
        for result in all_ensemble_results:
            logging.info(
                f"rank{result['rank']}: "
                f"Ranking confidence = {result['ranking_confidence']:.3f}"
            )

        # 更新输出文件信息
        logging.info(
            f"\n结果文件已保存:\n"
            f"- rank_1.cif 至 rank_{len(all_ensemble_results)}.cif (所有结构预测结果)\n"
            f"- ranking_scores.csv (排名信息)"
        )

        logging.info(f"原始结果文件位于: {result_zip}")
        return True

    except Exception as e:
        logging.error(f"处理结果时发生错误: {str(e)}")
        raise


ERR_CODE = 100


def main() -> int:
    """Main function
    Returns:
        int: 0 for success, 100 for failure
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
                arguments.get("ion"),
                arguments.get("recycle", 10),
                arguments.get("ensemble", 1),
                arguments.get("job_name", "complex"),
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
