import os
import json
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
import zipfile
import shutil
import csv
import sys

from dotenv import load_dotenv

from py_script_template.cli import get_cli_argument
from .utils import (
    set_logging_default_config,
    set_progress_value,
    load_fasta,
    load_ligands,
)


# 定义支持的氨基酸修饰
VALID_MODIFICATIONS = {
    "R": ["2MR", "AGM", "CIR"],
    "C": ["MCS", "P1L", "SNC"],
    "H": ["NEP", "HIP"],
    "K": ["ALY", "MLY", "M3L", "MLZ", "LYZ", "KCR", "YHA"],
    "N": ["AHB", "SNN"],
    "P": ["HYP", "HY3"],
    "S": ["SEP"],
    "T": ["TPO"],
    "W": ["TRF"],
    "Y": ["PTR"],
}

# 添加DNA修饰的定义
VALID_DNA_MODIFICATIONS = {
    "A": ["6MA", "3DR"],
    "C": ["5CM", "C34", "5HC", "1CC", "5FC", "3DR"],
    "T": ["3DR"],
    "G": ["6OG", "8OG", "3DR"],
}

# 添加RNA修饰的定义
VALID_RNA_MODIFICATIONS = {
    "A": ["A2M", "MA6", "6MZ"],
    "C": ["5MC", "OMC", "4OC", "RSQ"],
    "U": ["5MU", "OMU", "UR3", "PSU"],
    "G": ["2MG", "OMG", "7MG"],
}

# 添加支持的离子定义
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


def detect_sequence_type(sequence):
    """根据序列内容判断序列类型"""
    sequence = sequence.upper()
    protein_chars = set("ACDEFGHIKLMNPQRSTVWY")
    dna_chars = set("ATCG")
    rna_chars = set("AUCG")

    seq_chars = set(sequence)

    # 如果序列只包含AUCG，则可能是RNA
    if seq_chars.issubset(rna_chars):
        return "rna"
    # 如果序列只包含ATCG，则可能是DNA
    elif seq_chars.issubset(dna_chars):
        return "dna"
    # 如果序列包含蛋白质特有的氨基酸，则是蛋白质
    elif seq_chars.issubset(protein_chars):
        return "protein"
    else:
        raise ValueError(
            f"Unable to determine sequence type. Invalid characters: {seq_chars - (protein_chars | dna_chars | rna_chars)}"
        )


def parse_ion_string(ion_str):
    """解析离子输入字符串

    Args:
        ion_str: 形如 "MG:2,ZN,CU:3" 的字符串

    Returns:
        list: 包含(ion_ccd, count)元组的列表

    Examples:
        # 以下输入都是有效的：
        ion = "MG:2,ZN,CU:3"  # 2个MG离子，1个ZN离子，3个CU子
        ion = "MG,MG,ZN"      # 2个MG离子，1个ZN离子
        ion = "MG:20"         # 20个MG离子
    """
    if not ion_str:
        return []

    result = []
    # 分割不同的离子
    ion_parts = ion_str.strip().split(",")

    # 用于统计每子的出现次数
    ion_counts = {}

    for part in ion_parts:
        # 处理带有计数的情况 (如 "MG:2")
        if ":" in part:
            ccd, count = part.split(":")
            ccd = ccd.strip().upper()
            try:
                count = int(count)
            except ValueError:
                raise ValueError(f"离子计数必须是整数: {part}")
        else:
            # 不带计数的情况，默认为1
            ccd = part.strip().upper()
            # 如果这个离子已经出现过，增加计数
            if ccd in ion_counts:
                ion_counts[ccd] = ion_counts.get(ccd, 0) + 1
                continue
            count = 1

        if count < 1:
            raise ValueError(f"离子计数必须大于0: {part}")

        # 更新离子计数
        ion_counts[ccd] = ion_counts.get(ccd, 0) + count

    # 将累积的计数转换为结果
    for ccd, total_count in ion_counts.items():
        result.append((ccd, total_count))

    return result


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


def validate_protein_sequence(sequence, sid):
    """验证蛋白质序列"""
    if not sequence:
        raise ValueError(f"Empty sequence for {sid}")

    if len(sequence) > 2000:
        raise ValueError(f"Protein sequence length exceeds 2000: {sid}")

    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_aa = set(sequence) - valid_aa
    if invalid_aa:
        raise ValueError(
            f"Invalid amino acids found in protein sequence {sid}: {invalid_aa}"
        )

    return sequence


def validate_count(count, entity_type, max_count=50):
    """验证实体的count值"""
    if not isinstance(count, int):
        raise ValueError(f"{entity_type} count must be an integer, got {type(count)}")
    if count < 1:
        raise ValueError(f"{entity_type} count must be at least 1, got {count}")
    if count > max_count:
        raise ValueError(f"{entity_type} count cannot exceed {max_count}, got {count}")
    return count


# 在create_protein_entity中确保修饰被正确添加
def create_protein_entity(sequence, sid, modifications=None, count=1):
    """创建蛋白质实体"""
    print("=" * 80, file=sys.stderr)
    print(f"处理序列修饰: {modifications}", file=sys.stderr)
    print("=" * 80, file=sys.stderr)
    
    # 添加这段代码来显示完整序列
    print(f"完整序列: {sequence}", file=sys.stderr)
    print("序列位置对照:", file=sys.stderr)
    for i, aa in enumerate(sequence, 1):
        print(f"{i}: {aa}", end=' ', file=sys.stderr)
    print("\n", file=sys.stderr)
    
    if modifications:
        print("序列中的氨基酸位置:", file=sys.stderr)
        for mod in modifications:
            index = mod["index"]
            aa = sequence[index - 1]
            print(f"位置 {index}: {aa}", file=sys.stderr)
    
    sequence = validate_protein_sequence(sequence, sid)
    count = validate_count(count, "Protein")

    entity = {"type": "protein", "sequence": sequence, "count": count}

    if modifications:
        try:
            print("*" * 80, modifications)
            # 验证修饰信息
            for mod in modifications:
                # if mod["type"] != "residue_replace":
                #     raise ValueError(f"不支持的修饰类型: {mod['type']}")

                index = mod["index"]
                if not isinstance(index, int):
                    raise ValueError(f"修饰位置必须是整数，得到: {type(index)}")
                if index < 1 or index > len(sequence):
                    raise ValueError(
                        f"修饰位置 {index} 超出序列范围 (1-{len(sequence)})"
                    )

                aa = sequence[index - 1]
                ccd = mod["ccd"]

                if aa not in VALID_MODIFICATIONS:
                    raise ValueError(f"氨基酸 {aa} 在位置 {index} 不支持修饰")
                if ccd not in VALID_MODIFICATIONS[aa]:
                    raise ValueError(
                        f"氨基酸 {aa} 不支持修饰类型 {ccd}。"
                        f"支持的修饰类型: {', '.join(VALID_MODIFICATIONS[aa])}"
                    )

            # 添加修饰信息到实体
            entity["modification"] = modifications

        except ValueError as e:
            raise ValueError(f"蛋白质 {sid} 的修饰验证失败: {str(e)}")

    return entity


def validate_dna_sequence(sequence, sid):
    """验证DNA序列"""
    if not sequence:
        raise ValueError(f"Empty sequence for {sid}")

    if len(sequence) > 2000:
        raise ValueError(f"DNA sequence length exceeds 2000: {sid}")

    sequence = sequence.upper()
    valid_bases = set("ATCG")
    invalid_bases = set(sequence) - valid_bases
    if invalid_bases:
        raise ValueError(f"Invalid bases found in DNA sequence {sid}: {invalid_bases}")

    return sequence


def create_dna_entity(sequence, sid, modifications=None, count=1):
    """创建DNA实体"""
    sequence = validate_dna_sequence(sequence, sid)
    count = validate_count(count, "DNA")

    entity = {"type": "dna", "sequence": sequence, "count": count}

    if modifications:
        try:
            # 验证修饰信息
            for mod in modifications:
                if mod["type"] != "residue_replace":
                    raise ValueError(f"不支持的修饰类型: {mod['type']}")

                index = mod["index"]
                if not isinstance(index, int):
                    raise ValueError(f"修饰位置必须是整数，得到: {type(index)}")
                if index < 1 or index > len(sequence):
                    raise ValueError(
                        f"修饰位置 {index} 超出序列范围 (1-{len(sequence)})"
                    )

                base = sequence[index - 1]
                ccd = mod["ccd"]

                if base not in VALID_DNA_MODIFICATIONS:
                    raise ValueError(f"核苷酸 {base} 在位置 {index} 不支持修饰")
                if ccd not in VALID_DNA_MODIFICATIONS[base]:
                    raise ValueError(
                        f"核苷酸 {base} 不支持修饰类型 {ccd}。"
                        f"支持的修饰类型: {', '.join(VALID_DNA_MODIFICATIONS[base])}"
                    )

            # 添加修饰信息到实体
            entity["modification"] = [
                {
                    "type": "residue_replace",
                    "index": mod["index"],
                    "ccd": mod["ccd"],
                }
                for mod in modifications
            ]
        except ValueError as e:
            raise ValueError(f"DNA {sid} 的修饰验证失败: {str(e)}")

    return entity


def validate_rna_sequence(sequence, sid):
    """验证RNA序列"""
    if not sequence:
        raise ValueError(f"Empty sequence for {sid}")

    if len(sequence) > 2000:
        raise ValueError(f"RNA sequence length exceeds 2000: {sid}")

    sequence = sequence.upper().replace("T", "U")
    valid_bases = set("AUCG")
    invalid_bases = set(sequence) - valid_bases
    if invalid_bases:
        raise ValueError(f"Invalid bases found in RNA sequence {sid}: {invalid_bases}")

    return sequence


def create_rna_entity(sequence, sid, modifications=None, count=1):
    """创建RNA实体"""
    sequence = validate_rna_sequence(sequence, sid)
    count = validate_count(count, "RNA")

    entity = {"type": "rna", "sequence": sequence, "count": count}

    if modifications:
        try:
            # 验证修饰信息
            for mod in modifications:
                if mod["type"] != "residue_replace":
                    raise ValueError(f"不支持的修饰类型: {mod['type']}")

                index = mod["index"]
                if not isinstance(index, int):
                    raise ValueError(f"修饰位置必须是整数，得到: {type(index)}")
                if index < 1 or index > len(sequence):
                    raise ValueError(
                        f"修饰位置 {index} 超出序列范围 (1-{len(sequence)})"
                    )

                base = sequence[index - 1]
                ccd = mod["ccd"]

                if base not in VALID_RNA_MODIFICATIONS:
                    raise ValueError(f"核苷酸 {base} 在位置 {index} 不支持修饰")
                if ccd not in VALID_RNA_MODIFICATIONS[base]:
                    raise ValueError(
                        f"核苷酸 {base} 不支持修饰类型 {ccd}。"
                        f"支持的修饰类型: {', '.join(VALID_RNA_MODIFICATIONS[base])}"
                    )

            entity["modification"] = modifications
        except ValueError as e:
            raise ValueError(f"RNA {sid} 的修饰验证失败: {str(e)}")

    return entity


def create_ligand_entity(ligand, count=1):
    """创建配体实体"""
    count = validate_count(count, "Ligand")

    try:
        if isinstance(ligand, list):  # CCD codes
            if len(ligand) < 2:
                raise ValueError("Empty CCD code list")
            entity = {"type": "ligand", "ccd": ligand[1], "count": count}
        else:  # SMILES
            if not isinstance(ligand, str):
                raise ValueError(f"SMILES must be a string, got {type(ligand)}")
            if not ligand:
                raise ValueError("Empty SMILES string")
            # TODO: 添加对SMILES重核数量验证
            entity = {"type": "ligand", "smiles": ligand, "count": count}
        return entity
    except ValueError as e:
        raise ValueError(f"Invalid ligand: {str(e)}")


def create_ion_entity(ion_ccd, count=1):
    """创建离子实体"""
    count = validate_count(count, "Ion")

    try:
        if not isinstance(ion_ccd, str):
            raise ValueError(f"Ion CCD code must be a string, got {type(ion_ccd)}")
        if ion_ccd not in VALID_IONS:
            raise ValueError(
                f"Invalid ion CCD code: {ion_ccd}. "
                f"Valid options are: {', '.join(sorted(VALID_IONS))}"
            )
        return {"type": "ion", "ccd": ion_ccd, "count": count}
    except ValueError as e:
        raise ValueError(f"Invalid ion: {str(e)}")


def parse_modifications_file(file_path):
    print(f"开始解析修饰文件: {file_path}", file=sys.stderr)
    if not file_path or not os.path.exists(file_path):
        print(f"警告: 修饰文件不存在: {file_path}", file=sys.stderr)
        return {}

    modifications_dict = {}

    try:
        with open(file_path) as f:
            for line in f.readlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                try:
                    # 解析每行的三个字段
                    tmp = line.split(",")
                    if len(tmp) != 3:
                        raise ValueError(
                            f"每行必须包含三个字段 (序列编号,CCD编号,位置编号), 得到: {line}"
                        )

                    seq_index = int(tmp[0])
                    type_mdf = tmp[1].strip()
                    pos_aa = int(tmp[2])

                    # 验证格式
                    # if not (seq_index >= 1 and pos_aa >= 1 and len(type_mdf) == 3):
                    #     raise ValueError(f"Error format in modifications: {line}")

                    # 将修饰信息添加到对应序列列表中
                    if seq_index not in modifications_dict:
                        modifications_dict[seq_index] = []

                    # 创建修饰对象
                    mod_obj = {
                        "type": "residue_replace",
                        "index": pos_aa,
                        "ccd": type_mdf,
                    }
                    modifications_dict[seq_index].append(mod_obj)

                except ValueError as e:
                    raise ValueError(f"Error parsing modification: {str(e)}")

    except Exception as e:
        raise ValueError(f"读取修饰文件时出错: {str(e)}")

    return modifications_dict


def data_to_json(
    protein_file,
    dna_file,
    rna_file,
    ligand_file,
    ion,
    recycle,
    ensemble,
    modifications,
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

    # 解析修饰文件
    print(f"准备处理修饰文件: {modifications}", file=sys.stderr)
    
    modifications_dict = {}
    if modifications:
        try:
            print(f"正在从文件读取修饰信息: {modifications}", file=sys.stderr)
            modifications_dict = parse_modifications_file(modifications)
            print(f"解析到的修饰信息: {modifications_dict}", file=sys.stderr)
        except ValueError as e:
            print(f"错误: 处理修饰文件时出错: {str(e)}", file=sys.stderr)
            raise

    # 用于追踪当前处理的序列编号
    current_seq_num = 0  # 从0开始，在处理序列前递增

    # 处理所有序列文件
    for file_path, file_type in [
        (protein_file, "protein"),
        (dna_file, "dna"),
        (rna_file, "rna"),
    ]:
        if file_path:
            seqs = load_fasta(file_path)
            for sid, sequence in seqs:
                current_seq_num += 1  # 在处理序列前递增

                # 获取该序列的修饰信息
                seq_modifications = modifications_dict.get(current_seq_num, [])
                print(
                    f"Processing sequence {current_seq_num} with modifications: {seq_modifications}",
                    file=sys.stderr
                )  # 临时调试

                # 根据输入参数类型处理序列
                if file_type == "protein":
                    try:
                        entity = create_protein_entity(
                            sequence, sid, modifications=seq_modifications
                        )
                        if seq_modifications:  # 确保有修饰信息时才添加
                            entity["modification"] = seq_modifications
                        json_data["entities"].append(entity)
                    except ValueError as e:
                        print(f"Warning: Skipping sequence {sid}: {str(e)}", file=sys.stderr)
                        continue

                elif file_type == "dna":
                    try:
                        entity = create_dna_entity(
                            sequence, sid, modifications=seq_modifications
                        )
                        if seq_modifications:  # 确保有修饰信息时才添加
                            entity["modification"] = seq_modifications
                        json_data["entities"].append(entity)
                    except ValueError as e:
                        print(f"Warning: Skipping sequence {sid}: {str(e)}", file=sys.stderr)
                        continue

                elif file_type == "rna":
                    try:
                        entity = create_rna_entity(
                            sequence, sid, modifications=seq_modifications
                        )
                        json_data["entities"].append(entity)
                    except ValueError as e:
                        print(f"Warning: Skipping sequence {sid}: {str(e)}", file=sys.stderr)
                        continue

    # 处理配体
    if ligand_file:
        ligands = load_ligands(ligand_file)
        for ligand in ligands:
            entity = create_ligand_entity(ligand)
            json_data["entities"].append(entity)

    # 处理离子
    if ion:
        try:
            ion_list = parse_ion_string(ion)
            for ion_ccd, count in ion_list:
                entity = create_ion_entity(ion_ccd, count)
                json_data["entities"].append(entity)
        except ValueError as e:
            raise ValueError(f"处理离子输入时出错: {str(e)}")

    # 验证实体数量和总长度
    if not json_data["entities"]:
        raise ValueError("No valid entities found in input files")

    # 计算总token数量
    total_tokens = 0
    for entity in json_data["entities"]:
        if entity["type"] in ["protein", "dna", "rna"]:
            total_tokens += len(entity["sequence"])
        elif entity["type"] == "ligand" and "smiles" in entity:
            # 配体中的一个子算做一个token
            # TODO: 实现更准确的SMILES token计算
            total_tokens += len(entity["smiles"])  # 这是一个简化的计算方式

    if total_tokens > 2000:
        raise ValueError(f"Total sequence length exceeds 2000: {total_tokens}")

    # 在返回之前打印生成的JSON，方便调试
    print(f"Generated JSON data: {json.dumps(json_data, indent=2)}", file=sys.stderr)

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
        print(f"Error: Failed to run HelixFold3: {str(e)}", file=sys.stderr)
        return False


def get_results_single_run():
    try:
        # 从input.json读取job_name和ensemble数量
        with open("input.json") as f:
            input_data = json.load(f)
            job_name = input_data[0].get("job_name", "complex")

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
                    raise RuntimeError(f"在{ensemble_dir}中未找到需文件: {file}")

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
            fieldnames = ["Name", "Ranking_Score"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for result in all_ensemble_results:
                writer.writerow(
                    {
                        "Name": f"rank{result['rank']}",
                        "Ranking_Score": f"{result['ranking_confidence']:.3f}",
                    }
                )

        # 复制并重命名所有结构件，按rank排序
        for result in all_ensemble_results:
            src_cif = os.path.join(result["ensemble_dir"], "predicted_structure.cif")
            rank_num = result["rank"]

            # 复制并重命名CIF文件
            dst_cif = f"./rank_{rank_num}.cif"
            shutil.copy2(src_cif, dst_cif)

        # 出结果信息
        print(
            f"\n最佳结果(rank{best_result['rank']}):\n"
            f"PTM: {best_result['ptm']:.3f}\n"
            f"iPTM: {best_result['iptm']:.3f}\n"
            f"Mean pLDDT: {best_result['mean_plddt']:.3f}\n"
            f"Ranking confidence: {best_result['ranking_confidence']:.3f}",
            file=sys.stderr
        )

        # 输出所有ensemble结果的要信息
        print("\n所有预测结果排名:", file=sys.stderr)
        for result in all_ensemble_results:
            print(
                f"rank{result['rank']}: "
                f"Ranking confidence = {result['ranking_confidence']:.3f}",
                file=sys.stderr
            )

        # 更新输出文件信息
        print(
            f"\n结果文件已保存:\n"
            f"- rank_1.cif 至 rank_{len(all_ensemble_results)}.cif (所有结构预测结果)\n"
            f"- ranking_scores.csv (排名信息)",
            file=sys.stderr
        )

        print(f"原结果文件位置: {result_zip}", file=sys.stderr)
        return True

    except Exception as e:
        print(
            f"处理结果时发生错误: {str(e)} | 或计算失败意外导致结果异常, 请确认参数后重试",
            file=sys.stderr
        )
        raise


ERR_CODE = 100


def main() -> int:
    """Main function
    Returns:
        int: 0 for success, 100 for failure
    """
    try:
        load_dotenv()
        print("Starting HelixFold3 prediction pipeline...", file=sys.stderr)

        script_path = __file__
        script_dir = os.path.dirname(script_path)

        try:
            set_progress_value(5)
        except Exception as e:
            print(f"Error: Failed to set progress value: {e}", file=sys.stderr)
            return ERR_CODE

        # 获取参数
        try:
            config_file = os.path.join(script_dir, "..", "..", "cli_config.toml")
            arguments = get_cli_argument(config_file)
            print(f"Input Arguments: {arguments}", file=sys.stderr)
        except Exception as e:
            print(f"Error: Failed to get CLI arguments: {e}", file=sys.stderr)
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
                arguments.get("modification"),
                arguments.get("job_name", "complex"),
            )
            if not json_data:
                print("Error: Failed to generate input JSON: empty data", file=sys.stderr)
                return ERR_CODE
        except Exception as e:
            print(f"Error: Failed to generate input JSON: {e}", file=sys.stderr)
            return ERR_CODE

        try:
            set_progress_value(20)
        except Exception as e:
            print(f"Error: Failed to set progress value: {e}", file=sys.stderr)
            return ERR_CODE

        # 保存JSON
        try:
            with open("input.json", "w") as fout:
                json.dump(json_data, fout, indent=2)
            set_progress_value(46)
        except Exception as e:
            print(f"Error: Failed to save input JSON: {e}", file=sys.stderr)
            return ERR_CODE

        # 运行预
        try:
            if not run_hf3(json_data):
                print("Error: Failed to run HelixFold3", file=sys.stderr)
                return ERR_CODE
            set_progress_value(60)
        except Exception as e:
            print(f"Error: Failed to run HelixFold3: {e}", file=sys.stderr)
            return ERR_CODE

        # 处理结果
        try:
            get_results_single_run()
            set_progress_value(89)
        except RuntimeError as e:
            print(f"Error: Failed to process results: {e}", file=sys.stderr)
            return ERR_CODE
        except Exception as e:
            print(f"Error: Unexpected error while processing results: {e}", file=sys.stderr)
            return ERR_CODE

        try:
            set_progress_value(100)
        except Exception as e:
            print(f"Error: Failed to set final progress value: {e}", file=sys.stderr)
            return ERR_CODE

        print("Prediction completed successfully!", file=sys.stderr)
        return 0

    except Exception as e:
        print(f"Error: An unexpected error occurred: {e}", file=sys.stderr)
        return ERR_CODE
