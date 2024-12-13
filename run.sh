#!/bin/bash

# 设置执行路径
SCRIPT_PATH="/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py"

# 小型蛋白质测试
# echo "Running small protein test..."
# $SCRIPT_PATH --job_name pro_300AA --protein ./test_files/protein/pro_300AA.fasta

# 中型蛋白质测试
# echo "Running medium protein test..."
# $SCRIPT_PATH --job_name pro_600AA --protein ./test_files/protein/pro_600AA.fasta

# 蛋白质-DNA复合物测试
# echo "Running protein-DNA complex test..."
# $SCRIPT_PATH --job_name pro_DNA \
#     --protein ./test_files/protein-DNA/pro_DNA.fasta \
#     --dna ./test_files/protein-DNA/pro_DNA.fasta

# 蛋白质-RNA复合物测试
# echo "Running protein-RNA complex test..."
# $SCRIPT_PATH --job_name pro_RNA \
#     --protein ./test_files/protein-RNA/pro_RNA.fasta \
#     --rna ./test_files/protein-RNA/pro_RNA.fasta

# 蛋白质-配体复合物测试
# echo "Running protein-ligand complex test..."
# $SCRIPT_PATH --job_name pro_lig \
#     --protein ./test_files/protein-ligand/pro_lig.fasta \
#     --ligand ./test_files/protein-ligand/ligand.txt

# 蛋白质-配体-离子复合物测试
echo "Running protein-ligand-ion complex test..."
$SCRIPT_PATH --job_name pro_lig_ion \
    --protein ./test_files/potein-ligand-ion/pro_lig_ion.fasta \
    --ligand ./test_files/potein-ligand-ion/ligand.txt \
    --ion "ZN"

# 蛋白质-RNA-配体复合物测试
# echo "Running protein-RNA-ligand complex test..."
# $SCRIPT_PATH --job_name pro_RNA_lig \
#     --protein ./test_files/protein-RNA-ligand/pro_RNA_lig.fasta \
#     --rna ./test_files/protein-RNA-ligand/pro_RNA_lig.fasta \
#     --ligand ./test_files/protein-RNA-ligand/ligand.txt