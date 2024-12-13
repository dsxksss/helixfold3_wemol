from subprocess import run


def test_protein_small():
    protein_small_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--job_name",
        "pro_300AA",
        "--protein",
        "./test_files/protein/pro_300AA.fasta",
    ]
    assert run(protein_small_run_cmd).returncode == 0


def test_protein_medium():
    protein_medium_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--job_name",
        "pro_600AA",
        "--protein",
        "./test_files/protein/pro_600AA.fasta",
    ]
    assert run(protein_medium_run_cmd).returncode == 0


def test_complex_protein_dna():
    complex_protein_dna_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--job_name",
        "pro_DNA",
        "--protein",
        "./test_files/protein-DNA/pro_DNA.fasta",
        "--dna",
        "./test_files/protein-DNA/pro_DNA.fasta",
    ]
    assert run(complex_protein_dna_run_cmd).returncode == 0


def test_complex_protein_ligand():
    complex_protein_ligand_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--job_name",
        "pro_lig",
        "--protein",
        "./test_files/protein-ligand/pro_lig.fasta",
        "--ligand",
        "./test_files/protein-ligand/ligand.txt",  # 使用固定文件
    ]
    assert run(complex_protein_ligand_run_cmd).returncode == 0


def test_complex_protein_rna():
    complex_protein_rna_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--job_name",
        "pro_RNA",
        "--protein",
        "./test_files/protein-RNA/pro_RNA.fasta",
        "--rna",
        "./test_files/protein-RNA/pro_RNA.fasta",
    ]
    assert run(complex_protein_rna_run_cmd).returncode == 0


def test_complex_protein_ligand_ion():
    complex_protein_ligand_ion_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--job_name",
        "pro_lig_ion",
        "--protein",
        "./test_files/potein-ligand-ion/pro_lig_ion.fasta",
        "--ligand",
        "./test_files/potein-ligand-ion/ligand.txt",  # 使用固定文件
        "--ion",
        "ZN",
    ]
    assert run(complex_protein_ligand_ion_run_cmd).returncode == 0


def test_complex_protein_rna_ligand():
    complex_protein_rna_ligand_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--job_name",
        "pro_RNA_lig",
        "--protein",
        "./test_files/protein-RNA-ligand/pro_RNA_lig.fasta",
        "--rna",
        "./test_files/protein-RNA-ligand/pro_RNA_lig.fasta",
        "--ligand",
        "./test_files/protein-RNA-ligand/ligand.txt",  # 使用固定文件
    ]
    assert run(complex_protein_rna_ligand_run_cmd).returncode == 0
