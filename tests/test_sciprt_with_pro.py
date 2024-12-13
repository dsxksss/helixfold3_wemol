from subprocess import run


def test_protein_small():
    protein_small_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--job_name",
        "pro_300AA",
        "--protein",
        """>pro_300AA.A
GPVIRQGPVNQTVAVDGTFVLSCVATGSPVPTILWRKDGVLVSTQDSRIKQLENGVLQIRYAKLGDTGRYTCIASTPSGEATWSAYIEVQ
>pro_300AA.B
EVQLVESGGGVVQPGGSLKLSCAASGFTFSTYDMSWVRQTPDKRLELVATINSNGGSTYYPDSV
KGRFTSSRDNAKNILYLQMSSLKSEDTAMYYCAREALLRPPYYALDYWGQGTSVTVS
>pro_300AA.C
LDIQMTQSPASLSASVGETVTITCGASENIYGALTWYQRKQGKSPQLLIYGAINLADDKSSRFSGSGSGRQYSLKISSLHPDDVATYYCQNVLSTPFTFGSGTKLEIK""",
    ]
    assert run(protein_small_run_cmd, shell=True).returncode == 0


def test_protein_medium():
    protein_medium_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--job_name",
        "pro_600AA",
        "--protein",
        """>pro_600AA.A
VSGITALTVVVGTVIGAGIFFKPTAVYGAAGAPGLGLLAWFVAGIITIAGGLTVAEIGTIYPQTGGMMIYLEKVYGRWLGFLVGWAQMVIYYPANIAALAIIFATQFVNLFALSDSTIVPTAILTSIFLMGVNFLGTKYSGWIQTLATILKLIPLVVIIVAGLLYPGGGVIRLVPFSVETHPVLTSFGSALIATLFAYDGWINVGTLAGEMKNPGKMLPKVIIGGLSIVMAVYLLTNIAYLFVLDSSQLAGTDTPAALVASHLFEGIGSKLVTIGILISVFGGINGYIISGLRVPYALATQKMLPFSDWFARINPKTNLPINGGLVMLGIAIVMILTGQFNQLTDLIVFVIWFFITLTFIAVIILRKTQPDIERPYRVPFYPVIPLIAIIGGLYIIFNTLIVQPKNAFIGILLTLIGIPIYFY
CKKKYGS
>pro_600AA.B
QVQLVESGGGVVQAGGSLRLSCAASGRTFSSRAMGWFRQAPGEGREFVATISWSGSYTEYADSVKGRVTISRDNAKNTVYLQMNSLKPGDTAVYHCAAKNGGAASNYPNDYVYWGQGTQVTVSSHHHHHHE""",
    ]
    assert run(protein_medium_run_cmd, shell=True).returncode == 0


def test_complex_protein_dna():
    complex_protein_dna_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--protein",
        """>pro_DNA.A
ASSINPWILTGFADAEGSFGLYIINRNRGRIRYTTRLKFTITLHNKDKSILENIQSTWKVGI""",
        "--dna",
        """>pro_DNA.B
gggaatggcagtattcatccacaatg
>pro_DNA.C
ccattgtggatgaatactgccattcc""",
    ]
    assert run(complex_protein_dna_run_cmd, shell=True).returncode == 0


def test_complex_protein_ligand():
    complex_protein_ligand_run_cmd = [
        "python",
        "/data/PRG/tools/helixfold3/apps/helixfold3_wemol/src/py_script_template/__main__.py",
        "--protein",
        """>pro_lig.A
ADLKAFSKHIYNAYLKNFNMTKKKARSILTGKASHTAPFVIHDIETLWQAEKGLVWKQLVNGLPPYKEISVHVFYRCQCTTVETVRELTEFAKSIPSFSSLFLNDQVTLLKYGVHEAIFAMLASIVNKDGLLVANGSGFVTREFLRSLRKPFSDIIEPKFEFAVKFNALELDDSDLALFIAAIIILCGDRPGLMNVPRVEAIQDTILRALEFHLQANHPDAQYLFPKLLQKMADLRQLVTEHAQMMQRIKKTETETSLHPLLQEIYKDMY""",
        "--ligand",
        """CC(=O)OC1C[NH+]2CCC1CC2
CCD,ATP,HY3""",
    ]
    assert run(complex_protein_ligand_run_cmd, shell=True).returncode == 0
