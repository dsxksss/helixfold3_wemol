#!/alphafold3_venv/bin/python -B
import os
import json
import argparse
from subprocess import Popen, PIPE
import shutil
import sys
from glob import glob
import pandas as pd
from Bio.PDB import MMCIFParser,FastMMCIFParser
from Bio.PDB.PDBIO import PDBIO

def run_ext_cmder(cmds: list, query:str=None) -> str:
    ''' 功能: 执行外部命令行工具
    - cmds: 待执行命令, List格式[cmder, option, value...]
    - query: 可交互执行程序命令组合, 使用\\n标记回车
    '''
    fout = open("out.log","w")
    ferr = open("err.log","w")
    #process = Popen(cmds, stdout=PIPE, stderr=PIPE)
    #process = Popen(cmds, stdout=sys.stdout, stderr=sys.stderr)
    process = Popen(cmds, stdout=fout, stderr=ferr)
    #stdout, stderr = process.communicate() if query is None else process.communicate(query)
    retcode = process.wait()
    fout.close()
    ferr.close()
    if retcode:
        # raise RuntimeError('Run Ext Cmd Failed: \n{}\nSTDOUT:{}\nSTDERR:{}'.format(
        #     process.args, stdout.decode('utf8'), stderr.decode('utf8')
        # ))
        raise RuntimeError('Run Ext Cmd Faile')
    #return stdout.decode('utf8')
    

def cif_to_pdb(cif, pdb):
    # clean LIG_* in cif file genarated by AF3 on SMILES
    os.system(f"sed -i 's/LIG_[A-Z]/LIG  /g' {cif}")
    
    parser = MMCIFParser()
    structure = parser.get_structure("tmp", cif)

    io=PDBIO()
    io.set_structure(structure)
    io.save(pdb)
# cif_to_pdb("/app/test/test_3/rank_1.cif", "rank_1.pdb")
# exit(0)

def load_fasta(fst: str) -> list:
    fasta, sid, seq = [], None, ''
    for line in open(fst, 'r').readlines():
        line = line.strip()
        if line.startswith('>'):
            if len(seq) > 0:
                fasta.append((sid, seq))
                seq = ''
            sid = line[1:]
        else:
            seq += line.upper()
            continue
    if len(seq) > 0:
        fasta.append((sid, seq))
    return fasta

#
# batch format, the seq id should be unique
#
def load_fasta_batch(fst: str) -> list:
    fasta, sid, seq = [], None, ''
    unique_id = []
    for line in open(fst, 'r').readlines():
        line = line.strip()
        if line.startswith('>'):
            if len(seq) > 0:
                fasta.append((sid, seq))
                seq = ''
            sid = line[1:]
            if sid in unique_id:
                sys.stderr.write(f"The seq ID in '{fst}' should be unique in batch format: f{sid}")
            else:
                unique_id.append(sid)
        else:
            seq += line.upper()
            continue
    if len(seq) > 0:
        fasta.append((sid, seq))
    return fasta

#
# format: SMILES or CCD
# one ligand each line
#
'''
CC(=O)OC1C[NH+]2CCC1CC2
CCD,ATP,HY3,......
CCD,P1L
'''
def load_ligands(lig_file: str):
    ligands = []
    with open(lig_file) as fin:
        for line in fin.readlines():
            tmp = line.strip()
            if tmp.startswith("CCD,") or tmp.startswith("ccd,"):  # CCD codes
                tmp = tmp.split(",")
                ligands.append([t.strip() for t in tmp])
            else:
                ligands.append(tmp)
    return ligands
#
# format: UID:SMILES:CCD
# multiple ligands in one complex per line
#
'''
UID:CC(=O)OC1C[NH+]2CCC1CC2:CCD,ATP,HY3:CC1CC(CN2NC[N+](=C3CCC(OC(F)(F)F)CC3)C2=O)CC(C)C1OC(C)(C)C(=O)[O-]
'''
def load_ligands_batch(lig_file: str):
    ligands = {}
    with open(lig_file) as fin:
        for line in fin.readlines():
            ligs = []
            tmp = line.strip().split(":")
            uid = tmp[0].strip()
            for lig in tmp[1:]:
                lig = lig.strip()
                if lig.startswith("CCD,") or lig.startswith("ccd,"): # CCD codes
                    lig = lig.split(",")
                    ligs.append([t.strip() for t in lig])
                else:
                    ligs.append(lig)
            ligands[uid] = ligs
    return ligands

#
# format: Indexofseq,Type,Position
# one modification each line
'''
1,HY3,1
1,P1L,5
2,HY3,3
'''
def load_modification(mdf_file: str):
    mdfs = []
    with open(mdf_file) as fin:
        for line in fin.readlines():
            tmp = line.strip().split(",")
            # check format            
            assert len(tmp) == 3, f"Error format in modifications: {line}"
            seq_index = int(tmp[0])
            type_mdf  = tmp[1].strip()
            pos_aa    = int(tmp[2])
            assert seq_index >=1 and pos_aa >=1 and len(type_mdf)==3, f"Error format in modifications: {line}"
            mdfs.append((seq_index,type_mdf,pos_aa))
    return mdfs

#
# format: UID:Indexofseq,Type,Position:......
# multiple modification in one complex per line
'''
UID:1,HY3,1:1,P1L,5:2,HY3,3
'''
def load_modification_batch(mdf_file: str):
    mdfs_batch = {}
    with open(mdf_file) as fin:
        for line in fin.readlines():
            tmp = line.strip().split(":")
            uid = tmp[0].strip()
            mdfs = []
            for mdf in tmp[1:]:
                mdf = mdf.strip().split(",")
                # check format            
                assert len(mdf) == 3, f"Error format in modifications: {line}"
                seq_index = int(mdf[0])
                type_mdf  = mdf[1].strip()
                pos_aa    = int(mdf[2])
                assert seq_index >=1 and pos_aa >=1 and len(type_mdf)==3, f"Error format in modifications: {line}"
                mdfs.append((seq_index,type_mdf,pos_aa))
            mdfs_batch[uid] = mdfs
    return mdfs_batch

#
# format: Indexofseq,pos,atom;Indexofseq,pos,atom
# one covalent bond per line
#
'''
1,1,CA;2,1,CA
1,1,CA;3,1,CHA
'''
def load_bonds(bond_file: str):
    bonds = []
    with open(bond_file) as fin:
        for line in fin.readlines():
            tmp = line.strip().split(";")
            # check format            
            assert len(tmp) == 2, f"Error format in covalent bonds: {line}"
            left  = tmp[0].strip().split(",")
            right = tmp[1].strip().split(",")
            assert len(left) ==3 and len(right) ==3, f"Error format in covalent bonds: {line}"
            seq_index_l = int(left[0])
            pos_l       = int(left[1])
            atom_name_l = left[2].strip()
            seq_index_r = int(right[0])
            pos_r       = int(right[1])
            atom_name_r = right[2].strip() 
            assert seq_index_l >=1 and pos_l >=1 and seq_index_r>=1 and pos_r>=1, f"Error format in covalent bonds: {line}"
           
            atom_l = (seq_index_l,pos_l,atom_name_l)
            atom_r = (seq_index_r,pos_r,atom_name_r)
            bonds.append((atom_l,atom_r))
    return bonds
# print(load_bonds("bond"))
# exit(0)

#
# format: UID:Indexofseq,pos,atom;Indexofseq,pos,atom:......
# 
'''
UID:1,1,CA;2,1,CA:1,1,CA;3,1,CHA
'''
def load_bonds_batch(bond_file: str):
    bonds_batch = {}
    with open(bond_file) as fin:
        for line in fin.readlines():
            tmp = line.strip().split(":")
            uid = tmp[0].strip()
            bonds = []
            for bond in tmp[1:]:
                atoms = bond.strip().split(";")
                # check format            
                assert len(atoms) == 2, f"Error format in covalent bonds: {bond} with unique id '{uid}'"
                left  = atoms[0].strip().split(",")
                right = atoms[1].strip().split(",")
                assert len(left) ==3 and len(right) ==3, f"Error format in covalent bonds: {bond} with unique id '{uid}'"
                seq_index_l = int(left[0])
                pos_l       = int(left[1])
                atom_name_l = left[2].strip()
                seq_index_r = int(right[0])
                pos_r       = int(right[1])
                atom_name_r = right[2].strip() 
                assert seq_index_l >=1 and pos_l >=1 and seq_index_r>=1 and pos_r>=1, f"Error format in covalent bonds: {bond} with unique id '{uid}'"
            
                atom_l = (seq_index_l,pos_l,atom_name_l)
                atom_r = (seq_index_r,pos_r,atom_name_r)
                bonds.append((atom_l,atom_r))
            bonds_batch[uid] = bonds
    return bonds_batch
# print(load_bonds("bond"))
# exit(0)

chain_MARK = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
def data_to_json(protein, dna, rna, ligand, modification, bond):
    json_data         = {}
    json_data["name"] = "complex"
    json_data["modelSeeds"] = [1,5,39,73,248,333,1970,2024,20967,66666]
    json_data["dialect"]    = "alphafold3"
    json_data["version"]    = 1
    seqs_json = []
    chain_cnt = 0
    if protein:
        pro_seqs_batch = load_fasta(protein)
        for i, pro_seq in enumerate(pro_seqs_batch):
            seq_json = {}
            seq_json["protein"] = {}
            seq_json["protein"]["id"] = chain_MARK[chain_cnt]
            seq_json["protein"]["sequence"] = pro_seq[1]
            chain_cnt += 1
            seqs_json.append(seq_json)
    
    if dna:
        dna_seqs_batch = load_fasta(dna)  
        for i, dna_seq in enumerate(dna_seqs_batch):
            seq_json = {}
            seq_json["dna"] = {}
            seq_json["dna"]["id"] = chain_MARK[chain_cnt]
            seq_json["dna"]["sequence"] = dna_seq[1]
            chain_cnt += 1
            seqs_json.append(seq_json)
            
    if rna:
        rna_seqs_batch = load_fasta(rna)  
        for i, rna_seq in enumerate(rna_seqs_batch):
            seq_json = {}
            seq_json["rna"] = {}
            seq_json["rna"]["id"] = chain_MARK[chain_cnt]
            seq_json["rna"]["sequence"] = rna_seq[1]
            chain_cnt += 1
            seqs_json.append(seq_json)
    
    if ligand:
        ligands = load_ligands(ligand)
        for i, lig in enumerate(ligands):
            if isinstance(lig, list):  # CCD
                for ccd in lig[1:]:
                    seq_json = {}
                    seq_json["ligand"] = {}
                    seq_json["ligand"]["id"] = chain_MARK[chain_cnt]
                    seq_json["ligand"]["ccdCodes"] = [ccd]
                    chain_cnt += 1
                    seqs_json.append(seq_json)
            else:  # smiles
                seq_json = {}
                seq_json["ligand"] = {}
                seq_json["ligand"]["id"] = chain_MARK[chain_cnt]
                seq_json["ligand"]["smiles"] = lig
                chain_cnt += 1
                seqs_json.append(seq_json)
    
    if modification:
        mdfs = load_modification(modification)
        for mdf in mdfs:
            seq_index = mdf[0]
            type_mdf  = mdf[1]
            pos_aa    = mdf[2]
            
            sequence = seqs_json[seq_index-1]
            key = list(sequence.keys())[0]
            if "modifications" not in sequence[key]:
                sequence[key]["modifications"] = []
            
            if key == "protein":
                mdf_json = {}
                mdf_json["ptmType"]     = type_mdf
                mdf_json["ptmPosition"] = pos_aa
                
            elif key in ["dna","rna"]:
                mdf_json = {}
                mdf_json["modificationType"] = type_mdf
                mdf_json["basePosition"]     = pos_aa

            else:
                print(f"Error format in defined modifications :{mdf}")
                exit(-1)
            sequence[key]["modifications"].append(mdf_json)
    json_data["sequences"] = seqs_json    
    # covalent bonds
    if bond:
        cov_bonds = load_bonds(bond)
        bond_json = []
        for covb in cov_bonds:
            atom_l = covb[0]
            atom_r = covb[1]
            seq_index_l = atom_l[0]
            seq_index_r = atom_r[0]
            assert seq_index_l <= chain_cnt and seq_index_r <= chain_cnt, f"Error format in covalent bonds"
            atom_left  = [chain_MARK[seq_index_l-1], atom_l[1], atom_l[2]]
            atom_right = [chain_MARK[seq_index_r-1], atom_r[1], atom_r[2]]
            bond_json.append([atom_left, atom_right])
        json_data["bondedAtomPairs"] = bond_json
      
    return json_data
# json_data = data_to_json("proteins","dna","rna","ligands","modifications","bond")
# with open("input.json","w") as fout:
#     json.dump(json_data,fout,indent=2)
# exit(0)

#
# Batch format:
# 1, fasta files of protein, dna, rna, the sequences with same id will be complex
# 2, txt files of ligands, modifications, bonds, with same uid in the front 
#
# Ligands
'''
UID:CC(=O)OC1C[NH+]2CCC1CC2:CCD,ATP,HY3:CC(=O)OC1C[NH+]2CCC1CC2
'''
# Modifications
'''
UID:1,HY3,1:1,P1L,5:2,HY3,3
'''
# Bonds
'''
UID:1,1,CA;2,1,CA:1,1,CA;3,1,CHA
'''
def batch_fasta_to_json(protein, dna, rna, ligand, modification, bond):
    data_list = {}
    if protein:
        pro_seqs_batch = load_fasta_batch(protein)
        for i, record in enumerate(pro_seqs_batch):
            data_each = {}
            unique_id = record[0].strip()
            data_each["name"] = unique_id
            seqs   = record[1].split(":")
            chain_cnt = 0
            data_each["sequences"] = []
            for seq in seqs:
                seq_data = {}
                seq_data["protein"] = {}
                seq_data["protein"]["id"] = chain_MARK[chain_cnt]
                seq_data["protein"]["sequence"] = seq
                data_each["sequences"].append(seq_data)
                chain_cnt += 1
            
            data_each["modelSeeds"] = [1,5,39,73,248,333,1970,2024,20967,66666]
            data_each["dialect"] = "alphafold3"
            data_each["version"] = 1
            data_list[unique_id] = data_each
    if dna:
        dna_seqs_batch = load_fasta_batch(dna)
        for i, record in enumerate(dna_seqs_batch):
            unique_id = record[0].strip()
            if unique_id in data_list: # already exist
                data_each = data_list[unique_id]
            else:
                data_each = {}
                data_each["name"] = unique_id
                data_each["sequences"] = []
                data_each["modelSeeds"] = [1,5,39,73,248,333,1970,2024,20967,66666]
                data_each["dialect"] = "alphafold3"
                data_each["version"] = 1
                data_list[unique_id] = data_each
            chain_cnt = len(data_each["sequences"])
            seqs   = record[1].split(":")
            for seq in seqs:
                seq_data = {}
                seq_data["dna"] = {}
                seq_data["dna"]["id"] = chain_MARK[chain_cnt]
                seq_data["dna"]["sequence"] = seq
                data_each["sequences"].append(seq_data)
                chain_cnt += 1      
    if rna:
        rna_seqs_batch = load_fasta_batch(rna)
        for i, record in enumerate(rna_seqs_batch):
            unique_id = record[0].strip()
            if unique_id in data_list: # already exist
                data_each = data_list[unique_id]
            else:
                data_each = {}
                data_each["name"] = unique_id
                data_each["sequences"] = []
                data_each["modelSeeds"] = [1,5,39,73,248,333,1970,2024,20967,66666]
                data_each["dialect"] = "alphafold3"
                data_each["version"] = 1
                data_list[unique_id] = data_each
            chain_cnt = len(data_each["sequences"])
            seqs   = record[1].split(":")
            for seq in seqs:
                seq_data = {}
                seq_data["rna"] = {}
                seq_data["rna"]["id"] = chain_MARK[chain_cnt]
                seq_data["rna"]["sequence"] = seq
                data_each["sequences"].append(seq_data)
                chain_cnt += 1     
    if ligand:
        ligands_batch = load_ligands_batch(ligand)
        #print(ligands_batch)    
        for unique_id, ligands in ligands_batch.items():
            # print(unique_id)
            # print(ligands)
            if unique_id in data_list: # already exist
                data_each = data_list[unique_id]
            else:
                data_each = {}
                data_each["name"] = unique_id
                data_each["sequences"] = []
                data_each["modelSeeds"] = [1,5,39,73,248,333,1970,2024,20967,66666]
                data_each["dialect"] = "alphafold3"
                data_each["version"] = 1
                data_list[unique_id] = data_each
            chain_cnt = len(data_each["sequences"])
            for lig in ligands:
                if isinstance(lig, list):  # CCD
                    for ccd in lig[1:]:
                        #print(ccd)
                        seq_data = {}
                        seq_data["ligand"] = {}
                        seq_data["ligand"]["id"] = chain_MARK[chain_cnt]
                        seq_data["ligand"]["ccdCodes"] = [ccd]
                        data_each["sequences"].append(seq_data)
                        chain_cnt += 1  
                else: # smiles
                    seq_data = {}
                    seq_data["ligand"] = {}
                    seq_data["ligand"]["id"] = chain_MARK[chain_cnt]
                    seq_data["ligand"]["smiles"] = lig
                    data_each["sequences"].append(seq_data)
                    chain_cnt += 1      
                
    if modification:
        mdfs_batch = load_modification_batch(modification)
        for unique_id, mdfs in mdfs_batch.items():
            # print(unique_id)
            # print(mdfs)
            assert unique_id in data_list, f"Error on modifications, no unique id: '{unique_id}' found"
            data_each = data_list[unique_id]
            seqs_json = data_each["sequences"]
            for mdf in mdfs:
                seq_index = mdf[0]
                type_mdf  = mdf[1]
                pos_aa    = mdf[2]
                assert seq_index <= len(seqs_json), f"No available sequence found in defined modifications: {mdf} with unique id '{unique_id}'"
                sequence = seqs_json[seq_index-1]
                key = list(sequence.keys())[0]
                if "modifications" not in sequence[key]:
                    sequence[key]["modifications"] = []
                
                if key == "protein":
                    mdf_json = {}
                    mdf_json["ptmType"]     = type_mdf
                    mdf_json["ptmPosition"] = pos_aa
                elif key in ["dna","rna"]:
                    mdf_json = {}
                    mdf_json["modificationType"] = type_mdf
                    mdf_json["basePosition"]     = pos_aa
                else:
                    sys.stderr.write(f"No available sequence found in defined modifications: {mdf} with unique id '{unique_id}'")
                    exit(-1)
                sequence[key]["modifications"].append(mdf_json)          
    if bond:
        bonds_batch = load_bonds_batch(bond)
        for unique_id, bonds in bonds_batch.items():
            # print(unique_id)
            # print(bonds)
            assert unique_id in data_list, f"Error on bonds, no unique id: '{unique_id}' found"
            data_each = data_list[unique_id]
            seqs_json = data_each["sequences"]
            chain_cnt = len(seqs_json)
            bond_json = []
            for bond in bonds:
                atom_l = bond[0]
                atom_r = bond[1]
                seq_index_l = atom_l[0]
                seq_index_r = atom_r[0]
                assert seq_index_l <= chain_cnt and seq_index_r <= chain_cnt, f"Error format in covalent bonds: {bond} with unique id '{unique_id}'"
                atom_left  = [chain_MARK[seq_index_l-1], atom_l[1], atom_l[2]]
                atom_right = [chain_MARK[seq_index_r-1], atom_r[1], atom_r[2]]
                bond_json.append([atom_left, atom_right])
            data_each["bondedAtomPairs"] = bond_json
    return data_list
json_data = batch_fasta_to_json("batch_protein","batch_dna","batch_rna","batch_ligand","batch_modification","batch_bond")
with open("input.json","w") as fout:
    json.dump(json_data,fout,indent=2)
exit(0)

def run_af3(json_file, input_dir=None):
    cmd = []
    cmd.append("python")
    cmd.append("/app/alphafold/run_alphafold.py")
    if json_file:
        cmd.append(f"--json_path={json_file}")
    elif input_dir:
        cmd.append(f"--input_dir={input_dir}")
    else:
        print("Neither json file nor input directory was defined.")
        exit(-1)
    cmd.append("--model_dir=/data/PRG/dbs/af3.0.0/models")
    cmd.append("--db_dir=/data/PRG/dbs/af3.0.0/public_databases")
    cmd.append("--output_dir=./results")
    run_ext_cmder(cmd)
    #print(stdout)

def get_results_single_run(topN: int):
    
    fout = open("scores.csv","w")
    fout.write("Name,Ranking_Score\n")
    
    df = pd.read_csv("results/complex/ranking_scores.csv")
    df = df.sort_values(by='ranking_score', ascending=False, ignore_index=True) # 排序后重新index
    #print(df)
    df = df.iloc[:topN,:]
    df["seed"]   = df["seed"].astype("Int32")
    df["sample"] = df["sample"].astype("Int32")
    for index, row in df.iterrows():
        #print(row["seed"], row["sample"])
        cif_path = f"results/complex/seed-{row['seed']}_sample-{row['sample']}/model.cif"
        #shutil.copy(cif_path,f"rank_{index+1}.cif")
        cif_to_pdb(cif_path,f"rank_{index+1}.pdb")
        fout.write(f"rank_{index+1},{round(row['ranking_score'],4)}\n")
    fout.close()
# get_results_single_run(5)
# exit(0)

def get_results_batch_run(topN: int):
    for root, dirs, files in os.walk("results"):
        for dir in dirs:
            #print(dir)
            df = pd.read_csv(os.path.join("results",dir,"ranking_scores.csv"))
            df = df.sort_values(by='ranking_score', ascending=False, ignore_index=True) # 排序后重新index
            df = df.iloc[:topN,:]
            df["seed"]   = df["seed"].astype("Int32")
            df["sample"] = df["sample"].astype("Int32")

            # output to new path 
            new_path = os.path.join("final_results",dir)
            os.makedirs(new_path,exist_ok=True)
            fout = open(os.path.join(new_path,"scores.csv"),"w")
            fout.write("Name,Ranking_Score\n")
            for index, row in df.iterrows():
                #print(row["seed"], row["sample"])
                cif_path = os.path.join("results",dir,f"seed-{row['seed']}_sample-{row['sample']}","model.cif")
                #shutil.copy(cif_path,os.path.join(new_path,f"rank_{index+1}.cif"))
                cif_to_pdb(cif_path,os.path.join(new_path,f"rank_{index+1}.pdb"))
                fout.write(f"rank_{index+1},{round(row['ranking_score'],4)}\n")
            fout.close()
        break
    # get final_results.tar.gz file 
    cmd = []
    cmd.append("tar")
    cmd.append("-czvf")
    cmd.append("final_results.tar.gz")
    cmd.append("final_results")
    run_ext_cmder(cmd)

# get_results_batch_run(5)
# exit(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=''
    )
    # by json file
    parser.add_argument(
        '--json_file', type=str, default=None, help='input json file for prediction'
    )
    
    # single prediction    
    parser.add_argument(
        '--protein', type=str, default=None, help='fasta file of protein sequences'
    )
    parser.add_argument(
        '--dna', type=str, default=None, help='fasta file of dna sequences'
    )
    parser.add_argument(
        '--rna', type=str, default=None, help='fasta file of rna sequences'
    )
    parser.add_argument(
        '--ligand', type=str, default=None, help='txt file of ligands, one ligands each line'
    )
    parser.add_argument(
        '--modification', type=str, default=None, help='txt file of ligands, one ligands each line'
    )
    parser.add_argument(
        '--cov_bond', type=str, default=None, help='txt file of covalent bonds, one bond each line'
    )
    parser.add_argument(
        '--topN', type=int, default=5, help='topN structure for ranked output'
    )

    
    # batch prediction
    parser.add_argument(
        '--batch_pro', type=str, default=None, help='fasta file for batch  prediction'
    )   
    parser.add_argument(
        '--batch_dna', type=str, default=None, help='fasta file for batch  prediction'
    )   
    parser.add_argument(
        '--batch_rna', type=str, default=None, help='fasta file for batch  prediction'
    )   
    parser.add_argument(
        '--batch_ligand', type=str, default=None, help='file for batch prediction'
    )   
    parser.add_argument(
        '--batch_modification', type=str, default=None, help='file for batch prediction'
    ) 
    parser.add_argument(
        '--batch_cov_bond', type=str, default=None, help='file for batch prediction'
    )  
    parser.add_argument(
        '--max_num', type=int, default=100, help='max number of predictions per batch run'
    )
    
    args = parser.parse_args()
    

    # batch run
    #print("Prepare input data for prediction.")
    if args.batch_pro or args.batch_dna or args.batch_rna or args.batch_ligand: 
        input_data_dict = batch_fasta_to_json(args.batch_pro, args.batch_dna, args.batch_rna, args.batch_ligand, args.batch_modification, args.batch_cov_bond)
        cnt = 0
        os.mkdir("input")
        for unique_id in input_data_dict:
            with open(f"input/{unique_id}.json","w") as fin:
                json.dump(input_data_dict[unique_id], fin, indent=2)
            cnt += 1
            if cnt >= args.max_num: break
        run_af3(None, "input")
        get_results_batch_run(5)
    # single run
    else: 
        json_data = data_to_json(args.protein, args.dna, args.rna, args.ligand, args.modification, args.cov_bond)
        with open("input.json","w") as fout:
            json.dump(json_data,fout,indent=2)
        run_af3("input.json")
        
        # get results
        get_results_single_run(args.topN)
    
