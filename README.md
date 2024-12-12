# HelixFold3-WeMol

百度螺旋桨 PaddleHelix 团队研发的 HelixFold3，在常规的小分子配体、核酸分子（包括 DNA 和 RNA）以及蛋白质的结构预测精度上已与 AlphaFold3 相媲美

## 使用

### 命令行参数

您可以通过命令行运行此脚本，并传递以下参数：

- `--job_name`: 任务名称，默认值为complex
- `--protein`: 蛋白质序列，FASTA格式字符串
- `--dna`: DNA序列，FASTA格式字符串
- `--rna`: RNA序列，FASTA格式字符串
- `--ligand`: 配体和离子信息，每行一个
- `--recycle`: 模型推理的循环参数，默认值为10
- `--ensemble`: 设定结构预测数量，默认值为5

### 示例

以下是一些使用示例：

#### 预测蛋白质

```bash
python -m py_script_template --protein ">pro_300AA.A\nGPVIRQGPVNQTVAVDGTFVLSCVATGSPVPTILWRKDGVLVSTQDSRIKQLENGVLQIRYAKLGDTGRYTCIASTPSGEATWSAYIEVQ"
```

#### 预测复合物

```bash
python -m py_script_template --protein ">pro_lig.A\nADLKAFSKHIYNAYLKNFNMTKKKARSILTGKASHTAPFVIHDIETLWQAEKGLVWKQLVNGLPPYKEISVHVFYRCQCTTVETVRELTEFAKSIPSFSSLFLNDQVTLLKYGVHEAIFAMLASIVNKDGLLVANGSGFVTREFLRSLRKPFSDIIEPKFEFAVKFNALELDDSDLALFIAAIIILCGDRPGLMNVPRVEAIQDTILRALEFHLQANHPDAQYLFPKLLQKMADLRQLVTEHAQMMQRIKKTETETSLHPLLQEIYKDMY" --ligand "CC(=O)OC1C[NH+]2CCC1CC2\nCCD,ATP,HY3"
```
