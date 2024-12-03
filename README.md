# HelixFold3-WeMol

这是一个用于HelixFold3预测的WeMol端模块。用于处理蛋白质、DNA、RNA和配体数据，并通过HelixFold3进行预测。

## 使用

### 命令行参数

您可以通过命令行运行此脚本，并传递以下参数：

- `--protein`: 蛋白质序列，FASTA格式
- `--dna`: DNA序列，FASTA格式
- `--rna`: RNA序列，FASTA格式
- `--ligand`: 配体和离子信息，每行一个
- `--topN`: 输出排名前N的结构，默认值为5
- `--recycle`: 模型推理的循环参数，默认值为10
- `--ensemble`: 模型推理的集成参数，默认值为1

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
