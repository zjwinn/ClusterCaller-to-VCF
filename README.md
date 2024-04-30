# Introduction
ClusterCaller-to-VCF is a shell script written to take the output from ClusterCaller software and convert it into a compressed variant calling format file (vcf.gz). This VCF is compiant with formating requirments and is readable in other programs like [TASSEL](https://www.maizegenetics.net/tassel).

# Required inputs
ClusterCaller-to-VCF requires three inputs:
1. A ClusterCaller formatted output - [example](https://github.com/zjwinn/ClusterCaller-to-VCF/blob/main/clustercaller_example.txt)
2. A key file which relates ClusterCaller output to allelic states - [example](https://github.com/zjwinn/ClusterCaller-to-VCF/blob/main/keyfile_example.txt)
3. A string for output name (i.e., "output_example")

Both files provided to the ClusterCaller-to-VCF function must be tab delimited and marker names are case-sensitive and must match exactly. 

# Usage
To call on the ClusterCaller-to-VCF function, the user must call on the [clustercaller_to_vcf.sh](https://github.com/zjwinn/ClusterCaller-to-VCF/blob/main/clustercaller_to_vcf.sh) file directly using the following argument:
```bash
bash clustercaller_to_vcf.sh
```
This can either be in the directy where you have pulled the repository or a direct path to the location of the installation. To view the usage file for futher assistance use the following command
```bash
bash clustercaller_to_vcf.sh --help
```
