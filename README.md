# Introduction
ClusterCaller-to-VCF is a shell script written to take the output from ClusterCaller software and convert it into a compressed variant calling format file (vcf.gz) using base R data manipulation functions. This VCF is compliant with formatting requirements and is readable in other programs like [TASSEL](https://www.maizegenetics.net/tassel).

# Required Inputs
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
This can either be in the directory where you have pulled the repository or a direct path to the location of the installation. To view the usage file for further assistance use the following command:
```bash
bash clustercaller_to_vcf.sh --help
```
To run the provided example and look at textual output use the following command:
```bash
bash clustercaller_to_vcf.sh \
  -k keyfile_example.txt \
  -c clustercaller_example.txt \
  -o output_example \
  -v
```
This will result in a compressed VCF file that looks like [this example](https://github.com/zjwinn/ClusterCaller-to-VCF/blob/main/output_example.vcf.gz).

# Package Requirements
The bash script provided in this GitHub repository requires the following to run properly:
1. R statistical coding language - [link to page](https://www.r-project.org/)
2. samtools/bcftools - [link to page](https://samtools.github.io/bcftools/)

Please make sure you have both programs installed properly prior to running to avoid errors.
