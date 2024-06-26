# Introduction
KlusterCaller-to-VCF is a shell script written to take the output from KlusterCaller software and convert it into a compressed variant calling format file (vcf.gz) using base R data manipulation functions. This VCF is compliant with formatting requirements and is readable in other programs like [TASSEL](https://www.maizegenetics.net/tassel).

# Required Inputs
KlusterCaller-to-VCF requires three inputs:
1. A KlusterCaller formatted output - [example](https://github.com/zjwinn/KlusterCaller-to-VCF/blob/main/klustercaller_file_example.txt)
2. A keyfile which relates KlusterCaller output to allelic states - [example](https://github.com/zjwinn/KlusterCaller-to-VCF/blob/main/klustercaller_keyfile_example.txt)
3. A string for output name (i.e., "output_example")

Both files provided to the KlusterCaller-to-VCF function must be tab delimited and marker names are case-sensitive and must match exactly. 

# Usage
To call on the KlusterCaller-to-VCF function, the user must call on the [klustercaller_to_vcf.sh](https://github.com/zjwinn/ClusterCaller-to-VCF/blob/main/klustercaller_to_vcf.sh) file directly using the following argument:
```bash
bash klustercaller_to_vcf.sh
```
This can either be in the directory where you have pulled the repository or a direct path to the location of the installation. To view the usage file for further assistance use the following command:
```bash
bash klustercaller_to_vcf.sh --help
```
To run the provided example and look at textual output use the following command:
```bash
bash klustercaller_to_vcf.sh \
  -k klustercaller_keyfile_example.txt \
  -c klustercaller_file_example.txt \
  -o example_output \
  -v
```
This will result in a compressed VCF file that looks like [this example](https://github.com/zjwinn/KlusterCaller-to-VCF/blob/main/example.vcf.gz).

# Package Requirements
The bash script provided in this GitHub repository requires the following to run properly:
1. R statistical coding language - [link to page](https://www.r-project.org/)
2. samtools/bcftools - [link to page](https://samtools.github.io/bcftools/)

Please make sure you have both programs installed properly prior to running to avoid errors.
