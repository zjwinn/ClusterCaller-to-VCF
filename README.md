# Introduction
ClusterCaller-to-VCF is a shell script written to take the output from ClusterCaller software and convert it into a compressed variant calling format file (vcf.gz). This VCF is compiant with formating requirments and is readable in other programs like [TASSEL](https://www.maizegenetics.net/tassel).

# Required inputs
ClusterCaller-to-VCF requires three inputs:
1. A ClusterCaller formatted output - [example](https://github.com/zjwinn/ClusterCaller-to-VCF/blob/main/clustercaller_example.txt)
2. A key file which relates ClusterCaller output to allelic states - [example](https://github.com/zjwinn/ClusterCaller-to-VCF/blob/main/keyfile_example.txt)
3. A string for output name (i.e., "output_example")

Both files provided to the ClusterCaller-to-VCF function must be tab delimited and marker names are case-sensitive and must match exactly. 
