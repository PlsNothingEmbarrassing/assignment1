#!/bin/bash
# Set pipefail parameters
set -eou pipefail

help() {
	echo "Input structure: $0 
	-r1 R1.fastq.gz,
	-r2 R2.fastq.gz,
	-s sample,
	-r reference path,
	-o output directory,
	-a annovar db path"
}

# Set default vals for args (for debugging so I only have to define output dir)
ANNOVAR_DB=/home/scw1557/UNIX_5/annovar/humandb
SAMPLE=my_sample
REF_PATH=/home/scw1557/UNIX_5/BWA
R1=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR525/007/SRR5252327/SRR5252327_1.fastq.gz
R2=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR525/007/SRR5252327/SRR5252327_2.fastq.gz

declare -a READS # Array to hold R1, R2 ... file names
while getopts "r:s:o:a:" opt; do
	case "$opt" in
	r) READS+=$OPTARG ;; # R file paths	
	s) SAMPLE=$OPTARG ;; # sample path
	o) OUTDIR=$OPTARG ;; # Output directory
	a) ANNOVAR_DB=$OPTARG ;; # capture argument for different annovar path if wanted.	
	*) help ;;
	esac
done

# check inputs
if [ "${#READS[@]:-}" -eq 1 ]; # make sure at least 1 read file is provided
then
	echo "Use -r to specify at least 1 read file"
	exit 1
fi
# Create file structure
mkdir -p $OUTDIR/{data,mapping,variants,log} # pass directories as a list is more tidy
LOG_FILE="$OUTDIR/log"
touch $LOG_FILE

# Load ARCCA modules.
module load bwa
module load samtools
module load java

load_data() {
	i = 1 # iterator
	cd "$OUTDIR/data"
	# loop over data array and dowmload each
	for data in "${READS[@]}"; do
		echo -e "Downloading R$i...\n" | tee -a $LOG_FILE # Using -e for newline char (neater than using random echo calls)
		wget $data
		((i++))
	done
	# Set file permissions to read only	
	echo -e "Setting file permissions to read only...\n" | tee -a $LOG_FILE
	chmod 444 
	# Create symbolic links to mapping dir
	echo -e "Creating symbolic links...\n" | tee -a $LOG_FILE
	cd ../mapping
	j=1
	for path in "${READS[@]}"; do
		basename "$path"
		ln -s ../data/$path R$j.fastq.gz
		((j++))
	done
	
}

# Set file permissions to read only
# chmod 444 SRR5252327_1.fastq.gz
# chmod 444 SRR5252327_2.fastq.gz
# cd ..
# Create mapping dir and symbolic links to read only files
# mkdir mapping
# cd mapping
# ln -s ../data/SRR5252327_1.fastq.gz R1.fastq.gz
# ln -s ../data/SRR5252327_2.fastq.gz R2.fastq.gz

# echo
# gunzip -c R1.fastq.gz | grep -cP "^@\w+\.\w+"| xargs echo "Number of reads in R1:"
# gunzip -c R2.fastq.gz | grep -cP "^@\w+\.\w+"| xargs echo "Number of reads in R2:"
# echo
# echo
# echo "Mapping data..."
map_data() {
	# capture data list passed to function
	i=1
	local data_list=("$@")	
	for item in "${data_list[@]}"; do
		gunzip -c $item | grep -cP "^@\w+\.\w+"| xargs echo "Number of reads in R$i:"
		((i++))
	done
	
}
bwa mem -R "@RG\tID:$SAMPLE\tPL:illumina\tSM:$SAMPLE" \
$REF_PATH/hg38.fasta \
R1.fastq.gz \
R2.fastq.gz \
> $SAMPLE\.sam \
2>> ../log
echo "Done!"
echo
echo "Converting SAM to BAM..."
samtools view -bt $REF_PATH/hg38.fasta $SAMPLE\.sam | samtools sort - | tee $SAMPLE\.bam | samtools index - $SAMPLE\.bam.bai
echo
echo "Summarise mapping..."
samtools view $SAMPLE\.bam | cut -f3 | uniq -c | sort -n
cd ..
mkdir variants
cd variants
echo
echo "Run samtools mpileup..."
samtools mpileup -B -f $REF_PATH/hg38.fasta ../mapping/$SAMPLE\.bam > $SAMPLE\.pileup 2>> ../log
echo "Done!"
echo
echo "Summmarise number of base positions covered in each chromosome..."
cut -f1 $SAMPLE\.pileup | uniq -c
echo
echo "Representation according to base..."
cut -f3 $SAMPLE\.pileup | tr 'acgt' 'ACGT' | sort | uniq -c | sort -nr
echo
echo -e "depth\tlocations" > coverage_depth.txt
cut -f4 $SAMPLE\.pileup | sort | uniq -c | sort -n -k2 | awk ' { t = $1; $1 = $2; $2 = t; print $1 "\t" $2; } ' >> coverage_depth.txt
echo
echo "Call variants with varscan..."
java -jar /software/genomics/varscan/2.4.1/VarScan.jar mpileup2snp $SAMPLE\.pileup --output-vcf 1 > $SAMPLE\.vcf 2>> ../log
echo "Done!"
echo
echo "Total number of SNPs identified:"
grep -vc "^#" $SAMPLE\.vcf
echo "Creating a simple spreadsheet from the VCF file..."
grep -v "^##" $SAMPLE\.vcf | cut -f1,2,4,5 | sed 's/\t/,/g' > snps.csv
echo "Done!"
echo
echo "SNP locations:"
grep -v "^#" $SAMPLE\.vcf | cut -f1 | uniq -c | perl -ne '/\s+(\d+)\s+(\d+)/; print "$1 snps in chromosome $2\n"'
echo
echo "Annotate..."
/home/scw1557/UNIX_5/annovar/convert2annovar.pl $SAMPLE\.vcf -format vcf4old --includeinfo > $SAMPLE\.av 2>> ../log
/home/scw1557/UNIX_5/annovar/annotate_variation.pl -buildver hg38 -geneanno -dbtype knownGene $SAMPLE\.av $ANNOVAR_DB --outfile annovar_out
2>> ../log
echo -e "chr\tposition\tref\talt\tgene\texon\tcDNA pos\taa pos\tmutation type" > final.txt
cut -f2 annovar_out.exonic_variant_function | cut -f1 -d' ' > cols2.txt
cut -f3 annovar_out.exonic_variant_function | cut -f1,3,4,5 -d':' --output-delimiter=" " | tr ' ' '\t' > cols3.txt
cut -f4,5,7,8 annovar_out.exonic_variant_function > cols4578.txt
paste cols4578.txt cols3.txt cols2.txt >> final.txt
chmod 444 final.txt
cat final.txt
# Leave directory.
cd ..
# Leave working directory
cd ..
echo "Create tarball from working directory..."
tar czf $OUTDIR\.tar.gz $OUTDIR
echo
echo "Complete!"