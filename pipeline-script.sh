#!/bin/bash
# Set pipefail parameters
set -eou pipefail

help() {
	echo "Input structure: $0 
	-r1 R1.fastq.gz,
	-r2 R2.fastq.gz	
	-s sample,
	-r reference path,
	-o output directory,
	-a annovar db path"
}

# Set default vals for args (for debugging so I only have to define output dir)
ANNOVAR_DB=/home/scw1557/UNIX_5/annovar/humandb
SAMPLE=my_sample
REF_PATH=/home/scw1557/UNIX_5/BWA
OUTDIR="output"
R1=""
R2=""
#R1=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR525/007/SRR5252327/SRR5252327_1.fastq.gz
#R2=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR525/007/SRR5252327/SRR5252327_2.fastq.gz

# Parse arguments
while getopts "1:2:s:o:a:" opt; do
	case "$opt" in
		1) R1=$OPTARG ;;  # R1 filepath
		2) R2=$OPTARG ;;  # R2 filepath	
		s) SAMPLE=$OPTARG ;;  # sample path
		o) OUTDIR=$OPTARG ;;  # output directory
		a) ANNOVAR_DB=$OPTARG ;;  # optional annovar path
		*) echo "Invalid option" && exit 1 ;;
	esac
done

# Check if R1 is provided
if [[ -z "$R1" ]]; then
	echo "Error: Use -1 to specify at least one read file (R1)"
	exit 1
fi

# Load modules (adjust as necessary for your environment)
module load bwa
module load samtools
module load java

# Function to download data and set up symbolic links
load_data() {
	# Create file structure
	mkdir -p "$OUTDIR"/{data,mapping,variants}  # create necessary directories
	LOG_FILE="../log"
	touch "$LOG_FILE"  # Create the log file if it doesn't exist
	
	cd "$OUTDIR/data" || exit 1

	# Download each file
	echo "Downloading R1..." | tee -a "$LOG_FILE"
	wget "$R1" || { echo "Error downloading R1 file" && exit 1; }
	# If R2 is provided, download it
	if [[ -n "$R2" ]]; then
		echo "Downloading R2..." | tee -a "$LOG_FILE"
		wget "$R2" || { echo "Error downloading R2 file" && exit 1; }
	fi

	# Set file permissions to read-only
	echo -e "Setting file permissions to read-only...\n" | tee -a "$LOG_FILE"
	chmod 444 *.fastq.gz

	# Create symbolic links in the mapping directory
	echo -e "Creating symbolic links...\n" | tee -a "$LOG_FILE"
	cd ..
	cd "mapping" || exit 1
	R1_BASE=$(basename "$R1")
	ln -s ../data/$R1_BASE R1.fastq.gz

	if [[ -n "$R2" ]]; then
		R2_BASE=$(basename "$R2")	
		ln -s ../data/$R2_BASE R2.fastq.gz
	fi
	
}

map_data() {
	# Unzip data

	gunzip -c R1.fastq.gz | grep -cP "^@\w+\.\w+"| xargs echo "Number of reads in R1:"
	# If R2 is provided, unzip and count reads
	if [[ -n "$R2" ]]; then
		gunzip -c R2.fastq.gz | grep -cP "^@\w+\.\w+"| xargs echo "Number of reads in R2:"
	fi
	
	# Map data
	echo -e "\nMapping data..." | tee -a "$LOG_FILE"
	# If R2 is provided, map with paired-end reads
	if [[ -n "$R2" ]]; then
		bwa mem -R "@RG\tID:$SAMPLE\tPL:illumina\tSM:$SAMPLE" \
		$REF_PATH/hg38.fasta \
		R1.fastq.gz \
		R2.fastq.gz \
		> $SAMPLE\.sam \
		2>> "$LOG_FILE"
	else
		bwa mem -R "@RG\tID:$SAMPLE\tPL:illumina\tSM:$SAMPLE" \
		$REF_PATH/hg38.fasta \
		R1.fastq.gz \
		> $SAMPLE\.sam \
		2>> "$LOG_FILE"

	fi

	echo -e "Done!\n"
	# Convert SAM to BAM
	echo -e "Converting SAM to BAM...\n" | tee -a "$LOG_FILE"
	samtools view -bt $REF_PATH/hg38.fasta $SAMPLE\.sam | samtools sort - | tee $SAMPLE\.bam | samtools index - $SAMPLE\.bam.bai
	# Summarise mapping
	echo "Summarise mapping..."
	samtools view $SAMPLE\.bam | cut -f3 | uniq -c | sort -n

	# Check if the BAM index file is correct
	if samtools idxstats $SAMPLE\.bam > /dev/null 2>&1; then
		echo "BAM index file is correct." | tee -a "$LOG_FILE"
	else
		echo "Error: BAM index file is not correct." | tee -a "$LOG_FILE"
		exit 1
	fi


}

summarise_coverage(){
	cd ..
	cd variants
	echo -e "Running samtools mpileup..." | tee -a $LOG_FILE
	# Run samtools mpileup to generate a pileup file
	samtools mpileup -B -f $REF_PATH/hg38.fasta ../mapping/$SAMPLE\.bam > $SAMPLE\.pileup 2>> ../log
	echo "Done!" | tee -a $LOG_FILE
	
	echo -e "\nSummarise number of base positions covered in each chromosome..." | tee -a $LOG_FILE
	# Summarise the number of base positions covered in each chromosome
	cut -f1 $SAMPLE\.pileup | uniq -c
	
	echo -e "\nRepresentation according to base..." | tee -a $LOG_FILE
	# Summarise the representation of each base (A, C, G, T) in the pileup file
	cut -f3 $SAMPLE\.pileup | tr 'acgt' 'ACGT' | sort | uniq -c | sort -nr
	
	# Generate a coverage depth summary
	echo -e "depth\tlocations" > coverage_depth.txt
	cut -f4 $SAMPLE\.pileup | sort | uniq -c | sort -n -k2 | awk ' { t = $1; $1 = $2; $2 = t; print $1 "\t" $2; } ' >> coverage_depth.txt
}

annotate_variation(){
	cd ..
	cd variants
	echo "Call variants with varscan..."
	# Call variants using VarScan
	java -jar /software/genomics/varscan/2.4.1/VarScan.jar mpileup2snp $SAMPLE\.pileup --output-vcf 1 > $SAMPLE\.vcf 2>> ../log
	echo "Done!"
	
	# Count the total number of SNPs identified
	echo "Total number of SNPs identified:"
	grep -vc "^#" $SAMPLE\.vcf
	
	# Create a simple spreadsheet from the VCF file
	echo "Creating a simple spreadsheet from the VCF file..."
	grep -v "^##" $SAMPLE\.vcf | cut -f1,2,4,5 | sed 's/\t/,/g' > snps.csv
	echo "Done!"
	
	# Display SNP locations
	echo
	echo "SNP locations:"
	grep -v "^#" $SAMPLE\.vcf | cut -f1 | uniq -c | perl -ne '/\s+(\d+)\s+(\d+)/; print "$1 snps in chromosome $2\n"'
	echo
	
	# Annotate variants using ANNOVAR
	echo "Annotate..."
	/home/scw1557/UNIX_5/annovar/convert2annovar.pl $SAMPLE\.vcf -format vcf4old --includeinfo > $SAMPLE\.av 2>> ../log
	/home/scw1557/UNIX_5/annovar/annotate_variation.pl -buildver hg38 -geneanno -dbtype knownGene $SAMPLE\.av $ANNOVAR_DB --outfile annovar_out 2>> ../log
	
	# Create a final annotated file
	echo -e "chr\tposition\tref\talt\tgene\texon\tcDNA pos\taa pos\tmutation type" > final.txt
	cut -f2 annovar_out.exonic_variant_function | cut -f1 -d' ' > cols2.txt
	cut -f3 annovar_out.exonic_variant_function | cut -f1,3,4,5 -d':' --output-delimiter=" " | tr ' ' '\t' > cols3.txt
	cut -f4,5,7,8 annovar_out.exonic_variant_function > cols4578.txt
	paste cols4578.txt cols3.txt cols2.txt >> final.txt
	
	# Set final file permissions to read-only and display the content
	chmod 444 final.txt
	cat final.txt
	chmod 444 final.txt

	
}
compare_varscan(){
	# Run varscan with different min coverage values and compare outputs
	for min_coverage in 10 20 30; do
		echo "Running varscan with min coverage $min_coverage..."
		java -jar /software/genomics/varscan/2.4.1/VarScan.jar mpileup2snp $SAMPLE\.pileup --min-coverage $min_coverage --output-vcf 1 > $SAMPLE\_cov$min_coverage\.vcf 2>> ../log
		echo "Done!"
	done

	echo "Comparing outputs with diff..."
	for min_coverage in 10 20 30; do
		if [[ $min_coverage -ne 10 ]]; then
			prev_coverage=$((min_coverage - 10))
			echo "Comparing coverage $prev_coverage and $min_coverage..."
			diff $SAMPLE\_cov$prev_coverage\.vcf $SAMPLE\_cov$min_coverage\.vcf > diff_cov$prev_coverage\_cov$min_coverage\.txt
			echo "Done!"
		fi
	done
}
rm -r output
load_data
map_data 
summarise_coverage
annotate_variation




# # Leave directory.
# cd ..
# # Leave working directory
# cd ..
# echo "Create tarball from working directory..."
# tar czf $OUTDIR\.tar.gz $OUTDIR
# echo
# echo "Complete!"