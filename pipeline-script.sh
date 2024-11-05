#!/bin/bash
# Set pipefail parameters
set -eou pipefail

help() {
	echo "Input structure: $0 
	-1 R1.fastq.gz,
	-2 R2.fastq.gz (Optional),	
	-s sample (Optional, will use default if not provided),
	-r reference path (Optional, will use default if not provided),
	-o output directory (Optional, will use default if not provided),
	-a annovar db path (Optional, will use default if not provided)"
}

# Set default vals for args (for debugging so I only have to define input files)
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
		*) help ;;
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

annotate_variants(){
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
	# Get column containing info, discard irrelevant info. 
	# Split columns by ; delimiter, assign 1/1, 0/1 based on the values of the columns
	# Paste the final file with the temp file to make sure move is successful before overwriting
	# Add column name "Genotype" to snps.csv	
	# Generate the genotype data column (e.g., "1/1" or "0/1") and save it as a new column file
	grep -v "^#" my_sample.vcf | cut -f8 | \
	cut -f3,4 -d';' | \
	awk -F';' '{
		split($1, a, "="); 
		split($2, b, "="); 
		if (b[2] == 1) 
			print "1/1"; 
		else if (a[2] == 1) 
			print "0/1"; 
		else 
			print "0/1";
	}' > genotype_data.txt

	# Add the "Genotype" header to the genotype_data.txt file
	echo "Genotype" | cat - genotype_data.txt > temp_genotype_data.txt && mv temp_genotype_data.txt genotype_data.txt

	# Combine the header with the existing snps.csv using paste
	paste -d',' snps.csv genotype_data.txt > temp_snps.csv && mv temp_snps.csv snps.csv	
	
	# Set final file permissions to read-only and display the content
	chmod 444 final.txt
	cat final.txt
	chmod 444 final.txt

	
}
compare_varscan(){
	# Run varscan with different min coverage values and compare outputs
	for min_coverage in 5 10 15 20; do
		java -jar /software/genomics/varscan/2.4.1/VarScan.jar mpileup2snp $SAMPLE\.pileup --min-coverage $min_coverage --output-vcf 1 > $SAMPLE\_mincov$min_coverage\.vcf 2>> ../log
		echo "Done with min coverage $min_coverage"
	done
	
	# Compare the number of SNPs identified for each min coverage value
	for min_coverage in 5 10 15 20; do
		echo "Total number of SNPs identified with min coverage $min_coverage:"
		grep -vc "^#" $SAMPLE\_mincov$min_coverage\.vcf
	done

	reference_file=$SAMPLE.vcf

	# Loop through all VCF files in the current directory
	for vcf_file in *.vcf; do
		# Skip the reference file
		if [[ "$vcf_file" == "$reference_file" ]]; then
			continue
		fi

		echo "Comparing $vcf_file to $reference_file"
		vcf_content=$(grep -v "^#" $vcf_file)
		reference_content=$(grep -v "^#" $reference_file)
		diff_ouput=""
		# Compare the current VCF file to the reference and save differences to a temp file
		diff_output=$(diff <(echo "$reference_content") <(echo "$vcf_content") || true)
		
		# Check if there are differences
		if [ -z "$diff_output" ]; then
		 	echo "No differences found between $vcf_file and $reference_file."
		else
			echo "Differences found between $vcf_file and $reference_file. Saving to ${vcf_file%.vcf}_diff.txt"
			echo -e "$diff_output\n" > "${vcf_file%.vcf}_diff.txt"
		fi
	done
}

clean_up(){
	# remove unnecessary files, keep only the final output files
	rm cols2.txt cols3.txt cols4578.txt ../mapping/my_sample.sam my_sample.av my_sample.pileup annovar_out.exonic_variant_function annovar_out.variant_function genotype_data.txt
	#rm my_sample_mincov5.vcf my_sample_mincov10.vcf my_sample_mincov15.vcf my_sample_mincov20.vcf
	# Coverage depth file should be kept for quality control purposes
	# vcf file should be kept for comparison purposes
	# # Leave directory.
	cd ..
	# Leave working directory
	cd ..
	echo "Create tarball from working directory..."
	tar czf $OUTDIR\.tar.gz $OUTDIR
	echo
	echo "Complete!"
}
rm -r output
load_data
map_data 
summarise_coverage
annotate_variants

compare_varscan 
clean_up



