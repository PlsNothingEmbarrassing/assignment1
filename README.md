# Script writeup
## Task 1
### Opt Args for CLI inputs

For flexibility and reliability I decided to use "opt args", or "flags" for reducing the hardcoded elements of the script such as the input files or output directory.
```bash
# Parse arguments

while getopts "1:2:s:o:a:" opt; do
    case "$opt" in
        1) R1=$OPTARG ;;  # R1 filepath

        2) R2=$OPTARG ;;  # R2 filepath

        s) SAMPLE=$OPTARG ;;  # sample path

        o) OUTDIR=$OPTARG ;;  # output directory

        a) ANNOVAR_DB=$OPTARG ;;  # optional annovar path

        *) help ;;
    esac
done
```

This allows for arguments to be input in any order in the command line as long as the appropriate flag is set which could minimise errors from misplaced positional arguments in other approaches.

If the user inputs an incorrect flag, it calls the "help" function which outputs the flags, their format, and a short explanation of each to help the user fix their input.
```bash
help() {

    echo "Input structure: $0

    -1 R1.fastq.gz,

    -2 R2.fastq.gz (Optional),  

    -s sample (Optional, will use default if not provided),

    -r reference path (Optional, will use default if not provided),

    -o output directory (Optional, will use default if not provided),

    -a annovar db path (Optional, will use default if not provided)"
}
```
It also meant that I could set default values for the arguments, such as output folder and annovar db path.  This meant not having to specify them in my CLI input unless I wanted to use a different value than the default. This was very useful for debugging so I would not have to write out all the arguments.
This is also useful for anyone in the future wanting to edit use the script, as they could modify the default values to save themselves time.

I also considered having the script able to take an unlimited number of read files as I could pass them as a list on the command line, but realised that the number of input files was limited to 1-2 so decided that it would be unnecessary so limited the input files to a max of 2. Shown below was the original code. 

```bash
declare -a READS # Array to hold R1, R2 ... file names
while getopts "1:2:s:o:a:" opt; do
	`r) READS+=$OPTARG ;; # R file paths	
	`s) SAMPLE=$OPTARG ;; # sample path
	`o) OUTDIR=$OPTARG ;; # Output directory
	`a) ANNOVAR_DB=$OPTARG ;; # capture argument for different annovar path if wanted.
	`*) help ;;
```

---
### Splitting code into functions

I split the code into 6 separate functions. This made the code clearer, as each step of the process was distinct. Additionally it made debugging the code easier as each function could be run separately (provided that the necessary intermediary files were present) which saved time when testing each function as the previous steps as time consuming steps (such as mapping the data) would not need to be repeated each time.

1. **load_data** - This function contained the code for creating the necessary directories (data/mapping/variants), creating a log file and downloading the data. The only real changes I made to the script was adding conditional statements so that the script could run with 1 or 2 read files as shown below.  This function also contains the code for creating the symbolic links to the files in the mapping folder.
```bash
# If R2 is provided, download it

    if [[ -n "$R2" ]]; then

        echo "Downloading R2..." | tee -a "$LOG_FILE"

        wget "$R2" || { echo "Error downloading R2 file" && exit 1; }

    fi
```

2. **map_data** - This function contains the code for decompressing and mapping the data using bwa. Then converting the SAM file to BAM. Similarly to the previous function, most of the changes in this function are using conditional statements so that the script can run with 1 or 2 read files. Other changes are making sure any outputs are also being output to the log file.

3. **summarise_coverage** - This function contains code which is entirely unchanged other than the additional outputting to log file. The code analyses the coverage of the sequencing data.
4. **annotate_variants** - This function contains the majority of the changes, such as the extraction of "Genotype" from the generated vcf file. This will be discussed more in depth in the third task section.
5. **compare_varscan** - This function contains code for generating and comparing vcf files with different min_coverage settings and comparing the consequent outputs. (Code below showing simple code to output the number of snps for each run.)

```bash
# Compare the number of SNPs identified for each min coverage value

    for min_coverage in 5 10 15 20; do

        echo "Total number of SNPs identified with min coverage $min_coverage:"

        grep -vc "^#" $SAMPLE\_mincov$min_coverage\.vcf

    done
```

6. **clean_up** - This function contains the code to remove intermediary files and creating an archive file of the output.

![flowchart](https://github.com/PlsNothingEmbarrassing/assignment1/blob/master/flow-diagram.png?raw=true "flowchart")
---
### Sanity checking inputs

In order to make sure the user has provided at least 1 input file, there is a simple sanity check to make sure one has been supplied:
```bash
# Check if R1 is provided
if [[ -z "$R1" ]]; then
	echo "Error: Use -1 to specify at least one read file (R1)"
	exit 1
fi
```

---
## Task 2

Task 2 involved removing intermediary files generated during the script execution in order to minimise memory usage. Below is the command I used to remove what I judged to be intermediary files.
```bash
rm cols2.txt cols3.txt cols4578.txt ../mapping/my_sample.sam annovar_out.exonic_variant_function annovar_out.variant_function my_sample.pileup annovar_out.exonic_variant_function annovar_out.variant_function genotype_data.txt
```
- Extracted column text files - Should be removed as they are only used to generate final.txt and snps.csv. Data is redundant as it can be extracted again from the vcf file if needed.
- SAM file - Can be removed as we have generated BAM file which can be used for analysis and is a binary file so is more efficient and takes up less space.
- Pileup file - Pileup files are intermediary files used by samtools to create vcf file. Considering we do not need the raw file for further analysis it can be removed to save space.
- Annovar generated files - These files can be removed as any needed data has been extracted to final output files and vcf file contains all necessary annotations.

There could be an argument made for the removal of more files, so that if the user is doing another run of the script with different input files there would be no error, as when the script is rerun all necessary files are generated. As long as the archive file is conserved then the rest of the data could in theory be discarded. 
The downsides of that approach is that the user would not be able to verify the contents of any of the files post script execution for further analysis or error checking.

---
## Task 3

My approach to this task can be broken down into a series of steps:
1. Remove headers from vcf file and extract 8th column
```bash
grep -v "^#" my_sample.vcf | cut -f8 |
```
2. Split the column by ";" delimiter and extract 3rd and 4th columns.
```bash
cut -f3,4 -d';' 
```
3. Use awk to first split the remaining data into 2 data column with format: $data1=value1; $data2=value2.
```bash
awk -F';'
```
4. Within the awk command, spilt the data by the = sign and store in an array. (b for homozygous and a for heterozygous)
```bash
	split($1, a, "="); 
	split($2, b, "="); 
```
5. Use conditional statement to print the output to a text file, 1/1 for homozygous, 0/1 for heterozygous.
```bash
	if (b[2] == 1) 
			print "1/1"; 
		else if (a[2] == 1) 
			print "0/1"; 
		else 
			print "0/1";
```
6. Add genotype header to text file. (Use intermediary temp file first in case concat fails so data is preserved.)
```bash
# Add the "Genotype" header to the genotype_data.txt file
	echo "Genotype" | cat - genotype_data.txt > temp_genotype_data.txt && mv temp_genotype_data.txt genotype_data.txt
```
7. Combine genotype text file with snps.csv using paste for our final output.
```bash
# Combine the header with the existing snps.csv using paste
	paste -d',' snps.csv genotype_data.txt > temp_snps.csv && mv temp_snps.csv snps.csv
```

---

## Task 4

This task involved running varscan with different min_coverage settings and using diff to identify differences.

First I created a loop to run varscan with the different minimum coverage parameters and store them.
```bash
# Run varscan with different min coverage values and compare outputs
	for min_coverage in 5 10 15 20; do
		java -jar /software/genomics/varscan/2.4.1/VarScan.jar mpileup2snp $SAMPLE\.pileup --min-coverage $min_coverage --output-vcf 1 > $SAMPLE\_mincov$min_coverage\.vcf 2>> ../log
		echo "Done with min coverage $min_coverage"
	done
```

Then loop through all the generated vcf files and compare them to the original vcf file generated with default settings.  This first part ensures that the reference file does not compare itself to itself.

```bash
 Loop through all VCF files in the current directory

    for vcf_file in *.vcf; do

        # Skip the reference file

        if [[ "$vcf_file" == "$reference_file" ]]; then

            continue

        fi

```

Then the script strips the headers from the files and assigns each to a variable.
Then diff is run to compare the files. Using "|| true" stops the script from exiting if a difference is found.

```bash

        echo "Comparing $vcf_file to $reference_file"

        vcf_content=$(grep -v "^#" $vcf_file)

        reference_content=$(grep -v "^#" $reference_file)

        diff_ouput=""

        # Compare the current VCF file to the reference and save differences to a temp file

        diff_output=$(diff <(echo "$reference_content") <(echo "$vcf_content") || true)
 ```
 
If there is a difference they are saved to a text file.

```bash
        # Check if there are differences

        if [ -z "$diff_output" ]; then

            echo "No differences found between $vcf_file and $reference_file."

        else

            echo "Differences found between $vcf_file and $reference_file. Saving to ${vcf_file%.vcf}_diff.txt"

            echo -e "$diff_output\n" > "${vcf_file%.vcf}_diff.txt"

        fi

    done

```

A higher minimum coverage led to less SNPs as positions with fewer reads than that threshold are excluded from SNP analysis as they could be false positives.
