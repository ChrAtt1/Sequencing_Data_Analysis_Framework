#!/bin/bash

# Activate conda environment for samtools, bcftools, and whatshap
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

# Function to process reads: align, convert SAM to BAM, sort, index, and variant call
process_reads() {
    
    local file_base_name=$(basename "${file_r1/_R1/}" .fastq)
    local file_r2=${file_r1//_R1/_R2}    
    local bwa_output="$output_data_aligned/${file_base_name}.sam"
    local sam_output="$output_data_aligned/${file_base_name}.sam"
    local bam_output_base="$output_data_converted/${file_base_name}"
    local bam_output="$bam_output_base.bam"
    local bam_output_sorted="$bam_output_base_sorted.bam"
    local bcf_output="$output_data_variant_called/${file_base_name}_called.bcf"
    local bcf_output_highqual="$output_data_variant_called/${file_base_name}_called-highqual.bcf"
    local vcf_output_highqual="$output_data_variant_called/${file_base_name}_called-highqual.vcf"
    local vcf_output_highqual_gz="$vcf_output_highqual.gz"
    local whatshap_output_vcf="$output_data_phased/${file_base_name}_phased.vcf"
    local whatshap_output_vcf_gz="$whatshap_output_vcf.gz"
    local whatshap_output_gtf="$output_data_phased/${file_base_name}_called-highqual.gtf"
    local whatshap_output_stats="$output_data_phased/${file_base_name}_statistics.ods"
    local whatshap_output_haplotagged_bam="$output_data_phased/${file_base_name}_haplotagged.bam"
    local bcf_output_haplotype1_fasta="$output_data_haplotypes/${file_base_name}_haplotype1.fasta"
    local bcf_output_haplotype2_fasta="$output_data_haplotypes/${file_base_name}_haplotype2.fasta"
    
    conda activate /home/smink/anaconda3/envs/samtools
    if [ ! -f "$path_reference_sequence.bwt" ]; then
    bwa index $path_reference_sequence
    fi
    
    echo $file_r1
    echo $file_r2
    echo $bwa_output
    
    # Execute BWA mem, convert SAM to BAM, sort and index the BAM file
    bwa mem $path_reference_sequence $file_r1 $file_r2 > $bwa_output
    samtools view -S -b $bwa_output > "$bam_output"
    samtools sort "$bam_output" -o "$bam_output_sorted"
    samtools index "$bam_output_sorted"
    
    # Perform variant calling using bcftools
    bcftools mpileup --max-depth 1000 -f "$path_reference_sequence" "$bam_output_sorted" | bcftools call -P 1.1e-3 -mv -Ob -o "$bcf_output"
    bcftools view -i '%QUAL>=10' "$bcf_output" -o "$bcf_output_highqual"
    bcftools convert -Ov -o "$vcf_output_highqual" "$bcf_output_highqual"
    bgzip -c "$vcf_output_highqual" > "$vcf_output_highqual_gz"
    tabix "$vcf_output_highqual_gz"
    conda deactivate
    
    # Whatshap phase
    whatshap phase -o $whatshap_output_vcf --ignore-read-groups --reference=$path_reference_sequence $vcf_output_highqual $bam_output_sorted

    # Whatshap stats
    whatshap stats --gtf=$whatshap_output_gtf $whatshap_output_vcf
    
    conda activate /home/smink/anaconda3/envs/samtools

    # Bcftools convert to gz
    bcftools convert -Oz -o "$whatshap_output_vcf_gz" "$whatshap_output_vcf"

    # Index the gzipped VCF
    bcftools index "$whatshap_output_vcf_gz"
    
    conda deactivate

    # Whatshap haplotag
    whatshap haplotag --ignore-read-groups -o "$whatshap_output_haplotagged_bam" --reference="$path_reference_sequence" "$whatshap_output_vcf_gz" "$bam_output_sorted"
    
    conda activate /home/smink/anaconda3/envs/samtools
    
    # Samtools index the haplotagged BAM
    samtools index "$whatshap_output_haplotagged_bam"
    conda deactivate

    # Whatshap stats for phased VCF and generate statistics
    whatshap stats "$whatshap_output_vcf" --tsv="$whatshap_output_stats"
    
    conda activate /home/smink/anaconda3/envs/samtools
    # Bcftool calls
    bcftools consensus -H 1 -f $path_reference_sequence $whatshap_output_vcf_gz > $bcf_output_haplotype1_fasta
    bcftools consensus -H 2 -f $path_reference_sequence $whatshap_output_vcf_gz > $bcf_output_haplotype2_fasta
    conda deactivate
}

# Function to create directory if not exists
ensure_dir() {
    [ ! -d "$1" ] && mkdir -p "$1"
}

# Define base paths
# Specify the base path where the data is located
path_base_data="/path/to/your/base/data"
# Specify the path to the input data
input_data="/path/to/your/input/data"
# Specify the path to the output data
output_data="/path/to/your/output/data"

# Set output data paths
output_data_filtered="$output_data/1_filtered_dataset"
output_data_aligned="$output_data/2_mapped_dataset"
output_data_converted="$output_data/3_converted_dataset"
output_data_variant_called="$output_data/4_variant_called_dataset"
output_data_phased="$output_data/5_phased_dataset"
output_data_haplotypes="$output_data/6_haplotypes_dataset"

# Export reference sequence once
export path_reference_sequence

# Loop over each file for processing
ensure_dir "$output_data_filtered"
ensure_dir "$output_data_aligned"
ensure_dir "$output_data_converted"
ensure_dir "$output_data_variant_called"
ensure_dir "$output_data_phased"
ensure_dir "$output_data_haplotypes"

# Step 1 (optional): Filter data using fastq filter if the flag is true
# Check for flags since data filtering is optional
FILTER_DATA=false

# Filter data if it is set to true
if [ "$FILTER_DATA" = true ]; then
    echo "Filtering data using fastq filter"

    # Create the filtered data directory if it does not exist
    if [ ! -d "$output_data_filtered" ]; then
        mkdir -p "$output_data_filtered"
    fi

    for file in "$input_data"/*R1*.fastq; do
        if [ -e "$file" ]; then
            file_base=$(basename "$file" .fastq)
            file_R1="$output_data_filtered/${file_base}.fastq"
            file_R2="${file/R1/R2}"

            # Filter settings of fastq
            fastq-filter -l 250 -Q 30 -o "$file_R1" -o "${file_R1/R1/R2}" "$file" "$file_R2"
        else
            echo "$file does NOT exist."
        fi
    done
    
else
    # If skipping the filtering step, use the input data as filtered data for subsequent steps
    output_data_filtered="$input_data"
fi

# Process each fastq file in the input directory
for file_r1 in "$input_data"/*R1*.fastq; do
    if [ -f "$file_r1" ]; then
  
        process_reads "$file_r1"
        
    else
        echo "File does not exist: $file_r1"
    fi
done
