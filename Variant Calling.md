## Remove duplicates from BAM files

```console
for i in *.bam; do
    out="${i%.bam}.nodup.bam"
    samtools markdup -r -s "$i" "$out"
done
```
Compare the number of reads before and after removing duplicates. The ```pct_Remaining``` is the number of non-duplicated reads in the new BAM file.

```console
# Create header
echo -e "Sample\tMapped_Reads_Before\tMapped_Reads_After\tPct_Remaining" > read_counts.txt

# Loop through all original BAMs
for bam in *.bam; do
    # skip files that are already .nodup.bam
    if [[ "$bam" == *.nodup.bam ]]; then
        continue
    fi

    # Count mapped reads before
    before=$(samtools view -c -F 4 "$bam")

    # Find corresponding deduplicated BAM
    after_bam="${bam%.bam}.nodup.bam"
    if [[ -f "$after_bam" ]]; then
        after=$(samtools view -c -F 4 "$after_bam")
    else
        after="NA"
    fi

    # Calculate percentage remaining (skip if after is NA)
    if [[ "$after" != "NA" && "$before" -ne 0 ]]; then
        pct=$(awk -v a="$after" -v b="$before" 'BEGIN {printf "%.2f", (a/b)*100}')
    else
        pct="NA"
    fi

    # Append to file
    echo -e "$(basename $bam)\t$before\t$after\t$pct" >> read_counts.txt
done

```
## Pile up the reads

With newer versions of ```samtools```, the ```samtools mpileup``` from EE's GSF3231 analysis does not use the same command flag and it's recommended to use ```bcftools mpileup```. Below is the equivalent of the old ```samtools``` command.

Note: This is not really making a proper VCF file, it's only piling up the reads, but not yet calling variants.

```console
#!/bin/bash

# Reference genome
ref="assembly.fasta"

# Loop over all deduplicated BAM files
for bam in *.nodup.bam; do
    # Get sample name without extension
    sample=$(basename "$bam" .nodup.bam)
    outfile="${sample}_27Aug_rmdup.vcf"

    # Skip if output exists
    if [ -f "$outfile" ]; then
        echo "Skipping $bam, output $outfile already exists."
        continue
    fi

    # Generate pileup and annotate DP4, ignore indels
    bcftools mpileup -f "$ref" -a DP4 -I "$bam" | \
    bcftools call -mv -Ov -o "$outfile"

    echo "Variants generated for $bam -> $outfile"
done
```

```console
#!/bin/bash

# Loop over all rmdup VCF files
for vcf in *_27Aug_rmdup.vcf; do
    # Get sample name without .vcf
    sample=$(basename "$vcf" .vcf)

    # Define output file and log file
    outfile="${sample}_variant.csv"
    logfile="${sample}_variant.log"

    # Run the python script with nohup in background
    nohup python variant_updatedLHK073025.py \
        --v "$vcf" \
        --snp c.elegans.WS275.snps.sorted.bed \
        --o "$outfile" > "$logfile" 2>&1 &

    echo "Started processing $vcf -> $outfile (logging to $logfile)"
done
```
