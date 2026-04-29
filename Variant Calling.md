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
Run the already-called variants from above (in the quasi VCF files) through this pipeline to count the number of variants:

```console
#!/bin/bash

for vcf in *_27Aug_rmdup.vcf; do
    sample=$(basename "$vcf" .vcf)
    outfile="April29-2026-ELrerun/${sample}_variant-29Apr26.csv"
    logfile="April29-2026-ELrerun/${sample}_variant-29Apr26.log"

    nohup python3 variant_updatedEL080825.py \
        --v "$vcf" \
        --snp /N/slate/lhkelley/GSF4254/sailor/c.elegans.WS275.snps.nostrand.sorted.bed \
        --o "$outfile" > "$logfile" 2>&1 &

    while [ "$(jobs -r | wc -l)" -ge "$max_jobs" ]; do
        sleep 1
    done
done
```
