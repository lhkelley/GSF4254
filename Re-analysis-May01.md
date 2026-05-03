# Alignment

Copy ```annotation.gtf``` and ```assembly.fasta``` to ```May3_EWexact```.
Keeping fastq files in ```GSF4254``` and will pull from there.

```console
conda activate rnaseq
```

Create genome index:
```console
STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir genome_index \
  --genomeFastaFiles assembly.fasta \
  --sjdbGTFfile annotation.gtf \
  --genomeSAindexNbases 12
```

Align reads (only doing N2 first):
```console
for fq in /N/slate/lhkelley/GSF4254/*N2*.fastq; do
  prefix=$(basename "$fq" .fastq)

  STAR \
    --runThreadN 8 \
    --outFilterMultimapNmax 1 \
    --outFilterScoreMinOverLread 0.66 \
    --outFilterMismatchNmax 10 \
    --outFilterMismatchNoverLmax 0.3 \
    --runMode alignReads \
    --genomeDir genome_index \
    --readFilesIn "$fq" \
    --outFileNamePrefix "${prefix}_" \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate
done
```
Index reads:
```console
for BAM in *.bam; do
    samtools index "$BAM"
done
```

# Remove duplicated reads

Use ```samtools markdup``` with the ```-r``` to remove duplicate reads. The ```-s``` flag prints statistics.
```console
for i in *.bam; do
    out="${i%.bam}.nodup.bam"
    samtools markdup -r -s "$i" "$out"
done
```

Count the number of reads that were removed and compare to original number of reads:
```console
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

# Variant calling

```mpileup``` will pile up the reads at each variant sites. The ```-a DP4``` refers counting/including the number of reads for the reference and alternative nucleotide, as well as the forward and reverse strand. So, there are 4 categories: reference/forward, reference/reverse, alternative/forward, and alternative/reverse (in that order).

```call``` is actually doing the variant calling. The ```-m``` flag is saying to use the multiallelic calling model. ```-v``` means only include variant sites, not sites where all counts match the reference. The ```-A``` flag means to keep all alternative positions, even low frequency or weakly supported variants.

```console
# Loop over all deduplicated BAM files
for bam in *.nodup.bam; do
    # Get sample name without extension
    sample=$(basename "$bam" .nodup.bam)
    outfile="${sample}_rmdup.vcf"

    # Skip if output exists
    if [ -f "$outfile" ]; then
        echo "Skipping $bam, output $outfile already exists."
        continue
    fi

    # Generate pileup and annotate DP4, ignore indels
    bcftools mpileup -Ou -f "$ref" -a DP4 "$bam" | \
    bcftools call -mv -A -Ov -o "$outfile"

    echo "Variants generated for $bam -> $outfile"
done
```
