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

Output:
```console
samtools markdup: warning, unable to calculate estimated library size. Read pairs 0 should be greater than duplicate pairs 0, which should both be non zero.

COMMAND: samtools markdup -r -s GSF4254-N2-rep1_S10_R1_001_Aligned.sortedByCoord.out.bam GSF4254-N2-rep1_S10_R1_001_Aligned.sortedByCoord.out.nodup.bam
READ: 39755097
WRITTEN: 5025145
EXCLUDED: 0
EXAMINED: 39755097
PAIRED: 0
SINGLE: 39755097
DUPLICATE PAIR: 0
DUPLICATE SINGLE: 34729952
DUPLICATE PAIR OPTICAL: 0
DUPLICATE SINGLE OPTICAL: 0
DUPLICATE NON PRIMARY: 0
DUPLICATE NON PRIMARY OPTICAL: 0
DUPLICATE PRIMARY TOTAL: 34729952
DUPLICATE TOTAL: 34729952
ESTIMATED_LIBRARY_SIZE: 0

samtools markdup: warning, unable to calculate estimated library size. Read pairs 0 should be greater than duplicate pairs 0, which should both be non zero.

COMMAND: samtools markdup -r -s GSF4254-N2-rep2_S11_R1_001_Aligned.sortedByCoord.out.bam GSF4254-N2-rep2_S11_R1_001_Aligned.sortedByCoord.out.nodup.bam
READ: 32660880
WRITTEN: 4360273
EXCLUDED: 0
EXAMINED: 32660880
PAIRED: 0
SINGLE: 32660880
DUPLICATE PAIR: 0
DUPLICATE SINGLE: 28300607
DUPLICATE PAIR OPTICAL: 0
DUPLICATE SINGLE OPTICAL: 0
DUPLICATE NON PRIMARY: 0
DUPLICATE NON PRIMARY OPTICAL: 0
DUPLICATE PRIMARY TOTAL: 28300607
DUPLICATE TOTAL: 28300607
ESTIMATED_LIBRARY_SIZE: 0

samtools markdup: warning, unable to calculate estimated library size. Read pairs 0 should be greater than duplicate pairs 0, which should both be non zero.

COMMAND: samtools markdup -r -s GSF4254-N2-rep3_S12_R1_001_Aligned.sortedByCoord.out.bam GSF4254-N2-rep3_S12_R1_001_Aligned.sortedByCoord.out.nodup.bam
READ: 84246559
WRITTEN: 6889472
EXCLUDED: 0
EXAMINED: 84246559
PAIRED: 0
SINGLE: 84246559
DUPLICATE PAIR: 0
DUPLICATE SINGLE: 77357087
DUPLICATE PAIR OPTICAL: 0
DUPLICATE SINGLE OPTICAL: 0
DUPLICATE NON PRIMARY: 0
DUPLICATE NON PRIMARY OPTICAL: 0
DUPLICATE PRIMARY TOTAL: 77357087
DUPLICATE TOTAL: 77357087
ESTIMATED_LIBRARY_SIZE: 0
```


Count the number of reads that were removed and compare to original number of reads:
```console
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

Output:
```console
(rnaseq) [lhkelley@i73 May3_EWexact]$ cat read_counts.txt 
Sample	Mapped_Reads_Before	Mapped_Reads_After	Pct_Remaining
GSF4254-N2-rep1_S10_R1_001_Aligned.sortedByCoord.out.bam	39755097	5025145	12.64
GSF4254-N2-rep2_S11_R1_001_Aligned.sortedByCoord.out.bam	32660880	4360273	13.35
GSF4254-N2-rep3_S12_R1_001_Aligned.sortedByCoord.out.bam	84246559	6889472	8.18
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
