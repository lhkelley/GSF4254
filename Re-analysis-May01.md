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

Align reads:
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

