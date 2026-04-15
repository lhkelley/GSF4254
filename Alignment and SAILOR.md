## Align reads

Create genomic index:

```console
STAR \
  --runThreadN 4 \  # Number of CPU threads used in parallel
  --runMode genomeGenerate \  # Buuld a genome index (as opposed to aligning reads)
  --genomeDir genome/index \  # Directory where the genome index will be stored (create this before performing this command)
  --genomeFastaFiles genome/assembly.fasta \  # Reference genome FASTA file that you provide
  --sjdbGTFfile genome/annotation.gtf \  # Annotation file that you proivde; helps STAR handle exon-exon junctions
  --genomeSAindexNbases 12  # Length of the suffix array index; affects memory usage and speed
```

Align reads to genome:

```console

FASTQ=$GSF4254*

for FASTQ in ${FASTQ[@]}; do
  PREFIX=results/aligned/$(basename $FASTQ .fastq)_
  STAR \
    --runThreadN 8 \  # Uses 8 threads in parallel
    --outFilterMultimapNmax 1 \  # Set the maximum number of loci that a read is allowed to map to
    --outFilterScoreMinOverLread .66 \  # Set the minimum alignment score, scaled by read length
    --outFilterMismatchNmax 10 \  # Set the maximum number of mismatches allowed per read
    #--outFilterMismatchNoverLmax .3 \  # Set the maximum number of mismatches allowed per read, scaled by read length. I'm omitting it because Nmax supersedes it.
    --runMode alignReads \
    --genomeDir genome/index \
    --readFilesIn $FASTQ \
    --outFileNamePrefix $PREFIX \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate
done
```
## Decide if you should merge your BAM files or not

I will merge the replicates because I want to identify any possible site.

```console
for sample in GSF4254-709 GSF4254-712 GSF4254-adr2 GSF4254-N2
do
    # Collect all replicate BAMs for this sample
    bam_files=$(ls ${sample}-rep*_Aligned.sortedByCoord.out.bam)

    # Define output file
    outbam="${sample}_merged.bam"

    echo "Merging BAMs for ${sample} -> ${outbam}"

    # Merge with samtools
    samtools merge -@ 8 -o "$outbam" $bam_files

    # Index merged BAM
    samtools index "$outbam"
done
```

## Identify edit sites with SAILOR

### Wildtype

```run-sailor-N2.json``` :

```console
{
  "samples_path":"/N/slate/lhkelley/GSF4254/sailor/",
  "samples": [
    "GSF4254-N2.merged.bam",
  ],
  "reverse_stranded":true,
  "reference_fasta": "/N/slate/lhkelley/GSF4254/genome/assembly.fasta",
  "known_snps": "/N/slate/lhkelley/GSF4254/sailor/c.elegans.WS275.snps.nostrand.sorted.bed",
  "edit_type": "AG",
  "output_dir": "/N/slate/lhkelley/GSF4254/results/fromSAILOR_N2results"
}
```

Run: 

```console
snakemake \
  --snakefile /N/slate/lhkelley/GSF4254/sailor/workflow_sailor/Snakefile \
  --configfile /N/slate/lhkelley/GSF4254/sailor/run-sailor-N2.json \
  --use-singularity \
  --singularity-args "\ 
    --bind /N/slate/lhkelley/GSF4254 \
    --bind /N/slate/lhkelley/GSF4254/genome \
    --bind /N/slate/lhkelley/GSF4254/sailor/workflow_sailor/scripts \
    --bind /N/slate/lhkelley/GSF4254/sailor/workflow_sailor \
    --bind /N/slate/lhkelley/GSF4254/ranSailor_N2" \
  -j1
```

Run annotation:

```console
python3 annotator.sailor.py \
  --gtf c_elegans.PRJNA13758.WS275.canonical_geneset.gtf \
  --fwd /N/slate/lhkelley/GSF4254/results/fromSAILOR_N2results/7_scored_outputs/GSF4254-N2.merged.bam.fwd.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed \
  --rev /N/slate/lhkelley/GSF4254/results/fromSAILOR_N2results/7_scored_outputs/GSF4254-N2.merged.bam.rev.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed \
  --wb c.elegans.WS275.annotation.final.bed \
  --o /N/slate/lhkelley/GSF4254/results/fromSAILOR_N2results/N2.merged.FLAREannotated.sites.csv
```

```console
scp lhkelley@quartz.uits.iu.edu:/N/slate/lhkelley/GSF4254/results/fromSAILOR_N2results/N2.merged.FLAREannotated.sites.csv ~/Desktop/
```

### *adr-2(-)*

```run-sailor-adr2.json``` :

```console
{
  "samples_path":"/N/slate/lhkelley/GSF4254/sailor/",
  "samples": [
    "GSF4254-adr2.merged.bam",
  ],
  "reverse_stranded":true,
  "reference_fasta": "/N/slate/lhkelley/GSF4254/genome/assembly.fasta",
  "known_snps": "/N/slate/lhkelley/GSF4254/sailor/c.elegans.WS275.snps.nostrand.sorted.bed",
  "edit_type": "AG",
  "output_dir": "/N/slate/lhkelley/GSF4254/results/fromSAILOR"
}
```

Run: 

```console
snakemake \
  --snakefile /N/slate/lhkelley/GSF4254/sailor/workflow_sailor/Snakefile \
  --configfile /N/slate/lhkelley/GSF4254/sailor/run-sailor-adr2.json \
  --use-singularity \
  --singularity-args "\ 
    --bind /N/slate/lhkelley/GSF4254 \
    --bind /N/slate/lhkelley/GSF4254/genome \
    --bind /N/slate/lhkelley/GSF4254/sailor/workflow_sailor/scripts \
    --bind /N/slate/lhkelley/GSF4254/sailor/workflow_sailor \
    --bind /N/slate/lhkelley/GSF4254/ranSailor_adr2" \
  -j1
```

Run annotation:

```console
python3 annotator.sailor.py \
  --gtf c_elegans.PRJNA13758.WS275.canonical_geneset.gtf \
  --fwd /N/slate/lhkelley/GSF4254/results/fromSAILOR/7_scored_outputs/GSF4254-adr2.merged.bam.fwd.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed \
  --rev /N/slate/lhkelley/GSF4254/results/fromSAILOR/7_scored_outputs/GSF4254-adr2.merged.bam.rev.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed \
  --wb c.elegans.WS275.annotation.final.bed \
  --o /N/slate/lhkelley/GSF4254/results/fromSAILOR/adr2.merged.FLAREannotated.sites.csv
```

```console
scp lhkelley@quartz.uits.iu.edu:/N/slate/lhkelley/GSF4254/results/fromSAILOR/adr2.merged.FLAREannotated.sites.csv ~/Desktop/
```
