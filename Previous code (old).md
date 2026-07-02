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
## Shorten file names

```console
for f in GSF4254-N2-rep*_Aligned.sortedByCoord.out.nodup.bam; do
    rep=$(echo "$f" | grep -o 'rep[0-9]\+')
    mv "$f" "N2_${rep}_nodup.bam"
done
```

## Pile up the reads

```console
#!/bin/bash

ref="assembly.fasta"

for bam in *_nodup.bam; do
    sample=$(basename "$bam" _nodup.bam)
    outfile="${sample}_01May.vcf"
    mpileup_file="${sample}_01May.mpileup.vcf.gz"

    if [ -f "$outfile" ]; then
        echo "Skipping $bam, output $outfile already exists."
        continue
    fi

    echo "Processing $bam"

    # Step 1: generate and save  mpileup (for troubleshooting and sanity checks)
    bcftools mpileup -f "$ref" -a DP4 "$bam" -Ou \
        | bcftools view -Oz -o "$mpileup_file"

    bcftools index "$mpileup_file"

    # Step 2: call variants from mpileup
    bcftools call -mv -A -Ov "$mpileup_file" -o "$outfile"

    echo "Variants generated for $bam -> $outfile"
done
```

```console
   
```

## Variant calling - WRONG

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

Clean up and calculate variant frequency:
```console
#!/bin/bash

# Loop over all rmdup VCF files
for vcf in *.vcf; do
    # Get sample name without .vcf
    sample=$(basename "$vcf" .vcf)

    # Define output file and log file
    outfile="${sample}_variant.csv"
    logfile="${sample}_variant.log"

    # Run the python script --> if your are including EVERYTHING, not just variants, your files will be huge and this must be submitted as a job...
    /N/slate/lhkelley/GSF4254/varientcall/variant_updatedEL080825.py \
        --v "$vcf" \
        --snp /N/slate/lhkelley/GSF4254/sailor/c.elegans.WS275.snps.nostrand.sorted.bed \
        --o "$outfile" > "$logfile" 2>&1 &

    echo "Started processing $vcf -> $outfile (logging to $logfile)"
done
```
