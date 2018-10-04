This document contains all information required to reproduce the results in the paper. File paths have been shortened to filenames only for brevity.

# Software versions

Here are the versions and appropriate links to all software used for the papers benchmarks.

#### Variant callers

* DeepVariant 0.5.1
* FreeBayes: v1.1.0-56-ga180635
* GATK: v4.0.0.0
* Octopus: v0.4-beta
* Platypus: 0.8.1
* Strelka2: v2.9.1
* LoFreq: v2.1.3.1
* VarDict: 1.5.2 (java)
* Lancet: v1.0.6

#### Other tools

* BWA: 0.7.17-r1188
* Samtools: 1.7 (using htslib 1.7)
* Bcftools: 1.7 (using htslib 1.7)
* RTG Tools: 3.9.1


#### Data

* GIAB v3.3.2
* GATK resource budle

# Germline WGS

## Data

#### Platinum genomes NA12878

We downloaded raw FASTQ files from Google cloud:

```shell
$ wget https://storage.googleapis.com/genomics-public-data/platinum-genomes/fastq/ERR194147_1.fastq.gz
$ wget https://storage.googleapis.com/genomics-public-data/platinum-genomes/fastq/ERR194147_2.fastq.gz
```

#### Precision FDA truth NA12878

Paired FASTQs downloaded from private urls. Login required.

#### Precision FDA truth HG002

Paired FASTQs downloaded from private urls. Login required.

#### Precision FDA consistency NA12878

Paired FASTQs downloaded from private urls. Login required.

#### GIAB HG002 10X

```shell
$ wget ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/10XGenomics/NA24385_phased_possorted_bam.bam
```

#### GIAB HG005 HiSeq

We downloaded two libraries of the HiSeq 300x run:

```shell
$ wget -r ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_250bps_fastqs/150430_HG005_Homogeneity_03_FCB-22310288
$ wget -r ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_250bps_fastqs/150506_HG005_Homogeneity_04_FCA-22365346
```

We then merged the FASTQs in each library, e.g:

```shell
$ find . -name \*_R1_*.fastq.gz -print0 | sort -z | xargs -r0 cat > 150506_HG005_Homogeneity_04_FCA-22365346_cat_R1.fastq.gz
$ find . -name \*_R2_*.fastq.gz -print0 | sort -z | xargs -r0 cat > 150506_HG005_Homogeneity_04_FCA-22365346_cat_R2.fastq.gz
```

Each produces reads of around 15x.

## Commands

Here we list the commands used to call and evaluate variants for each variant caller. We just show the commands for the platinum genomes read set, since the commands for all the others are identical (other than obvious file name changes).

#### BWA

All tools require pre-mapped reads. We used BWA-MEM:

```shell
$ bwa mem -t 15 -R "@RG\tID:1\tSM:NA12878\tLB:platinum\tPU:illumina" \
      hs37d5.fa \
      ERR194147_1.fastq.gz ERR194147_2.fastq.gz \
      | samtools view -bh | samtools sort -@ 5 -o NA12878.platinum.b37.bwa-mem.bam -
$ samtools index NA12878.platinum.b37.bwa-mem.bam
```

#### DeepVariant

Set environment variables:

```shell
$ BUCKET="gs://deepvariant"
$ BIN_VERSION="0.5.1"
$ MODEL_VERSION="0.5.0"
$ MODEL_CL="182548131"
$ BIN_BUCKET="${BUCKET}/binaries/DeepVariant/${BIN_VERSION}/DeepVariant-${BIN_VERSION}+cl-*"
$ MODEL_BUCKET="${BUCKET}/models/DeepVariant/${MODEL_VERSION}/DeepVariant-inception_v3-${MODEL_VERSION}+cl-${MODEL_CL}.data-wgs_standard"

$ BIN_DIR="deepvariant"
$ REF=hs37d5.fa
$ BAM=NA12878.platinum.b37.bwa-mem.bam
$ N_SHARDS=50
$ EXAMPLES="NA12878.platinum.b37.bwa-mem.tfrecord@${N_SHARDS}.gz"
$ GVCF_TFRECORDS="NA12878.platinum.b37.bwa-mem.gvcf.tfrecord@${N_SHARDS}.gz"
$ CALL_VARIANTS_OUTPUT="NA12878.platinum.b37.bwa-mem.cvo.tfrecord.gz"
$ OUTPUT_VCF="$deepvariant.NA12878.platinum.b37.bwa-mem.vcf.gz"
$ LOG_DIR=log
```

Make examples:

```shell
$ ( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog "${LOG_DIR}/log" --res "${LOG_DIR}" \
    python "${BIN_DIR}"/make_examples.zip \
      --mode calling \
      --ref "${REF}" \
      --reads "${BAM}" \
      --examples "${EXAMPLES}" \
      --gvcf "${GVCF_TFRECORDS}" \
      --task {}
) >"${LOG_DIR}/make_examples.log" 2>&1
```

Call variants:

```shell
$ ( time python "${BIN_DIR}"/call_variants.zip \
    --outfile "${CALL_VARIANTS_OUTPUT}" \
    --examples "${EXAMPLES}" \
    --checkpoint "${MODEL}" \
    --batch_size 32
) >"${LOG_DIR}/call_variants.log" 2>&1
```

Postprocess variants:

```shell
$ ( time python "${BIN_DIR}"/postprocess_variants.zip \
    --ref "${REF}" \
    --infile "${CALL_VARIANTS_OUTPUT}" \
    --outfile "${OUTPUT_VCF}"
) >"${LOG_DIR}/postprocess_variants.log" 2>&1
```

Evaluate variants:

```shell
$ rtg vcfeval \
    -t hs37d5_SDF \
    -b HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
    --evaluation-regions HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
    --ref-overlap \
    -c deepvariant.NA12878.platinum.b37.bwa-mem.vcf.gz \
    -o deepvariant.NA12878.platinum.b37.bwa-mem.eval
```

#### FreeBayes

Mark duplicates:

```shell
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12878.platinum.b37.bwa-mem.bam \
     -O NA12878.platinum.b37.bwa-mem.dedup.bam \
     -M NA12878.platinum.b37.bwa-mem.dedup.metrics
```

Call variants:

```shell
$ freebayes \
    -f hs37d5.fa \
    -b NA12878.platinum.b37.bwa-mem.dedup.bam \
    -t hs37d5.chromosomes.bed \
    -= | bgzip > freebayes.NA12878.platinum.b37.bwa-mem.dedup.vcf.gz
$ tabix freebayes.NA12878.platinum.b37.bwa-mem.dedup.vcf.gz
```

Filter variants:

```shell
$ bcftools filter -s FAIL \
    -i 'QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' \
    -Oz -o freebayes.NA12878.platinum.b37.bwa-mem.dedup.filtered.vcf.gz \
    freebayes.NA12878.platinum.b37.bwa-mem.dedup.vcf.gz
$ tabix freebayes.NA12878.platinum.b37.bwa-mem.dedup.vcf.gz
```

Evaluate variants:

```shell
$ rtg vcfeval \
    -t hs37d5_SDF \
    -b HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
    --evaluation-regions HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
    --ref-overlap \
    -c freebayes.NA12878.platinum.b37.bwa-mem.dedup.filtered.vcf.gz \
    -o freebayes.NA12878.platinum.b37.bwa-mem.dedup.filtered.eval
```

#### GATK HaplotypeCaller

Mark duplicates:

```shell
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12878.platinum.b37.bwa-mem.bam \
     -O NA12878.platinum.b37.bwa-mem.dedup.bam \
     -M NA12878.platinum.b37.bwa-mem.dedup.metrics
```

Recalibrate base qualities:

```shell
$ gatk --java-options -Xmx4G BaseRecalibrator \
    -R hs37d5.fa \
    -I NA12878.platinum.b37.bwa-mem.dedup.bam \
    --known-sites gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    --known-sites gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --known-sites gatk_bundle/b37/1000G_phase1.indels.b37.vcf.gz \
    -O NA12878.platinum.b37.bwa-mem.dedup.recal.table
$ gatk --java-options -Xmx4G ApplyBQSR
    -R hs37d5.fa \
    -I NA12878.platinum.b37.bwa-mem.dedup.bam \
    --bqsr-recal-file NA12878.platinum.b37.bwa-mem.dedup.recal.table \
    -O NA12878.platinum.b37.bwa-mem.dedup.recal.bam
```

Call variants:

```shell
$ gatk --java-options -Xmx12G HaplotypeCaller \
    -R hs37d5.fa \
    -I NA12878.platinum.b37.bwa-mem.dedup.recal.bam \
    -O gatk.NA12878.platinum.b37.bwa-mem.dedup.recal.vcf.gz \
    -L hs37d5.chromosomes.bed 
```

Evaluate variants:

```shell
$ rtg vcfeval \
    -t hs37d5_SDF \
    -b HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
    --evaluation-regions HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
    --ref-overlap \
    -c gatk.NA12878.platinum.b37.bwa-mem.dedup.recal.vcf.gz \
    -f QUAL \
    -o gatk.NA12878.platinum.b37.bwa-mem.dedup.recal.eval
```

#### Octopus

Call variants:

```shell
$ octopus \
    -R hs37d5.fa \
    -I NA12878.platinum.b37.bwa-mem.bam \
    -t hs37d5.chromosomes.bed \
    -o octopus.NA12878.platinum.b37.bwa-mem.vcf.gz \
    --forest octopus-germline.forest \
    --legacy
```

Evaluate variants:

```shell
$ rtg vcfeval \
    -t hs37d5_SDF \
    -b HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
    --evaluation-regions HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
    --ref-overlap \
    -c octopus.NA12878.platinum.b37.bwa-mem.legacy.vcf.gz \
    -f INFO.RFQUAL \
    -o octopus.NA12878.platinum.b37.bwa-mem.eval
```

#### Platypus

Call variants:

```shell
$ python Platypus.py callVariants \
    --refFile hs37d5.fa \
    --bamFiles NA12878.platinum.b37.bwa-mem.bam \
    --regions hs37d5.chromosomes.bed \
    --output platypus.NA12878.platinum.b37.bwa-mem.vcf
$ bgzip platypus.NA12878.platinum.b37.bwa-mem.vcf
$ tabix platypus.NA12878.platinum.b37.bwa-mem.vcf.gz
```

Evaluate variants:

```shell
$ rtg vcfeval \
    -t hs37d5_SDF \
    -b HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
    --evaluation-regions HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
    --ref-overlap \
    -c platypus.NA12878.platinum.b37.bwa-mem.vcf.gz \
    -o platypus.NA12878.platinum.b37.bwa-mem.eval
```

#### Strelka2

Call variants:

```shell
$ configureStrelkaGermlineWorkflow.py \
    --referenceFasta hs37d5.fa \
    --bam NA12878.platinum.b37.bwa-mem.bam \
    --callRegions hs37d5.chromosomes.bed.gz \
    --runDir strelka_tmp
$ strelka_tmp/runWorkflow.py -m local -j 10
$ mv strelka_tmp/results/variants/variants.vcf.gz strelka2.NA12878.platinum.b37.bwa-mem.vcf.gz
$ mv strelka_tmp/results/variants/variants.vcf.gz.tbi strelka2.NA12878.platinum.b37.bwa-mem.vcf.gz.tbi
```

Evaluate variants:

```shell
$ rtg vcfeval \
    -t hs37d5_SDF \
    -b HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
    --evaluation-regions HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
    --ref-overlap \
    -c strelka2.NA12878.platinum.b37.bwa-mem.vcf.gz \
    -f FORMAT.GQX \
    -o strelka2.NA12878.platinum.b37.bwa-mem.eval
```

# Platinum genomes trio

## Commands

#### DeepVariant

Make gvcfs compatible with GATK4:

```shell
$ bioawk -tHc vcf '{gsub("<\\*>","<NON_REF>",$alt); print}' \
    deepvariant.NA12878.platinum.b37.bwa-mem.g.vcf.gz \
    | bgzip > deepvariant.NA12878.platinum.b37.bwa-mem.non_ref.g.vcf.gz
$ tabix deepvariant.NA12878.platinum.b37.bwa-mem.non_ref.g.vcf.gz
$ bioawk -tHc vcf '{gsub("<\\*>","<NON_REF>",$alt); print}' \
    deepvariant. NA12891.platinum.b37.bwa-mem.g.vcf.gz \
    | bgzip > deepvariant. NA12891.platinum.b37.bwa-mem.non_ref.g.vcf.gz
$ tabix deepvariant.NA12891.platinum.b37.bwa-mem.non_ref.g.vcf.gz
$ bioawk -tHc vcf '{gsub("<\\*>","<NON_REF>",$alt); print}' \
    deepvariant. NA12892.platinum.b37.bwa-mem.g.vcf.gz \
    | bgzip > deepvariant. NA12892.platinum.b37.bwa-mem.non_ref.g.vcf.gz
$ tabix deepvariant.NA12892.platinum.b37.bwa-mem.non_ref.g.vcf.gz
```

Combine gvcfs:

```shell
$ gatk --java-options -Xmx12G CombineGVCFs \
    -R hs37d5.fa \
    -V deepvariant.NA12878.platinum.b37.bwa-mem.non_ref.g.vcf.gz \
    -V deepvariant.NA12891.platinum.b37.bwa-mem.non_ref.g.vcf.gz \
    -V deepvariant.NA12892.platinum.b37.bwa-mem.non_ref.g.vcf.gz \
    -L hs37d5.chromosomes.bed \
    -O deepvariant.ceu-trio.platinum.b37.bwa-mem.g.vcf
```

Genotype gvcf:

```shell
$ gatk --java-options -Xmx12G GenotypeGVCFs \
    -R hs37d5.fa \
    -V deepvariant.ceu-trio.platinum.b37.bwa-mem.g.vcf \
    -L hs37d5.chromosomes.bed \
    -D gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    -O deepvariant.ceu-trio.platinum.b37.bwa-mem.vcf.gz
```

Annotate *de novo* mutations:

**Note:** This operation was not possible with version 4.0.0.0 (or the later 4.0.5.1), so we had to revert to version 3.8.1.

```shell
$ java -jar GenomeAnalysisTK.jar -T VariantAnnotator \
    -R hs37d5.fa \
    -V deepvariant.ceu-trio.platinum.b37.bwa-mem.vcf.gz \
    -A PossibleDeNovo \
    -ped ceu-trio.ped \
    -o deepvariant.ceu-trio.platinum.b37.bwa-mem.denovo_annotated.vcf.gz
```

Filter *de novos*:

```shell
$ bcftools view \
    -i 'INFO/hiConfDeNovo="NA12878"' \
    -Oz -o deepvariant.ceu-trio.platinum.b37.bwa-mem.denovo.vcf.gz \
    deepvariant.ceu-trio.platinum.b37.bwa-mem.denovo_annotated.vcf.gz
```

#### GATK HaplotypeCaller

Mark duplicates:

```shell
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12878.platinum.b37.bwa-mem.bam \
     -O NA12878.platinum.b37.bwa-mem.dedup.bam \
     -M NA12878.platinum.b37.bwa-mem.dedup.metrics
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12891.platinum.b37.bwa-mem.bam \
     -O NA12891.platinum.b37.bwa-mem.dedup.bam \
     -M NA12891.platinum.b37.bwa-mem.dedup.metrics
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12892.platinum.b37.bwa-mem.bam \
     -O NA12892.platinum.b37.bwa-mem.dedup.bam \
     -M NA12892.platinum.b37.bwa-mem.dedup.metrics
```

Recalibrate base qualities:

```shell
$ gatk --java-options -Xmx4G BaseRecalibrator \
    -R hs37d5.fa \
    -I NA12878.platinum.b37.bwa-mem.dedup.bam \
    --known-sites gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    --known-sites gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --known-sites gatk_bundle/b37/1000G_phase1.indels.b37.vcf.gz \
    -O NA12878.platinum.b37.bwa-mem.dedup.recal.table
$ gatk --java-options -Xmx4G ApplyBQSR
    -R hs37d5.fa \
    -I NA12878.platinum.b37.bwa-mem.dedup.bam \
    --bqsr-recal-file NA12878.platinum.b37.bwa-mem.dedup.recal.table \
    -O NA12878.platinum.b37.bwa-mem.dedup.recal.bam
$ gatk --java-options -Xmx4G BaseRecalibrator \
    -R hs37d5.fa \
    -I NA12891.platinum.b37.bwa-mem.dedup.bam \
    --known-sites gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    --known-sites gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --NA12891-sites gatk_bundle/b37/1000G_phase1.indels.b37.vcf.gz \
    -O NA12878.platinum.b37.bwa-mem.dedup.recal.table
$ gatk --java-options -Xmx4G ApplyBQSR
    -R hs37d5.fa \
    -I NA12891.platinum.b37.bwa-mem.dedup.bam \
    --bqsr-recal-file NA12878.platinum.b37.bwa-mem.dedup.recal.table \
    -O NA12891.platinum.b37.bwa-mem.dedup.recal.bam
$ gatk --java-options -Xmx4G BaseRecalibrator \
    -R hs37d5.fa \
    -I NA12892.platinum.b37.bwa-mem.dedup.bam \
    --known-sites gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    --known-sites gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --known-sites gatk_bundle/b37/1000G_phase1.indels.b37.vcf.gz \
    -O NA12892.platinum.b37.bwa-mem.dedup.recal.table
$ gatk --java-options -Xmx4G ApplyBQSR
    -R hs37d5.fa \
    -I NA12892.platinum.b37.bwa-mem.dedup.bam \
    --bqsr-recal-file NA12878.platinum.b37.bwa-mem.dedup.recal.table \
    -O NA12892.platinum.b37.bwa-mem.dedup.recal.bam
```

Make gvcfs:

```shell
$ gatk --java-options -Xmx12G HaplotypeCaller \
    -R hs37d5.fa \
    -I NA12878.platinum.b37.bwa-mem.dedup.recal.bam \
    -L hs37d5.chromosomes.bed \
    -ERC GVCF \
    -O gatk.NA12878.platinum.b37.bwa-mem.dedup.recal.g.vcf
$ gatk --java-options -Xmx12G HaplotypeCaller \
    -R hs37d5.fa \
    -I NA12891.platinum.b37.bwa-mem.dedup.recal.bam \
    -L hs37d5.chromosomes.bed \
    -ERC GVCF \
    -O NA12891.platinum.b37.bwa-mem.dedup.recal.g.vcf
$ gatk --java-options -Xmx12G HaplotypeCaller \
    -R hs37d5.fa \
    -I NA12892.platinum.b37.bwa-mem.dedup.recal.bam \
    -L hs37d5.chromosomes.bed \
    -ERC GVCF \
    -O gatk.NA12892.platinum.b37.bwa-mem.dedup.recal.g.vcf
```

Combine gvcfs:

```shell
$ gatk --java-options -Xmx12G CombineGVCFs \
    -R hs37d5.fa \
    -V gatk.NA12878.platinum.b37.bwa-mem.dedup.recal.g.vcf \
    -V gatk.NA12891.platinum.b37.bwa-mem.dedup.recal.g.vcf \
    -V gatk.NA12892.platinum.b37.bwa-mem.dedup.recal.g.vcf \
    -L hs37d5.chromosomes.bed \
    -O gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.g.vcf
```

Genotype gvcf:

```shell
$ gatk --java-options -Xmx12G GenotypeGVCFs \
    -R hs37d5.fa \
    -V gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.g.vcf \
    -L hs37d5.chromosomes.bed \
    -D gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    -O gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vcf.gz
```

Run VQSR:

```shell
$ gatk VariantRecalibrator \
    -R hs37d5.fa \
    -V gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vcf.gz \
    --resource hapmap,known=false,training=true,truth=true,prior=15.0:gatk_bundle/b37/hapmap_3.3.b37.vcf.gz \
    --resource omni,known=false,training=true,truth=false,prior=12.0:gatk_bundle/b37/1000G_omni2.5.b37.vcf.gz \
    --resource 1000G,known=false,training=true,truth=false,prior=10.0:gatk_bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -mode SNP \
    --max-attempts 4 \
    -tranche 99.5 \
    --tranches-file gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.snp.tranches \
    -L hs37d5.chromosomes.bed \
    -O gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.recal 
$ gatk ApplyVQSR \
    -R hs37d5.fa \
    -V gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vcf.gz \
    -ts-filter-level 99.5 \
    --tranches-file gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.snp.tranches \
    --recal-file gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.recal \
    -mode SNP \
    -L hs37d5.chromosomes.bed \
    -O gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.recal_snv.raw_indels.vcf.gz
$ gatk VariantRecalibrator \
    -R hs37d5.fa \
    -V gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.recal_snv.raw_indels.vcf.gz \
    --resource mills,known=false,training=true,truth=true,prior=12.0:gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
    -mode INDEL \
    --max-attempts 4 \
    --max-gaussians 4 \
    -tranche 99.0 \
    --tranches-file gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.recal_snv.raw_indels.trances \
    -L hs37d5.chromosomes.bed \
    -O gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.recal_snv.raw_indels.indel.recal
$ gatk ApplyVQSR \
    -R hs37d5.fa \
    -V gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.recal_snv.raw_indels.vcf.gz \
    -ts-filter-level 99.0 \
    --tranches-file gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.recal_snv.raw_indels.tranches \
    --recal-file /gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.recal_snv.raw_indels.indel.recal \
    -mode INDEL \
    -L hs37d5.chromosomes.bed \
    -O gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.vcf.gz
```

Calculate genotype posteriors:

**Note:** This command did not output compressed VCF when requested. It also output unsorted VCF. We therefore had to manually sort and compress the output.

```shell
$ gatk CalculateGenotypePosteriors \
    -V gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.vcf.gz \
    -ped ceu-trio.ped \
    --de-novo-prior 1.3e-8 \
    -O gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.vcf    
$ bcftools sort \
    -Oz -o gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.sorted.vcf.gz \
    gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.vcf
$ mv gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.sorted.vcf.gz \
    gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.vcf.gz
$ tabix gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.vcf.gz                 
```

Annotate *de novo* mutations:

**Note:** This operation was not possible with version 4.0.0.0 (or the later 4.0.5.1), so we had to revert to version 3.8.1.

```shell
$ java -jar GenomeAnalysisTK.jar -T VariantAnnotator \
    -R hs37d5.fa \
    -V gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.vcf.gz \
    -A PossibleDeNovo \
    -ped ceu-trio.ped \
    -o gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.denovo_annotated.vcf.gz
```

Filter *de novos*:

```shell
$ bcftools view \
    -i 'INFO/hiConfDeNovo="NA12878"' \
    -Oz -o gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.denovo.vcf.gz \
    gatk.ceu-trio.platinum.b37.bwa-mem.dedup.recal.vqsr.cgp.denovo_annotated.vcf.gz
```

#### FreeBayes

Mark duplicates:

```shell
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12878.platinum.b37.bwa-mem.bam \
     -O NA12878.platinum.b37.bwa-mem.dedup.bam \
     -M NA12878.platinum.b37.bwa-mem.dedup.metrics
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12891.platinum.b37.bwa-mem.bam \
     -O NA12891.platinum.b37.bwa-mem.dedup.bam \
     -M NA12891.platinum.b37.bwa-mem.dedup.metrics
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12892.platinum.b37.bwa-mem.bam \
     -O NA12892.platinum.b37.bwa-mem.dedup.bam \
     -M NA12892.platinum.b37.bwa-mem.dedup.metrics
```

Call variants:

```shell
$ freebayes \
    -f hs37d5.fa \
    -b NA12878.platinum.b37.bwa-mem.dedup.bam \
       NA12891.platinum.b37.bwa-mem.bam \
       NA12892.platinum.b37.bwa-mem.bam \
    -t hs37d5.chromosomes.bed \
    -= | bgzip > freebayes.ceu-trio.platinum.b37.bwa-mem.dedup.vcf.gz
$ tabix freebayes.ceu-trio.platinum.b37.bwa-mem.dedup.vcf.gz
```

Filter variants:

```shell
$ bcftools filter -s FAIL \
    -i 'QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' \
    -Oz -o freebayes.ceu-trio.platinum.b37.bwa-mem.dedup.filtered.vcf.gz \
    freebayes.ceu-trio.platinum.b37.bwa-mem.dedup.vcf.gz
$ tabix freebayes.ceu-trio.platinum.b37.bwa-mem.dedup.vcf.gz
```

Annotate *de novo* mutations:

```shell
$ java -jar GenomeAnalysisTK.jar -T VariantAnnotator \
    -R hs37d5.fa \
    -V freebayes.ceu-trio.platinum.b37.bwa-mem.dedup.vcf.gz \
    -A PossibleDeNovo \
    -ped ceu-trio.ped \
    -o freebayes.ceu-trio.platinum.b37.bwa-mem.dedup.denovo_annotated.vcf.gz
```

Filter *de novos*:

```shell
$ bcftools view \
    -i 'INFO/hiConfDeNovo="NA12878"' \
    -Oz -o freebayes.ceu-trio.platinum.b37.bwa-mem.dedup.denovo.vcf.gz \
    freebayes.ceu-trio.platinum.b37.bwa-mem.dedup.denovo_annotated.vcf.gz
```

#### Octopus

Call *de novos*:

```shell
$ octopus
    -R hs37d5.fa \
    -I NA12878.platinum.b37.bwa-mem.bam \
       NA12891.platinum.b37.bwa-mem.bam \
       NA12892.platinum.b37.bwa-mem.bam \
    -t hs37d5.chromosomes.bed \
    --ped ceu-trio.ped \
    --denovos-only \
    -o octopus.ceu-trio.platinum.b37.bwa-mem.denovo.vcf.gz
```

#### Platypus

Call variants:

```shell
$ python Platypus.py callVariants \
    --refFile hs37d5.fa \
    --bamFiles NA12878.platinum.b37.bwa-mem.bam \
               NA12891.platinum.b37.bwa-mem.bam \
               NA12892.platinum.b37.bwa-mem.bam \
    --regions hs37d5.chromosomes.bed \
    --output platypus.ceu-trio.platinum.b37.bwa-mem.vcf
$ bgzip platypus.ceu-trio.platinum.b37.bwa-mem.vcf
$ tabix platypus.ceu-trio.platinum.b37.bwa-mem.vcf.gz 
```

Filter *de novos*

```shell
$ bcftools view \
    platypus.ceu-trio.platinum.b37.bwa-mem.vcf.gz \
    | python findDeNovoMutations.py 11 9 10 10 \
    | bgzip > platypus.ceu-trio.platinum.b37.bwa-mem.denovo.vcf.gz
$ tabix platypus.ceu-trio.platinum.b37.bwa-mem.denovo.vcf.gz
```

#### Strelka2

Call variants:

```shell
$ configureStrelkaGermlineWorkflow.py \
    --referenceFasta hs37d5.fa \
    --bam NA12878.platinum.b37.bwa-mem.bam \
    --bam NA12891.platinum.b37.bwa-mem.bam \
    --bam NA12892.platinum.b37.bwa-mem.bam \
    --callRegions hs37d5.chromosomes.bed.gz \
    --runDir strelka_tmp
$ strelka_tmp/runWorkflow.py -m local -j 10
$ mv strelka_tmp/results/variants/variants.vcf.gz strelka2.ceu-trio.platinum.b37.bwa-mem.vcf.gz
$ mv strelka_tmp/results/variants/variants.vcf.gz.tbi strelka2.ceu-trio.platinum.b37.bwa-mem.vcf.gz.tbi
```

Remove 'uncalled' and bad ploidy variants (we needed to do this otherwise the following step errors):

```shell
$ bcftools view \
    -a -U -c 1 \
    -e 'GT="0"|GT="1"' \
    -Oz -o strelka2.ceu-trio.platinum.b37.bwa-mem.corrected.vcf.gz \
    strelka2.ceu-trio.platinum.b37.bwa-mem.vcf.gz.tbi
$ tabix strelka2.ceu-trio.platinum.b37.bwa-mem.corrected.vcf.gz
```

Annotate *de novo* mutations:

```shell
$ java -jar GenomeAnalysisTK.jar -T VariantAnnotator \
    -R hs37d5.fa \
    -V strelka2.ceu-trio.platinum.b37.bwa-mem.corrected.vcf.gz \
    -A PossibleDeNovo \
    -ped ceu-trio.ped \
    -o strelka2.ceu-trio.platinum.b37.bwa-mem.corrected.denovo_annotated.vcf.gz
```

Filter *de novos*:

```shell
$ bcftools view \
    -i 'INFO/hiConfDeNovo="NA12878"' \
    -Oz -o strelka2.ceu-trio.platinum.b37.bwa-mem.denovo.vcf.gz \
    strelka2.ceu-trio.platinum.b37.bwa-mem.corrected.denovo_annotated.vcf.gz
```

Filter passing calls with no sample filters:

```shell
$ bcftools view \
    strelka2.ceu-trio.platinum.b37.bwa-mem.denovo.vcf.gz \
    | awk -F $"\t" \
        'BEGIN {OFS = FS} {if (NF>1 && $7 != "PASS" && $9~/:FT:/) $7="FAIL"; print} \
        | bgzip > strelka.aw_sc_trio.wgs500.bwa.b37.denovo.filter.vcf.gz
$ tabix strelka.aw_sc_trio.wgs500.bwa.b37.denovo.filter.vcf.gz
```

# Synthetic tumour

## Commands

#### GATK Mutect2

Mark duplicates

```shell
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12878.TUMOUR.60x.skin.b37.bwa-mem.bam \
     -O NA12878.TUMOUR.60x.skin.b37.bwa-mem.dedup.bam \
     -M NA12878.TUMOUR.60x.skin.b37.bwa-mem.dedup.metrics
$ gatk --java-options -Xmx12G MarkDuplicates \
     -I NA12878.NORMAL.30x.b37.bwa-mem.bam \
     -O NA12878.NORMAL.30x.b37.bwa-mem.dedup.bam \
     -M NA12878.NORMAL.30x.b37.bwa-mem.dedup.metrics
```

Recalibrate base qualities:

```shell
$ gatk --java-options -Xmx4G BaseRecalibrator \
    -R hs37d5.fa \
    -I NA12878.TUMOUR.60x.skin.b37.bwa-mem.dedup.bam \
    --known-sites gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    --known-sites gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --known-sites gatk_bundle/b37/1000G_phase1.indels.b37.vcf.gz \
    -O NA12878.TUMOUR.60x.skin.b37.bwa-mem.dedup.recal.table
$ gatk --java-options -Xmx4G ApplyBQSR
    -R hs37d5.fa \
    -I NA12878.TUMOUR.60x.skin.b37.bwa-mem.dedup.bam \
    --bqsr-recal-file NA12878.TUMOUR.60x.skin.b37.bwa-mem.dedup.table \
    -O NA12878.TUMOUR.60x.skin.b37.bwa-mem.dedup.recal.bam
$ gatk --java-options -Xmx4G BaseRecalibrator \
    -R hs37d5.fa \
    -I NA12878.NORMAL.30x.b37.bwa-mem.dedup.bam \
    --known-sites gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    --known-sites gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
    --known-sites gatk_bundle/b37/1000G_phase1.indels.b37.vcf.gz \
    -O NA12878.NORMAL.30x.b37.bwa-mem.dedup.recal.table
$ gatk --java-options -Xmx4G ApplyBQSR
    -R hs37d5.fa \
    -I NA12878.NORMAL.30x.b37.bwa-mem.dedup.bam \
    --bqsr-recal-file NA12878.NORMAL.30x.b37.bwa-mem.dedup.recal.table \
    -O NA12878.NORMAL.30x.b37.bwa-mem.dedup.recal.bam
```

Call variants:

```shell
$ gatk --java-options -Xmx12G Mutect2 \
    -R hs37d5.fa \
    -I NA12878.TUMOUR.60x.skin.b37.bwa-mem.dedup.recal.bam \
    -I NA12878.NORMAL.30x.b37.bwa-mem.dedup.recal.bam \
    -tumor NA12878.TUMOUR \
    -normal NA12878.NORMAL \
    --germline-resource gatk_bundle/b37/af-only-gnomad.raw.sites.b37.vcf.gz \
    --af-of-alleles-not-in-resource 0.0000025 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    -L hs37d5.chromosomes.bed \
    -O mutect2.raw.NA12878.syntumour.skin.bwa.dedup.recal.b37.vcf.gz
```

Get pileup summaries:

```shell
$ gatk GetPileupSummaries \
    -I NA12878.TUMOUR.60x.skin.b37.bwa-mem.dedup.recal.bam \
    -V gatk_bundle/b37/small_exac_common_3_b37.vcf.gz \
    -O tumor_getpileupsummaries.table
```

Calculate contamination:

```shell
$ gatk CalculateContamination \
    -I tumor_getpileupsummaries.table \
    -O tumor_calculatecontamination.table
```

Filter calls:

```shell
$ gatk FilterMutectCalls \
    -V mutect2.raw.NA12878.syntumour.skin.bwa-mem.dedup.recal.b37.vcf.gz \
    --contamination-table tumor_calculatecontamination.table \
    -O mutect2.NA12878.syntumour.skin.bwa-mem.dedup.recal.b37.vcf.gz
```

#### Lancet

Call variants:

```shell
for chrom in `seq 1 22` X; do
    lancet \
        --ref hs37d5.fa \
        --normal NA12878.NORMAL.30x.bwa.b37.bam \
        --tumor NA12878.TUMOUR.60x.skin.bwa-mem.b37.bam \
        --reg $chrom \
        --num-threads 10 \
        | bgzip > lancet.NA12878.syntumour.skin.bwa-mem.b37.chr${chrom}.vcf.gz
        tabix lancet.NA12878.syntumour.skin.bwa-mem.b37.chr${chrom}.vcf.gz
done
```

Merge results:

```shell
$ 
```

#### LoFreq

Call variants

```shell
$ lofreq somatic \
    -f hs37d5.fa \
    -n NA12878.NORMAL.30x.b37.bwa-mem.bam \
    -t NA12878.TUMOUR.60x.skin.b37.bwa-mem.bam \
    --call-indels \
    -l hs37d5.chromosomes.bed \
    -d gatk_bundle/b37/dbsnp_138.b37.vcf.gz \
    --threads 10 \
    -o lofreq.NA12878.syntumour.skin.b37.bwa-mem.
```

Merge SNVs and indels:

```shell
$ bcftools concat \
-a -Oz \
-o lofreq.NA12878.syntumour.skin.b37.bwa-mem.somatic_final_minus-dbsnp.snvs.indels.vcf.gz \
lofreq.NA12878.syntumour.skin.b37.bwa-mem.somatic_final_minus-dbsnp.snvs.vcf.gz \lofreq.NA12878.syntumour.skin.b37.bwa-mem.somatic_final_minus-dbsnp.indels.vcf.gz
$ tabix lofreq.NA12878.syntumour.skin.b37.bwa-mem.somatic_final_minus-dbsnp.snvs.indels.vcf.gz
```

#### Octopus

Call somatics:

```shell
$ octopus \
    -R hs37d5.fa
    -I NA12878.NORMAL.30x.b37.bwa-mem.bam \
       NA12878.TUMOUR.60x.skin.b37.bwa-mem.bam \
    -t hs37d5.chromosomes.bed \
    -N NA12878.NORMAL \
    --forest octopus-somatic.forest \
    --somatics-only \
    --threads 10 \
    -o octopus.NA12878.syntumour.skin.b37.bwa-mem.somatic.vcf.gz
```

#### Platypus

Call Variants:

```shell
$ python Platypus.py callVariants \
    --refFile hs37d5.fa \
    --bamFiles NA12878.NORMAL.30x.bwa.b37.bam \
               NA12878.TUMOUR.60x.skin.b37.bwa-mem.bam \
    --regions hs37d5.chromosomes.bed \
    --output platypus.NA12878.TUMOUR.60x.skin.b37.bwa-mem.vcf
```

Filter somatics:

```shell
$ findSomaticMutationsInTumour.py \
    --inputVCF platypus.NA12878.TUMOUR.60x.skin.b37.bwa-mem.vcf \
    --outputVCF platypus.NA12878.TUMOUR.60x.skin.b37.bwa-mem.somatic.vcf \
    --tumourSample NA12878.TUMOUR \
    --normalSample NA12878.NORMAL \
    --minPosterior 1
$ bgzip platypus.NA12878.TUMOUR.60x.skin.b37.bwa-mem.somatic.vcf
$ tabix platypus.NA12878.TUMOUR.60x.skin.b37.bwa-mem.somatic.vcf.gz
```

#### Strelka2

Call variants:

```shell
$ configureStrelkaSomaticWorkflow.py \
    --referenceFasta hs37d5.fa \
    --normalBam NA12878.NORMAL.30x.bwa.b37.bam \
    --tumorBam NA12878.TUMOUR.60x.skin.b37.bwa-mem.bam \
    --callRegions hs37d5.chromosomes.bed.gz \
    --runDir strelka_tmp
$ strelka_tmp/runWorkflow.py -m local -j 10
```

Merge SNVs and indels:

```shell
$ bcftools concat -a \
    -Oz -o strelka2.NA12878.syntumour.skin.b37.bwa-mem.vcf.gz \
    strelka_tmp/results/variants/somatic.snvs.vcf.gz \
    strelka_tmp/results/variants/somatic.indels.vcf.gz
$ tabix strelka2.NA12878.syntumour.skin.b37.bwa-mem.vcf.gz
```

#### VarDict

Call variants:

```shell
$ VarDict \
    -G hs37d5.fa \
    -f 0.01 \
    -N NA12878.TUMOUR \
    -b "NA12878.TUMOUR.60x.skin.bwa-mem.b37.bam|NA12878.NORMAL.30x.bwa.b37.bam" \
    -z 0 -c 1 -S 2 -E 3 \
    hs37d5.chromosomes.bed \
    | testsomatic.R \
    | var2vcf_paired.pl \
      -N "NA12878.TUMOUR|NA12878.NORMAL" \
      -f 0.01 \
      | bgzip > vardict.NA12878.syntumour.skin.b37.bwa-mem.vcf.gz
$ tabix vardict.NA12878.syntumour.skin.b37.bwa-mem.vcf.gz
```

Filter somatics:

```shell
$ bcftools view \
    -i 'INFO/STATUS="StrongSomatic"' \
    -f PASS \
    -Oz -o vardict.NA12878.syntumour.skin.b37.bwa-mem.somatic.PASS.vcf.gz \
    vardict.NA12878.syntumour.skin.b37.bwa-mem.vcf.gz
```

Apply bcbio recommended filtera:

```shell
$ bcftools filter \
    -s "BCBIO" \
    -e '((INFO/AF * INFO/DP < 6) && ((MQ < 55.0 && NM > 1.0) || (MQ < 60.0 && NM > 2.0) || (FORMAT/DP < 10) || (QUAL < 45)))' \
    -Oz -o VarDict.vcf.gz \
    vardict.NA12878.syntumour.skin.b37.bwa-mem.vcf.gz
```

## Tumour-only

#### Pisces

Make genome size file:

```shell
$ dotnet CreateGenomeSizeFile.dll \
    -g hs37d5 \
    -s "Homo sapien (b37)" \
    -o hs37d5
```

Call somatic and germline variants:

```shell
$ dotnet Pisces.dll \
    -bam NA12878.TUMOUR.65x.breast.bwa-mem.b37.bam \
    -g hs37d5 \
    -i hs37d5.chromosomes.bed \
    -CallMNVs false \
    -gVCF false \
    -RMxNFilter 5,9,0.35 \
    -OutFolder pisces.NA12878.syntumour.skin.unpaired.b37.bwa-mem.vcf.gz
```
Call germline variants:

```shell
$ dotnet Pisces.dll \
    -bam NA12878.TUMOUR.65x.breast.bwa-mem.b37.bam \
    -g hs37d5 \
    -i hs37d5.chromosomes.bed \
    -CallMNVs false \
    -gVCF false \
    -RMxNFilter 5,9,0.35 \
    -ploidy diploid \
    -OutFolder pisces.NA12878.syntumour.skin.unpaired.germline.b37.bwa-mem.vcf.gz
```

Subtract germline to get somatic calls:

```shell
$ dotnet VennVcf.dll \
    
```

## Miscellaneous commands

One-liner to identify all PASSing microinversions:

```shell
$ bcftools view -v mnps -f PASS in.vcf | bioawk -tHc vcf '{ if ($ref==revcomp($alt)) print; }'
```

One-liner to identify PASSing microinversions longer than 3 bases:

```shell
$ bcftools view -v mnps -f PASS in.vcf | bioawk -tHc vcf '{ if (length($ref) > 3 && $ref==revcomp($alt)) print; }'
```

Remove `chr` from a VCF:

```shell
$ bcftools view ${chr_vcf} \
    | awk '{if(match($0,/(##contig=<ID=)(.*)/,m)) gsub("chr", ""); else gsub(/^chr/,""); print}' \
    | bgzip > ${no_chr_vcf}
```

#### Identifying likely cell line artefacts

```shell
$  bcftools isec Lancet/fp.vcf.gz LoFreq/fp.vcf.gz Mutect2/fp.vcf.gz Platypus/fp.vcf.gz Strelka2/fp.vcf.gz VarDict/fp.vcf.gz octopus/fp.vcf.gz -n +5 -w 6 | bioawk -tc vcf '{print $chrom, $pos-1, $pos+1}' > artefacts.bed
``` 

#### Identifying BAMSurgeon missed spike-ins:

Find unanimous false negatives and their regions:

```shell
$ bcftools isec Lancet/fn.vcf.gz LoFreq/fn.vcf.gz Mutect2/fn.vcf.gz Platypus/fn.vcf.gz Strelka2/fn.vcf.gz VarDict/fn.vcf.gz octopus/fn.vcf.gz -n7 -w1 | bgzip > unanimous_fn.vcf.gz
$ bioawk -tc vcf '{print $chrom, $pos-50, $pos+50}' unanimous_fn.vcf.gz > unanimous_fn_padded_regions.bed

```

Select those regions that have two or more truth variants:

```python
import pysam as ps
import csv

truth_vcf = ps.VariantFile("truth.vcf.gz")

with open("unanimous_fn_padded_regions.bed", encoding='latin-1') as in_bed:
    with open("fns_probably_not_spiked_in.bed", 'w') as out_bed: 
        bedreader = csv.reader(in_bed, delimiter='\t')
        bedwriter = csv.writer(out_bed, delimiter='\t')
        for row in bedreader:
            if (sum(1 for _ in truth_vcf.fetch(row[0], int(row[1]), int(row[2]))) > 1):
                bedwriter.writerow(row)
```

Copy the unanimous false negatives that fall in these regions and make a padded BED around them:

```shell
$ bcftools view -T fns_probably_not_spiked_in.bed unanimous_fn.vcf.gz > missing_truth_variants.vcf
$ bioawk -tc vcf '{print $chrom, $pos-1, $pos+1}' missing_truth_variants.vcf > missing_truth_regions.bed
```

Merge BEDs:

```shell
$ cat a.bed b.bed | bedtools sort -i stdin -g hs37d5.chromosomes.genome | bedtools merge -i stdin > c.bed
```
Complement BED:

```shell
$ bedtools complement -i remove.bed -g hs37d5.chromosomes.genome > final.bed
```
