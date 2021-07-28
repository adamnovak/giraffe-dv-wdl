#!/bin/bash
module load singularity samtools picard bcftools

WORKDIR="$HOME/hprc_workdir"
mkdir -p $WORKDIR && cd $WORKDIR

## DOWNLOAD INPUT DATA
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/debug_pgrc_graph/july7_version_analysis/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7.fna .
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/debug_pgrc_graph/july7_version_analysis/GRCh38-f1g-90-mc-jul7.xg .
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/debug_pgrc_graph/july7_version_analysis/GRCh38-f1g-90-mc-jul7.gbwt .
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/debug_pgrc_graph/july7_version_analysis/GRCh38-freeze1.cactus.f1g-90-mc-jul7.k29.w11.paths.min .
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/debug_pgrc_graph/july7_version_analysis/GRCh38-freeze1.cactus.f1g-90-mc-jul7.dist .
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/debug_pgrc_graph/july7_version_analysis/HG003.novaseq.pcr-free.35x.R1.fastq.gz .
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/debug_pgrc_graph/july7_version_analysis/HG003.novaseq.pcr-free.35x.R2.fastq.gz .
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/benchmark_data/hg002_cohort/hg003/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz .
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/benchmark_data/hg002_cohort/hg003/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi .
wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/benchmark_data/hg002_cohort/hg003/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed .

for REGION in "GRCh38_MHC" "GRCh38_alllowmapandsegdupregions" "GRCh38_alldifficultregions" "GRCh38_notinalldifficultregions" "GRCh38_notinalllowmapandsegdupregions" "GRCh38_lowmappabilityall" "GRCh38_nonunique_l100_m2_e1" "GRCh38_nonunique_l250_m0_e0" "GRCh38_notinlowmappabilityall" "GRCh38_chainSelf" "GRCh38_chainSelf_gt10kb" "GRCh38_gt5segdups_gt10kb_gt99percidentity" "GRCh38_notinchainSelf" "GRCh38_notinchainSelf_gt10kb" "GRCh38_notinsegdups" "GRCh38_notinsegdups_gt10kb" "GRCh38_segdups" "GRCh38_segdups_gt10kb" "GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5" "GRCh38_AllTandemRepeats_gt10000bp_slop5" "GRCh38_notinAllTandemRepeatsandHomopolymers_slop5" "GRCh38_SimpleRepeat_homopolymer_4to6_slop5" "GRCh38_SimpleRepeat_quadTR_20to50_slop5" "GRCh38_SimpleRepeat_triTR_51to200_slop5" "GRCh38_AllTandemRepeats_201to10000bp_slop5" "GRCh38_AllTandemRepeats_gt100bp_slop5" "GRCh38_SimpleRepeat_diTR_11to50_slop5" "GRCh38_SimpleRepeat_homopolymer_7to11_slop5" "GRCh38_SimpleRepeat_quadTR_51to200_slop5" "GRCh38_SimpleRepeat_triTR_gt200_slop5" "GRCh38_AllTandemRepeats_51to200bp_slop5" "GRCh38_AllTandemRepeats_lt51bp_slop5" "GRCh38_SimpleRepeat_diTR_51to200_slop5" "GRCh38_SimpleRepeat_homopolymer_gt11_slop5" "GRCh38_SimpleRepeat_quadTR_gt200_slop5" "GRCh38_AllTandemRepeatsandHomopolymers_slop5" "GRCh38_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5" "GRCh38_SimpleRepeat_diTR_gt200_slop5" "GRCh38_SimpleRepeat_imperfecthomopolgt10_slop5" "GRCh38_SimpleRepeat_triTR_15to50_slop5" ; do
    wget -c https://storage.googleapis.com/cmarkell-vg-wdl-dev/benchmark_data/hg002_cohort/hg003/HG003_GRCh38_1_22_v4.2.1.${REGION}.chr.bed
done

samtools faidx GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7.fna
time java -Xmx100g -XX:ParallelGCThreads=32 -jar $PICARDJARPATH/picard.jar \
CreateSequenceDictionary \
R=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7.fna \
O=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7.dict

## RUN GIRAFFE ALIGNMENT
singularity exec -H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/jmonlong/vg:beea35e \
vg giraffe \
-x GRCh38-f1g-90-mc-jul7.xg \
-H GRCh38-f1g-90-mc-jul7.gbwt \
-m GRCh38-freeze1.cactus.f1g-90-mc-jul7.k29.w11.paths.min \
-d GRCh38-freeze1.cactus.f1g-90-mc-jul7.dist \
-g GRCh38-f1g-90-mc-jul7.gg \
-f HG003.novaseq.pcr-free.35x.R1.fastq.gz \
-f HG003.novaseq.pcr-free.35x.R2.fastq.gz \
--output-format GAM \
-p \
-t 56 > giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.gam

singularity exec -H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/jmonlong/vg:beea35e \
vg surject \
giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.gam \
-x GRCh38-f1g-90-mc-jul7.xg \
-b -i \
-t 70 > giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.bam

samtools sort --threads 32 -O BAM giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.bam > giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.positionsorted.bam && rm giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.bam
samtools index -@ 32 giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.positionsorted.bam

time java -Xmx100g -XX:ParallelGCThreads=32 -jar $PICARDJARPATH/picard.jar \
ReorderSam \
VALIDATION_STRINGENCY=SILENT \
INPUT=giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.positionsorted.bam \
OUTPUT=giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.positionsorted.reordered.bam \
SEQUENCE_DICTIONARY=/data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7.dict && rm giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.positionsorted.bam

## indel realign BAMs
cd ${WORKDIR}
samtools index -@ 32 giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.positionsorted.reordered.bam
mkdir -p ${WORKDIR}/grch38_pgrc_decoys_july7_HG003_bam_by_chr
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    samtools view -@ 32 -b giraffe_all_wgs.grch38_pgrc_decoys_july7.HG003.positionsorted.reordered.bam GRCh38.chr${CHR} > ${WORKDIR}/grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.bam
    samtools index -@ 32 ${WORKDIR}/grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.bam
done

cd ${WORKDIR}
rm indel_realignment_giraffe_grch38_pgrc_decoys_july7.HG003.swarm
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    echo -e "module load singularity GATK/3.8-1 samtools \
    && cd ${WORKDIR} \
    && samtools addreplacerg -@ 32 -O BAM -r ID:1 -r LB:lib1 -r SM:HG003 -r PL:illumina -r PU:unit1 grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.bam > grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.gatk_ready.bam \
    && samtools index -@ 32 grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.gatk_ready.bam \
    && java -jar \$GATK_JAR -T RealignerTargetCreator -R ${WORKDIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7.fna -L GRCh38.chr${CHR} -I grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.gatk_ready.bam -o grch38_pgrc_decoys_july7_HG003_bam_by_chr/${CHR}.intervals \
    && awk -F '[:-]' 'BEGIN { OFS = \"\t\" } { if( \$3 == \"\") { print \$1, \$2-1, \$2 } else { print \$1, \$2-1, \$3}}' grch38_pgrc_decoys_july7_HG003_bam_by_chr/${CHR}.intervals > grch38_pgrc_decoys_july7_HG003_bam_by_chr/${CHR}.intervals.bed \
    && singularity run -H \${PWD}:\${HOME} --pwd \${HOME} docker://dceoy/abra2:latest --targets grch38_pgrc_decoys_july7_HG003_bam_by_chr/${CHR}.intervals.bed --in grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.gatk_ready.bam --out grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.indel_realigned.bam --ref GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7.fna --threads 16 && samtools index -@ 32 grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.indel_realigned.bam" >> indel_realignment_giraffe_grch38_pgrc_decoys_july7.HG003.swarm
done

indel_realign_swarm_jobid=$(swarm -f indel_realignment_giraffe_grch38_pgrc_decoys_july7.HG003.swarm -t 16 -g 100 --partition=norm --time=12:00:00 --module singularity,GATK/3.8-1,samtools)

## SETUP AND RUN DEEPVARIANT
cd ${WORKDIR}
rm deepvariant_calling.minmapq0.swarm
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    echo -e "module load singularity; cd ${WORKDIR}; singularity run --nv -H \${PWD}:\${HOME} --pwd \${HOME} docker://google/deepvariant:1.1.0-gpu /opt/deepvariant/bin/run_deepvariant --make_examples_extra_args 'min_mapping_quality=0' --model_type=WGS --regions GRCh38.chr${CHR} --ref=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7.fna --reads=grch38_pgrc_decoys_july7_HG003_bam_by_chr/raw.chr${CHR}.indel_realigned.bam --output_vcf=grch38_pgrc_decoys_july7_HG003_bam_by_chr/chr${CHR}.deepvariant.vcf.gz --output_gvcf=grch38_pgrc_decoys_july7_HG003_bam_by_chr/chr${CHR}.g.vcf.gz --intermediate_results_dir=grch38_pgrc_decoys_july7_HG003_bam_by_chr/tmp_deepvariant_chr${CHR} --num_shards=16" >> deepvariant_calling.minmapq0.swarm
done

cd ${WORKDIR}
swarm -f deepvariant_calling.minmapq0.swarm --dependency=afterok:$indel_realign_swarm_jobid -t 16 -g 50 --partition=gpu --gres=gpu:k80:1 --time=12:00:00 --module singularity)




