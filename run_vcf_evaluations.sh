#!/bin/bash
module load samtools bcftools singularity

WORKDIR="$HOME/hprc_workdir"
mkdir -p $WORKDIR && cd $WORKDIR

## COLLECT AND EVALUATE DEEPVARIANT CALLS
cd ${WORKDIR}/grch38_pgrc_decoys_july7_HG003_bam_by_chr
bcftools concat -O v *deepvariant.vcf.gz > giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.vcf
bcftools sort giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.vcf -O v > giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.vcf && rm giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.vcf
bgzip giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.vcf
tabix -p vcf giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.vcf.gz
mv giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.vcf.gz giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.vcf.gz.tbi ../

cd ${WORKDIR}
bgzip -d -f giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.vcf.gz
awk '{gsub(/GRCh38.chr/,"chr"); print}' giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.vcf > giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.normal_contig_names.vcf
bgzip -f giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.normal_contig_names.vcf
tabix -p vcf giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.normal_contig_names.vcf.gz

## EVALUATE ALL CONFIDENT BENCHMARK REGIONS
cd ${WORKDIR}
singularity run -H ${PWD}:${HOME} --pwd ${HOME} docker://jmcdani20/hap.py:v0.3.12 \
time /opt/hap.py/bin/hap.py \
HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.normal_contig_names.vcf.gz \
-f HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
--reference GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
--threads 32 \
--engine=vcfeval \
-o hg003_grch38_hprc_deepvariant_happy_vcfeval


## EVALUATE DIFFICULT REGIONS
cd ${WORKDIR}

ANALYSIS_DIR="${WORKDIR}/OverviewRegions"
mkdir -p ${ANALYSIS_DIR}
rm overview_region.happy.GIRAFFE_PGRC_HG003.swarm
for REGION in "GRCh38_MHC" "GRCh38_alllowmapandsegdupregions" "GRCh38_alldifficultregions" "GRCh38_notinalldifficultregions" "GRCh38_notinalllowmapandsegdupregions" ; do
    echo -e "module load singularity \
    && cd ${WORKDIR} \
    && singularity run -H \${PWD}:\${HOME} --pwd \${HOME} docker://jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
    HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.normal_contig_names.vcf.gz \
    -f HG003_GRCh38_1_22_v4.2.1.${REGION}.chr.bed \
    --reference GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
    --threads 32 \
    --engine=vcfeval \
    -o OverviewRegions/happy_GIRAFFE_PGRC_HG003_${REGION}" >> overview_region.happy.GIRAFFE_PGRC_HG003.swarm
done

cd ${WORKDIR}
echo "running overview regions"
swarm -f overview_region.happy.GIRAFFE_PGRC_HG003.swarm -t 32 -g 60 --partition=norm --time=4:00:00 --module singularity

ANALYSIS_DIR="${WORKDIR}/low_mappability"
mkdir -p ${ANALYSIS_DIR}
rm low_map_region.happy.GIRAFFE_PGRC_HG003.swarm
for REGION in "GRCh38_lowmappabilityall" "GRCh38_nonunique_l100_m2_e1" "GRCh38_nonunique_l250_m0_e0" "GRCh38_notinlowmappabilityall" ; do
    echo -e "module load singularity \
    && cd ${WORKDIR} \
    && singularity run -H \${PWD}:\${HOME} --pwd \${HOME} docker://jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
    HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.normal_contig_names.vcf.gz \
    -f HG003_GRCh38_1_22_v4.2.1.${REGION}.chr.bed \
    --reference GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
    --threads 32 \
    --engine=vcfeval \
    -o low_mappability/happy_GIRAFFE_PGRC_HG003_${REGION}" >> low_map_region.happy.GIRAFFE_PGRC_HG003.swarm
done

cd ${WORKDIR}
echo "running low mappable regions"
swarm -f low_map_region.happy.GIRAFFE_PGRC_HG003.swarm -t 32 -g 60 --partition=norm --time=4:00:00 --module singularity


ANALYSIS_DIR="${WORKDIR}/seg_duplicate_regions"
mkdir -p ${ANALYSIS_DIR}
rm seg_dup_region.happy.GIRAFFE_PGRC_HG003.swarm
for REGION in "GRCh38_chainSelf" "GRCh38_chainSelf_gt10kb" "GRCh38_gt5segdups_gt10kb_gt99percidentity" "GRCh38_notinchainSelf" "GRCh38_notinchainSelf_gt10kb" "GRCh38_notinsegdups" "GRCh38_notinsegdups_gt10kb" "GRCh38_segdups" "GRCh38_segdups_gt10kb" ; do
    echo -e "module load singularity \
    && cd ${WORKDIR} \
    && singularity run -H \${PWD}:\${HOME} --pwd \${HOME} docker://jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
    HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.normal_contig_names.vcf.gz \
    -f HG003_GRCh38_1_22_v4.2.1.${REGION}.chr.bed \
    --reference GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
    --threads 32 \
    --engine=vcfeval \
    -o seg_duplicate_regions/happy_GIRAFFE_PGRC_HG003_${REGION}" >> seg_dup_region.happy.GIRAFFE_PGRC_HG003.swarm
done

cd ${WORKDIR}
echo "running seg dup regions"
swarm -f seg_dup_region.happy.GIRAFFE_PGRC_HG003.swarm -t 32 -g 60 --partition=norm --time=4:00:00 --module singularity

ANALYSIS_DIR="${WORKDIR}/low_complex_regions"
mkdir -p ${ANALYSIS_DIR}
rm low_complexity_region.happy.GIRAFFE_PGRC_HG003.swarm
for REGION in "GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5" "GRCh38_AllTandemRepeats_gt10000bp_slop5" "GRCh38_notinAllTandemRepeatsandHomopolymers_slop5" "GRCh38_SimpleRepeat_homopolymer_4to6_slop5" "GRCh38_SimpleRepeat_quadTR_20to50_slop5" "GRCh38_SimpleRepeat_triTR_51to200_slop5" "GRCh38_AllTandemRepeats_201to10000bp_slop5" "GRCh38_AllTandemRepeats_gt100bp_slop5" "GRCh38_SimpleRepeat_diTR_11to50_slop5" "GRCh38_SimpleRepeat_homopolymer_7to11_slop5" "GRCh38_SimpleRepeat_quadTR_51to200_slop5" "GRCh38_SimpleRepeat_triTR_gt200_slop5" "GRCh38_AllTandemRepeats_51to200bp_slop5" "GRCh38_AllTandemRepeats_lt51bp_slop5" "GRCh38_SimpleRepeat_diTR_51to200_slop5" "GRCh38_SimpleRepeat_homopolymer_gt11_slop5" "GRCh38_SimpleRepeat_quadTR_gt200_slop5" "GRCh38_AllTandemRepeatsandHomopolymers_slop5" "GRCh38_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5" "GRCh38_SimpleRepeat_diTR_gt200_slop5" "GRCh38_SimpleRepeat_imperfecthomopolgt10_slop5" "GRCh38_SimpleRepeat_triTR_15to50_slop5" ; do
    echo -e "module load singularity \
    && cd ${WORKDIR} \
    && singularity run -H \${PWD}:\${HOME} --pwd \${HOME} docker://jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
    HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    giraffe_grch38_pgrc_decoys_july7_HG003.merged.deepvariant.sorted.normal_contig_names.vcf.gz \
    -f HG003_GRCh38_1_22_v4.2.1.${REGION}.chr.bed \
    --reference GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
    --threads 32 \
    --engine=vcfeval \
    -o low_complex_regions/happy_GIRAFFE_PGRC_HG003_${REGION}" >> low_complexity_region.happy.GIRAFFE_PGRC_HG003.swarm
done

cd ${WORKDIR}
echo "running low complex regions"
swarm -f low_complexity_region.happy.GIRAFFE_PGRC_HG003.swarm -t 32 -g 60 --partition=norm --time=4:00:00 --module singularity




