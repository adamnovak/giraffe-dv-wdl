version 1.0


workflow generateTrainingData {
    input {
        Array[File] INPUT_READ_FILE_1                   # Input samples' 1st read pair fastq.gz files
        Array[File] INPUT_READ_FILE_2                   # Input samples' 2nd read pair fastq.gz files
        Array[String] SAMPLE_NAME                       # The sample names
        # VG Container used in the pipeline (e.g. quay.io/vgteam/vg:v1.16.0)
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.36.0"
        Int READS_PER_CHUNK = 20000000                  # Number of reads contained in each mapping chunk (20000000 for wgs)
        String? GIRAFFE_OPTIONS                         # (OPTIONAL) extra command line options for Giraffe mapper
        Array[String]+? CONTIGS                         # (OPTIONAL) Desired reference genome contigs, which are all paths in the XG index.
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index, to use instead of CONTIGS. If neither is given, paths are extracted from the XG and subset to chromosome-looking paths.
        File XG_FILE                                    # Path to .xg index file
        File GBWT_FILE                                  # Path to .gbwt index file
        File GGBWT_FILE                                 # Path to .gg index file
        File DIST_FILE                                  # Path to .dist index file
        File MIN_FILE                                   # Path to .min index file
        Int REALIGNMENT_EXPANSION_BASES = 160           # Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions.
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_DISK = 10
        Int MAP_CORES = 16
        Int MAP_DISK = 10
        Int MAP_MEM = 50
        Int REALIGN_DISK = 40
    }

    # Prepare the reference
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from xg file if PATH_LIST_FILE input not provided
            call extractPathNames {
                input:
                    in_xg_file=XG_FILE,
                    in_vg_container=VG_CONTAINER,
                    in_extract_disk=MAP_DISK,
                    in_extract_mem=MAP_MEM
            }
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call subsetPathNames {
                input:
                    in_path_list_file=extractPathNames.output_path_list_file
            }
        }
    } 
    if (defined(CONTIGS)) {
        # Put the paths in a file to use later. We know the value is defined,
        # but WDL is a bit low on unboxing calls for optionals so we use
        # select_first.
        File written_path_names_file = write_lines(select_first([CONTIGS]))
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, subsetPathNames.output_path_list_file, written_path_names_file])
    
    # To make sure that we have a FASTA reference with a contig set that
    # exactly matches the graph, we generate it ourselves, from the graph.
    call extractReference {
        input:
            in_xg_file=XG_FILE,
            in_vg_container=VG_CONTAINER,
            in_extract_disk=MAP_DISK,
            in_extract_mem=MAP_MEM
    }
    File reference_file = extractReference.reference_file
    
    call indexReference {
        input:
            in_reference_file=reference_file,
            in_index_disk=MAP_DISK,
            in_index_mem=MAP_MEM
    }
    File reference_index_file = indexReference.reference_index_file
    File reference_dict_file = indexReference.reference_dict_file


    # Create samples in a useful per-sample internal format.
    # TODO: use a struct?
    Array[Pair[String, Pair[File, File]]] samples = zip(SAMPLE_NAME, zip(INPUT_READ_FILE_1, INPUT_READ_FILE_2))
    
    scatter(sample in samples) {
        # TODO: make this a per-sample workflow and share it with the calling workflow?
        # Would need to get set up somewhere that supports hosting workflows with imports (Dockstore?)

        # Unpack the sample
        String sample_name = sample.left
        File sample_read1 = sample.right.left
        File sample_read2 = sample.right.right

        # Split input reads into chunks for parallelized mapping
        call splitReads as firstReadPair {
            input:
                in_read_file=sample_read1,
                in_pair_id="1",
                in_vg_container=VG_CONTAINER,
                in_reads_per_chunk=READS_PER_CHUNK,
                in_split_read_cores=SPLIT_READ_CORES,
                in_split_read_disk=SPLIT_READ_DISK
        }
        call splitReads as secondReadPair {
            input:
                in_read_file=sample_read2,
                in_pair_id="2",
                in_vg_container=VG_CONTAINER,
                in_reads_per_chunk=READS_PER_CHUNK,
                in_split_read_cores=SPLIT_READ_CORES,
                in_split_read_disk=SPLIT_READ_DISK
        }

       
        ################################################################
        # Distribute vg mapping opperation over each chunked read pair #
        ################################################################
        Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
        scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
            call runVGGIRAFFE {
                input:
                    in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                    in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                    in_vg_container=VG_CONTAINER,
                    in_giraffe_options=GIRAFFE_OPTIONS,
                    in_xg_file=XG_FILE,
                    in_gbwt_file=GBWT_FILE,
                    in_ggbwt_file=GGBWT_FILE,
                    in_dist_file=DIST_FILE,
                    in_min_file=MIN_FILE,
                    in_ref_dict=reference_dict_file,
                    in_sample_name=sample_name,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
            call sortBAMFile {
                input:
                    in_sample_name=sample_name,
                    in_bam_chunk_file=runVGGIRAFFE.chunk_bam_file,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM,
            }
        }
        Array[File] alignment_chunk_bam_files = select_all(sortBAMFile.sorted_chunk_bam)

        call mergeAlignmentBAMChunks {
            input:
                in_sample_name=sample_name,
                in_alignment_bam_chunk_files=alignment_chunk_bam_files,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
        
        # Split merged alignment by contigs list
        call splitBAMbyPath { 
            input:
                in_sample_name=sample_name,
                in_merged_bam_file=mergeAlignmentBAMChunks.merged_bam_file,
                in_merged_bam_file_index=mergeAlignmentBAMChunks.merged_bam_file_index,
                in_path_list_file=pipeline_path_list_file,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
        scatter (deepvariant_caller_input_files in zip(splitBAMbyPath.bam_contig_files, splitBAMbyPath.bam_contig_files_index)) {
            # Just left-shift each read individually
            call leftShiftBAMFile {
                input:
                in_sample_name=sample_name,
                in_bam_file=deepvariant_caller_input_files.left,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_realign_disk=REALIGN_DISK
            }
            # This tool can't make an index itself so we need to re-index the BAM
            call indexBAMFile {
            input:
                in_sample_name=sample_name,
                in_bam_file=leftShiftBAMFile.left_shifted_bam,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
            } 
            call runGATKRealignerTargetCreator {
                input:
                    in_sample_name=sample_name,
                    in_bam_file=leftShiftBAMFile.left_shifted_bam,
                    in_bam_index_file=indexBAMFile.bam_index,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_reference_dict_file=reference_dict_file,
                    in_realign_disk=REALIGN_DISK
            }
            if (REALIGNMENT_EXPANSION_BASES != 0) {
                # We want the realignment targets to be wider
                call widenRealignmentTargets {
                    input:
                        in_target_bed_file=runGATKRealignerTargetCreator.realigner_target_bed,
                        in_reference_index_file=reference_index_file,
                        in_expansion_bases=REALIGNMENT_EXPANSION_BASES,
                        in_realign_disk=REALIGN_DISK
                }
            }
            File target_bed_file = select_first([widenRealignmentTargets.output_target_bed_file, runGATKRealignerTargetCreator.realigner_target_bed])
            call runAbraRealigner {
                input:
                    in_sample_name=sample_name,
                    in_bam_file=leftShiftBAMFile.left_shifted_bam,
                    in_bam_index_file=indexBAMFile.bam_index,
                    in_target_bed_file=target_bed_file,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_realign_disk=REALIGN_DISK
            }
            call sortBAMFile as sortRealignedBAMFile {
                input:
                    in_sample_name=sample_name,
                    in_bam_chunk_file=runAbraRealigner.indel_realigned_bam,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM,
            }
        }
        Array[File] shifted_chunk_bam_files = select_all(leftShiftBAMFile.left_shifted_bam)
        Array[File] realignment_chunk_bam_files = select_all(sortRealignedBAMFile.sorted_chunk_bam)
        
        call mergeAlignmentBAMChunks as mergeShiftedBAMChunks {
            input:
                in_sample_name=sample_name,
                in_alignment_bam_chunk_files=shifted_chunk_bam_files,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
        call mergeAlignmentBAMChunks as mergeRealignmentBAMChunks {
            input:
                in_sample_name=sample_name,
                in_suffix=".realigned",
                in_alignment_bam_chunk_files=realignment_chunk_bam_files,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }

    }
    
    # Collect the files we are interested in
    Array[File] output_aligned_bam = select_all(mergeAlignmentBAMChunks.merged_bam_file)
    Array[File] output_aligned_bam_index = select_all(mergeAlignmentBAMChunks.merged_bam_file_index)
    Array[File] output_shifted_bam = select_all(mergeShiftedBAMChunks.merged_bam_file)
    Array[File] output_shifted_bam_index = select_all(mergeShiftedBAMChunks.merged_bam_file_index)
    Array[File] output_realigned_bam = select_all(mergeRealignmentBAMChunks.merged_bam_file)
    Array[File] output_realigned_bam_index = select_all(mergeRealignmentBAMChunks.merged_bam_file_index)
   
    output {
        Array[Pair[File, File]] output_aligned_indexed_bam = zip(output_aligned_bam, output_aligned_bam_index)
        Array[Pair[File, File]] output_shifted_indexed_bam = zip(output_shifted_bam, output_shifted_bam_index)
        Array[Pair[File, File]] output_realigned_indexed_bam = zip(output_realigned_bam, output_realigned_bam_index)
    }   
}

########################
### TASK DEFINITIONS ###
########################

task splitReads {
    input {
        File in_read_file
        String in_pair_id
        String in_vg_container
        Int in_reads_per_chunk
        Int in_split_read_cores
        Int in_split_read_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        CHUNK_LINES=$(( ~{in_reads_per_chunk} * 4 ))
        gzip -cd ~{in_read_file} | split -l $CHUNK_LINES --filter='pigz -p ~{in_split_read_cores} > ${FILE}.fq.gz' - "fq_chunk_~{in_pair_id}.part."
    >>>
    output {
        Array[File] output_read_chunks = glob("fq_chunk_~{in_pair_id}.part.*")
    }
    runtime {
        preemptible: 2
        time: 120
        cpu: in_split_read_cores
        memory: "2 GB"
        disks: "local-disk " + in_split_read_disk + " SSD"
        docker: "quay.io/glennhickey/pigz:2.3.1"
    }
}

task extractPathNames {
    input {
        File in_xg_file
        String in_vg_container
        Int in_extract_disk
        Int in_extract_mem
    }

    command {
        set -eux -o pipefail

        vg paths \
            --list \
            --xg ${in_xg_file} > path_list.txt
    }
    output {
        File output_path_list_file = "path_list.txt"
    }
    runtime {
        preemptible: 2
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: in_vg_container
    }
}

task subsetPathNames {
    input {
        File in_path_list_file
    }

    command <<<
        set -eux -o pipefail

        grep -v _decoy ~{in_path_list_file} | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM > path_list.txt
    >>>
    output {
        File output_path_list_file = "path_list.txt"
    }
    runtime {
        preemptible: 2
        memory: "1 GB"
        disks: "local-disk 10 SSD"
        docker: "ubuntu:20.04"
    }
}

task extractReference {
    input {
        File in_xg_file
        String in_vg_container
        Int in_extract_disk
        Int in_extract_mem
    }

    command {
        set -eux -o pipefail

        vg paths \
            --extract-fasta \
            --xg ${in_xg_file} > ref.fa
    }
    output {
        File reference_file = "ref.fa"
    }
    runtime {
        preemptible: 2
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: in_vg_container
    }
}

task indexReference {
    input {
        File in_reference_file
        Int in_index_mem
        Int in_index_disk
    }

    command <<<
        set -eux -o pipefail
        
        ln -s ~{in_reference_file} ref.fa
        
        samtools faidx ref.fa
        
        # Save a reference copy by making the dict now
        java -jar /usr/picard/picard.jar CreateSequenceDictionary \
          R=ref.fa \
          O=ref.dict
    >>>
    output {
        File reference_index_file = "ref.fa.fai"
        File reference_dict_file = "ref.dict"
    }
    runtime {
        preemptible: 2
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task runVGGIRAFFE {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File in_gbwt_file
        File in_ggbwt_file
        File in_dist_file
        File in_min_file
        File in_ref_dict
        String in_vg_container
        String? in_giraffe_options
        String in_sample_name
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        READ_CHUNK_ID=($(ls ~{in_left_read_pair_chunk_file} | awk -F'.' '{print $(NF-2)}'))
        vg giraffe \
          --progress \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --sample "~{in_sample_name}" \
          --output-format BAM \
          ~{in_giraffe_options} \
          --ref-paths ~{in_ref_dict} \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -x ~{in_xg_file} \
          -H ~{in_gbwt_file} \
          -g ~{in_ggbwt_file} \
          -d ~{in_dist_file} \
          -m ~{in_min_file} \
          -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.bam
    >>>
    output {
        File chunk_bam_file = glob("*bam")[0]
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
    }
}

task sortBAMFile {
    input {
        String in_sample_name
        File in_bam_chunk_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        samtools sort \
          --threads ~{in_map_cores} \
          ~{in_bam_chunk_file} \
          -O BAM > ~{in_sample_name}.positionsorted.bam
    >>>
    output {
        File sorted_chunk_bam = "~{in_sample_name}.positionsorted.bam"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task indexBAMFile {
    input {
        String in_sample_name
        File in_bam_file
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        # Never use Samtools 1.4.1 here! See https://github.com/samtools/samtools/issues/687
        samtools index \
          ~{in_bam_file} \
          index.bai
    >>>
    output {
        File bam_index = "index.bai"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: in_map_mem + " GB"
        cpu: 1
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task mergeAlignmentBAMChunks {
    input {
        String in_sample_name
        String? in_suffix
        Array[File] in_alignment_bam_chunk_files
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        samtools merge \
          -f -p -c --threads ~{in_map_cores} \
          ~{in_sample_name}_merged.positionsorted~{in_suffix}.bam \
          ~{sep=" " in_alignment_bam_chunk_files} \
        && samtools index \
          ~{in_sample_name}_merged.positionsorted~{in_suffix}.bam
    >>>
    output {
        File merged_bam_file = "~{in_sample_name}_merged.positionsorted~{in_suffix}.bam"
        File merged_bam_file_index = "~{in_sample_name}_merged.positionsorted~{in_suffix}.bam.bai"
    }
    runtime {
        preemptible: 2
        time: 240
        memory: 5 + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task splitBAMbyPath {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        set -eux -o pipefail

        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai

        while read -r contig; do
            samtools view \
              -@ ~{in_map_cores} \
              -h -O BAM \
              input_bam_file.bam ${contig} \
              -o ~{in_sample_name}.${contig}.bam \
            && samtools index \
              ~{in_sample_name}.${contig}.bam
        done < "~{in_path_list_file}"
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai")
    }
    runtime {
        preemptible: 2
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task runGATKRealignerTargetCreator {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_index_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_realign_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command 
        # to exit with a non-zero status, or zero if all commands of the pipeline exit 
        set -o pipefail
        # cause a bash script to exit immediately when a command fails 
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately 
        set -u
        # echo each line of the script to stdout so we can see what is happening 
        set -o xtrace
        #to turn off echo do 'set +o xtrace' 

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s "~{in_reference_file}" reference.fa
        ln -f -s "~{in_reference_index_file}" reference.fa.fai
        # And the dict must be adjacent to both
        ln -f -s "~{in_reference_dict_file}" reference.dict

        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          --remove_program_records \
          -drf DuplicateRead \
          --disable_bam_indexing \
          -nt "32" \
          -R reference.fa \
          -L ${CONTIG_ID} \
          -I input_bam_file.bam \
          --out forIndelRealigner.intervals

        awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{in_sample_name}.${CONTIG_ID}.intervals.bed
    >>>
    output {
        File realigner_target_bed = glob("*.bed")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: 20 + " GB"
        cpu: 16
        disks: "local-disk " + in_realign_disk + " SSD"
        docker: "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b"
    }
}
task widenRealignmentTargets {
    input {
        File in_target_bed_file
        File in_reference_index_file
        Int in_expansion_bases
        Int in_realign_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        BASE_NAME=($(ls ~{in_target_bed_file} | rev | cut -f1 -d'/' | rev | sed s/.bed$//g))

        # Widen the BED regions, but don't escape the chromosomes
        bedtools slop -i "~{in_target_bed_file}" -g "~{in_reference_index_file}" -b "~{in_expansion_bases}" > "${BASE_NAME}.widened.bed"
    >>>
    output {
        File output_target_bed_file = glob("*.widened.bed")[0]
    }
    runtime {
        preemptible: 2
        memory: 4 + " GB"
        cpu: 1
        disks: "local-disk " + in_realign_disk + " SSD"
        docker: "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    }
}
task runAbraRealigner {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_index_file
        File in_target_bed_file
        File in_reference_file
        File in_reference_index_file
        Int in_realign_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai

        java -Xmx20G -jar /opt/abra2/abra2.jar \
          --targets ~{in_target_bed_file} \
          --in input_bam_file.bam \
          --out ~{in_sample_name}.${CONTIG_ID}.indel_realigned.bam \
          --ref reference.fa \
          --index \
          --threads 32
    >>>
    output {
        File indel_realigned_bam = glob("~{in_sample_name}.*.indel_realigned.bam")[0]
        File indel_realigned_bam_index = glob("~{in_sample_name}.*.indel_realigned.bai")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: 20 + " GB"
        cpu: 16
        disks: "local-disk " + in_realign_disk + " SSD"
        # This used to be docker: "dceoy/abra2:latest" but they moved the tag
        # and it stopped working. A known good version has been rehosted on
        # Quay in case Docker Hub deletes it.
        docker: "quay.io/adamnovak/dceoy-abra2@sha256:43d09d1c10220cfeab09e2763c2c5257884fa4457bcaa224f4e3796a28a24bba"
    }
}

task leftShiftBAMFile {
    input {
        String in_sample_name
        File in_bam_file
        File in_reference_file
        File in_reference_index_file
        Int in_realign_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai
        
        bamleftalign \
            <~{in_bam_file} \
            >~{in_sample_name}.${CONTIG_ID}.left_shifted.bam \
            --fasta-reference reference.fa \
            --compressed
    >>>
    output {
        File left_shifted_bam = glob("~{in_sample_name}.*.left_shifted.bam")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: 20 + " GB"
        cpu: 1
        disks: "local-disk " + in_realign_disk + " SSD"
        docker: "biocontainers/freebayes:v1.2.0-2-deb_cv1"
    }
}

