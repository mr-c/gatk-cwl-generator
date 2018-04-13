id: DepthOfCoverage
cwlVersion: v1.0
baseCommand:
- java
- -jar
- /usr/GenomeAnalysisTK.jar
- --analysis_type
- DepthOfCoverage
class: CommandLineTool
doc: |2-


   <p>
   This tool processes a set of bam files to determine coverage at different levels of partitioning and
   aggregation. Coverage can be analyzed per locus, per interval, per gene, or in total; can be partitioned by
   sample, by read group, by technology, by center, or by library; and can be summarized by mean, median, quartiles,
   and/or percentage of bases covered to or beyond a threshold. Additionally, reads and bases can be filtered by
   mapping or base quality score.
   </p>

   <h3>Input</h3>
   <ul>
       <li>One or more bam files (with proper headers) to be analyzed for coverage statistics</li>
       <li>(Optional) A REFSEQ file to aggregate coverage to the gene level (for information about creating the REFSEQ Rod, please consult the online documentation)</li>
   </ul>

   <h3>Output</h3>
   <p>
   Tables pertaining to different coverage summaries. Suffix on the table files declares the contents:
   </p>
   <ul>
       <li>no suffix: per locus coverage</li>
       <li>_summary: total, mean, median, quartiles, and threshold proportions, aggregated over all bases</li>
       <li>_statistics: coverage histograms (# locus with X coverage), aggregated over all bases</li>
       <li>_interval_summary: total, mean, median, quartiles, and threshold proportions, aggregated per interval</li>
       <li>_interval_statistics: 2x2 table of # of intervals covered to >= X depth in >=Y samples</li>
       <li>_gene_summary: total, mean, median, quartiles, and threshold proportions, aggregated per gene</li>
       <li>_gene_statistics: 2x2 table of # of genes covered to >= X depth in >= Y samples</li>
       <li>_cumulative_coverage_counts: coverage histograms (# locus with >= X coverage), aggregated over all bases</li>
       <li>_cumulative_coverage_proportions: proprotions of loci with >= X coverage, aggregated over all bases</li>
   </ul>

   <h3>Usage example</h3>
   <pre>
   java -jar GenomeAnalysisTK.jar \
     -T DepthOfCoverage \
     -R reference.fasta \
     -o file_name_base \
     -I input_bams.list
     [-geneList refSeq.sorted.txt] \
     [-pt readgroup] \
     [-ct 4 -ct 6 -ct 10] \
     [-L my_capture_genes.interval_list]
   </pre>
requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - |-
    /**
     * File of functions to be added to cwl files
     */

    function generateGATK4BooleanValue(){
        /**
         * Boolean types in GATK 4 are expressed on the command line as --<PREFIX> "true"/"false",
         * so patch here
         */
        if(self === true || self === false){
            return self.toString()
        }

        return self;
    }

    function applyTagsToArgument(prefix, tags){
        /**
         * Function to be used in the field valueFrom of File objects to add gatk tags.
         */

        if(!self){
            return null;
        }
        else if(!tags){
            return generateArrayCmd(prefix);
        }
        else{
            function addTagToArgument(tagObject, argument){
                var allTags = Array.isArray(tagObject) ? tagObject.join(",") : tagObject;

                return [prefix + ":" + allTags, argument];
            }

            if(Array.isArray(self)){
                if(!Array.isArray(tags) || self.length !== tags.length){
                    throw new TypeError("Argument '" + prefix + "' tag field is invalid");
                }

                var value = self.map(function(element, i) {
                    return addTagToArgument(tags[i], element);
                }).reduce(function(a, b){return a.concat(b)})

                return value;
            }
            else{
                return addTagToArgument(tags, self);
            }
        }
    }

    function generateArrayCmd(prefix){
        /**
         * Function to be used in the field valueFrom of array objects, so that arrays are optional
         * and prefixes are handled properly.
         *
         * The issue that this solves is documented here:
         * https://www.biostars.org/p/258414/#260140
         */
        if(!self){
            return null;
        }

        if(!Array.isArray(self)){
            self = [self];
        }

        var output = [];
        self.forEach(function(element) {
            output.push(prefix);
            output.push(element);
        })

        return output;
    }
- class: DockerRequirement
  dockerPull: broadinstitute/gatk3:3.5-0
inputs:
- doc: Reference sequence file
  id: reference_sequence
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--reference_sequence", inputs['reference_sequence_tags']))
  secondaryFiles:
  - .fai
  - ^.dict
- doc: Reference sequence file
  id: reference_sequence
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--reference_sequence", inputs['reference_sequence_tags']))
  secondaryFiles:
  - .fai
  - ^.dict
- doc: Input file containing sequence data (BAM or CRAM)
  id: input_file
  type:
  - 'null'
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--input_file", inputs['input_file_tags']))
  secondaryFiles: $(self.basename + self.nameext.replace('m','i'))
- doc: Input file containing sequence data (BAM or CRAM)
  id: input_file
  type:
  - 'null'
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--input_file", inputs['input_file_tags']))
  secondaryFiles: $(self.basename + self.nameext.replace('m','i'))
- doc: Ignore warnings about base quality score encoding
  id: allow_potentially_misencoded_quality_scores
  type: boolean?
  inputBinding:
    prefix: --allow_potentially_misencoded_quality_scores
- doc: Ignore warnings about base quality score encoding
  id: allow_potentially_misencoded_quality_scores
  type: boolean?
  inputBinding:
    prefix: --allow_potentially_misencoded_quality_scores
- doc: Compression level to use for writing BAM files (0 - 9, higher is more compressed)
  id: bam_compression
  type: int?
  inputBinding:
    prefix: --bam_compression
- doc: Compression level to use for writing BAM files (0 - 9, higher is more compressed)
  id: bam_compression
  type: int?
  inputBinding:
    prefix: --bam_compression
- doc: Type of BAQ calculation to apply in the engine
  id: baq
  type:
  - 'null'
  - type: enum
    symbols:
    - OFF
    - CALCULATE_AS_NECESSARY
    - RECALCULATE
  inputBinding:
    prefix: --baq
- doc: Type of BAQ calculation to apply in the engine
  id: baq
  type:
  - 'null'
  - type: enum
    symbols:
    - OFF
    - CALCULATE_AS_NECESSARY
    - RECALCULATE
  inputBinding:
    prefix: --baq
- doc: BAQ gap open penalty
  id: baqGapOpenPenalty
  type: double?
  inputBinding:
    prefix: --baqGapOpenPenalty
- doc: BAQ gap open penalty
  id: baqGapOpenPenalty
  type: double?
  inputBinding:
    prefix: --baqGapOpenPenalty
- doc: Input covariates table file for on-the-fly base quality score recalibration
  id: BQSR
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--BQSR", inputs['BQSR_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'BQSR'
  id: BQSR_tags
- doc: Input covariates table file for on-the-fly base quality score recalibration
  id: BQSR
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--BQSR", inputs['BQSR_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'BQSR'
  id: BQSR_tags
- doc: Disable both auto-generation of index files and index file locking
  id: disable_auto_index_creation_and_locking_when_reading_rods
  type: boolean?
  inputBinding:
    prefix: --disable_auto_index_creation_and_locking_when_reading_rods
- doc: Disable both auto-generation of index files and index file locking
  id: disable_auto_index_creation_and_locking_when_reading_rods
  type: boolean?
  inputBinding:
    prefix: --disable_auto_index_creation_and_locking_when_reading_rods
- doc: Turn off on-the-fly creation of indices for output BAM/CRAM files.
  id: disable_bam_indexing
  type: boolean?
  inputBinding:
    prefix: --disable_bam_indexing
- doc: Disable printing of base insertion and deletion tags (with -BQSR)
  id: disable_indel_quals
  type: boolean?
  inputBinding:
    prefix: --disable_indel_quals
- doc: Disable printing of base insertion and deletion tags (with -BQSR)
  id: disable_indel_quals
  type: boolean?
  inputBinding:
    prefix: --disable_indel_quals
- doc: Read filters to disable
  id: disable_read_filter
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--disable_read_filter"))
- doc: Read filters to disable
  id: disable_read_filter
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--disable_read_filter"))
- doc: Target coverage threshold for downsampling to coverage
  id: downsample_to_coverage
  type: int?
  inputBinding:
    prefix: --downsample_to_coverage
- doc: Target coverage threshold for downsampling to coverage
  id: downsample_to_coverage
  type: int?
  inputBinding:
    prefix: --downsample_to_coverage
- doc: Fraction of reads to downsample to
  id: downsample_to_fraction
  type: double?
  inputBinding:
    prefix: --downsample_to_fraction
- doc: Fraction of reads to downsample to
  id: downsample_to_fraction
  type: double?
  inputBinding:
    prefix: --downsample_to_fraction
- doc: Type of read downsampling to employ at a given locus
  id: downsampling_type
  type:
  - 'null'
  - type: enum
    symbols:
    - NONE
    - ALL_READS
    - BY_SAMPLE
  inputBinding:
    prefix: --downsampling_type
- doc: Type of read downsampling to employ at a given locus
  id: downsampling_type
  type:
  - 'null'
  - type: enum
    symbols:
    - NONE
    - ALL_READS
    - BY_SAMPLE
  inputBinding:
    prefix: --downsampling_type
- doc: Emit the OQ tag with the original base qualities (with -BQSR)
  id: emit_original_quals
  type: boolean?
  inputBinding:
    prefix: --emit_original_quals
- doc: Emit the OQ tag with the original base qualities (with -BQSR)
  id: emit_original_quals
  type: boolean?
  inputBinding:
    prefix: --emit_original_quals
- doc: One or more genomic intervals to exclude from processing
  id: excludeIntervals
  type:
  - 'null'
  - type: array
    items:
    - File
    - string
    inputBinding:
      valueFrom: $(null)
  - File
  - string
  inputBinding:
    valueFrom: $(applyTagsToArgument("--excludeIntervals", inputs['excludeIntervals_tags']))
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'excludeIntervals'
  id: excludeIntervals_tags
- doc: One or more genomic intervals to exclude from processing
  id: excludeIntervals
  type:
  - 'null'
  - type: array
    items:
    - File
    - string
    inputBinding:
      valueFrom: $(null)
  - File
  - string
  inputBinding:
    valueFrom: $(applyTagsToArgument("--excludeIntervals", inputs['excludeIntervals_tags']))
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'excludeIntervals'
  id: excludeIntervals_tags
- doc: Fix mis-encoded base quality scores
  id: fix_misencoded_quality_scores
  type: boolean?
  inputBinding:
    prefix: --fix_misencoded_quality_scores
- doc: Fix mis-encoded base quality scores
  id: fix_misencoded_quality_scores
  type: boolean?
  inputBinding:
    prefix: --fix_misencoded_quality_scores
- doc: GATK key file required to run with -et NO_ET
  id: gatk_key
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--gatk_key", inputs['gatk_key_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'gatk_key'
  id: gatk_key_tags
- doc: GATK key file required to run with -et NO_ET
  id: gatk_key
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--gatk_key", inputs['gatk_key_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'gatk_key'
  id: gatk_key_tags
- doc: Enable on-the-fly creation of md5s for output BAM files.
  id: generate_md5
  type: boolean?
  inputBinding:
    prefix: --generate_md5
- doc: Global Qscore Bayesian prior to use for BQSR
  id: globalQScorePrior
  type: double?
  inputBinding:
    prefix: --globalQScorePrior
- doc: Global Qscore Bayesian prior to use for BQSR
  id: globalQScorePrior
  type: double?
  inputBinding:
    prefix: --globalQScorePrior
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'input_file'
  id: input_file_tags
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'input_file'
  id: input_file_tags
- doc: Interval merging rule for abutting intervals
  id: interval_merging
  type:
  - 'null'
  - type: enum
    symbols:
    - ALL
    - OVERLAPPING_ONLY
  inputBinding:
    prefix: --interval_merging
- doc: Interval merging rule for abutting intervals
  id: interval_merging
  type:
  - 'null'
  - type: enum
    symbols:
    - ALL
    - OVERLAPPING_ONLY
  inputBinding:
    prefix: --interval_merging
- doc: Amount of padding (in bp) to add to each interval
  id: interval_padding
  type: int?
  inputBinding:
    prefix: --interval_padding
- doc: Amount of padding (in bp) to add to each interval
  id: interval_padding
  type: int?
  inputBinding:
    prefix: --interval_padding
- doc: Set merging approach to use for combining interval inputs
  id: interval_set_rule
  type:
  - 'null'
  - type: enum
    symbols:
    - UNION
    - INTERSECTION
  inputBinding:
    prefix: --interval_set_rule
- doc: Set merging approach to use for combining interval inputs
  id: interval_set_rule
  type:
  - 'null'
  - type: enum
    symbols:
    - UNION
    - INTERSECTION
  inputBinding:
    prefix: --interval_set_rule
- doc: One or more genomic intervals over which to operate
  id: intervals
  type:
  - 'null'
  - type: array
    items:
    - File
    - string
    inputBinding:
      valueFrom: $(null)
  - File
  - string
  inputBinding:
    valueFrom: $(applyTagsToArgument("--intervals", inputs['intervals_tags']))
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'intervals'
  id: intervals_tags
- doc: One or more genomic intervals over which to operate
  id: intervals
  type:
  - 'null'
  - type: array
    items:
    - File
    - string
    inputBinding:
      valueFrom: $(null)
  - File
  - string
  inputBinding:
    valueFrom: $(applyTagsToArgument("--intervals", inputs['intervals_tags']))
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'intervals'
  id: intervals_tags
- doc: Keep program records in the SAM header
  id: keep_program_records
  type: boolean?
  inputBinding:
    prefix: --keep_program_records
- doc: Keep program records in the SAM header
  id: keep_program_records
  type: boolean?
  inputBinding:
    prefix: --keep_program_records
- doc: Set the logging location
  id: log_to_file
  type: string?
  inputBinding:
    prefix: --log_to_file
- doc: Set the logging location
  id: log_to_file
  type: string?
  inputBinding:
    prefix: --log_to_file
- doc: Set the minimum level of logging
  id: logging_level
  type: string?
  inputBinding:
    prefix: --logging_level
- doc: Set the minimum level of logging
  id: logging_level
  type: string?
  inputBinding:
    prefix: --logging_level
- doc: Stop execution cleanly as soon as maxRuntime has been reached
  id: maxRuntime
  type: long?
  inputBinding:
    prefix: --maxRuntime
- doc: Stop execution cleanly as soon as maxRuntime has been reached
  id: maxRuntime
  type: long?
  inputBinding:
    prefix: --maxRuntime
- doc: Unit of time used by maxRuntime
  id: maxRuntimeUnits
  type:
  - 'null'
  - type: enum
    symbols:
    - NANOSECONDS
    - MICROSECONDS
    - MILLISECONDS
    - SECONDS
    - MINUTES
    - HOURS
    - DAYS
  inputBinding:
    prefix: --maxRuntimeUnits
- doc: Unit of time used by maxRuntime
  id: maxRuntimeUnits
  type:
  - 'null'
  - type: enum
    symbols:
    - NANOSECONDS
    - MICROSECONDS
    - MILLISECONDS
    - SECONDS
    - MINUTES
    - HOURS
    - DAYS
  inputBinding:
    prefix: --maxRuntimeUnits
- doc: Enable threading efficiency monitoring
  id: monitorThreadEfficiency
  type: boolean?
  inputBinding:
    prefix: --monitorThreadEfficiency
- doc: Enable threading efficiency monitoring
  id: monitorThreadEfficiency
  type: boolean?
  inputBinding:
    prefix: --monitorThreadEfficiency
- doc: Always output all the records in VCF FORMAT fields, even if some are missing
  id: never_trim_vcf_format_field
  type: boolean?
  inputBinding:
    prefix: --never_trim_vcf_format_field
- doc: Always output all the records in VCF FORMAT fields, even if some are missing
  id: never_trim_vcf_format_field
  type: boolean?
  inputBinding:
    prefix: --never_trim_vcf_format_field
- doc: Use a non-deterministic random seed
  id: nonDeterministicRandomSeed
  type: boolean?
  inputBinding:
    prefix: --nonDeterministicRandomSeed
- doc: Use a non-deterministic random seed
  id: nonDeterministicRandomSeed
  type: boolean?
  inputBinding:
    prefix: --nonDeterministicRandomSeed
- doc: Number of CPU threads to allocate per data thread
  id: num_cpu_threads_per_data_thread
  type: int?
  inputBinding:
    prefix: --num_cpu_threads_per_data_thread
- doc: Number of CPU threads to allocate per data thread
  id: num_cpu_threads_per_data_thread
  type: int?
  inputBinding:
    prefix: --num_cpu_threads_per_data_thread
- doc: Number of data threads to allocate to this analysis
  id: num_threads
  type: int?
  inputBinding:
    prefix: --num_threads
- doc: Number of data threads to allocate to this analysis
  id: num_threads
  type: int?
  inputBinding:
    prefix: --num_threads
- doc: Pedigree files for samples
  id: pedigree
  type:
  - 'null'
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--pedigree", inputs['pedigree_tags']))
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'pedigree'
  id: pedigree_tags
- doc: Pedigree files for samples
  id: pedigree
  type:
  - 'null'
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--pedigree", inputs['pedigree_tags']))
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'pedigree'
  id: pedigree_tags
- doc: Pedigree string for samples
  id: pedigreeString
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--pedigreeString"))
- doc: Pedigree string for samples
  id: pedigreeString
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--pedigreeString"))
- doc: Validation strictness for pedigree information
  id: pedigreeValidationType
  type:
  - 'null'
  - type: enum
    symbols:
    - STRICT
    - SILENT
  inputBinding:
    prefix: --pedigreeValidationType
- doc: Validation strictness for pedigree information
  id: pedigreeValidationType
  type:
  - 'null'
  - type: enum
    symbols:
    - STRICT
    - SILENT
  inputBinding:
    prefix: --pedigreeValidationType
- doc: Write GATK runtime performance log to this file
  id: performanceLog
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--performanceLog", inputs['performanceLog_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'performanceLog'
  id: performanceLog_tags
- doc: Write GATK runtime performance log to this file
  id: performanceLog
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--performanceLog", inputs['performanceLog_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'performanceLog'
  id: performanceLog_tags
- doc: Run reporting mode
  id: phone_home
  type:
  - 'null'
  - type: enum
    symbols:
    - NO_ET
    - AWS
    - STDOUT
  inputBinding:
    prefix: --phone_home
- doc: Run reporting mode
  id: phone_home
  type:
  - 'null'
  - type: enum
    symbols:
    - NO_ET
    - AWS
    - STDOUT
  inputBinding:
    prefix: --phone_home
- doc: Don't recalibrate bases with quality scores less than this threshold (with
    -BQSR)
  id: preserve_qscores_less_than
  type: int?
  inputBinding:
    prefix: --preserve_qscores_less_than
- doc: Don't recalibrate bases with quality scores less than this threshold (with
    -BQSR)
  id: preserve_qscores_less_than
  type: int?
  inputBinding:
    prefix: --preserve_qscores_less_than
- doc: Quantize quality scores to a given number of levels (with -BQSR)
  id: quantize_quals
  type: int?
  inputBinding:
    prefix: --quantize_quals
- doc: Quantize quality scores to a given number of levels (with -BQSR)
  id: quantize_quals
  type: int?
  inputBinding:
    prefix: --quantize_quals
- doc: Number of reads per SAM file to buffer in memory
  id: read_buffer_size
  type: int?
  inputBinding:
    prefix: --read_buffer_size
- doc: Number of reads per SAM file to buffer in memory
  id: read_buffer_size
  type: int?
  inputBinding:
    prefix: --read_buffer_size
- doc: Filters to apply to reads before analysis
  id: read_filter
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--read_filter"))
- doc: Filters to apply to reads before analysis
  id: read_filter
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--read_filter"))
- doc: Exclude read groups based on tags
  id: read_group_black_list
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--read_group_black_list"))
- doc: Exclude read groups based on tags
  id: read_group_black_list
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--read_group_black_list"))
- doc: Reduce NDN elements in CIGAR string
  id: refactor_NDN_cigar_string
  type: boolean?
  inputBinding:
    prefix: --refactor_NDN_cigar_string
- doc: Reduce NDN elements in CIGAR string
  id: refactor_NDN_cigar_string
  type: boolean?
  inputBinding:
    prefix: --refactor_NDN_cigar_string
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'reference_sequence'
  id: reference_sequence_tags
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'reference_sequence'
  id: reference_sequence_tags
- doc: Reference window stop
  id: reference_window_stop
  type: int?
  inputBinding:
    prefix: --reference_window_stop
- doc: Reference window stop
  id: reference_window_stop
  type: int?
  inputBinding:
    prefix: --reference_window_stop
- doc: Remove program records from the SAM header
  id: remove_program_records
  type: boolean?
  inputBinding:
    prefix: --remove_program_records
- doc: Remove program records from the SAM header
  id: remove_program_records
  type: boolean?
  inputBinding:
    prefix: --remove_program_records
- doc: Rename sample IDs on-the-fly at runtime using the provided mapping file
  id: sample_rename_mapping_file
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--sample_rename_mapping_file", inputs['sample_rename_mapping_file_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'sample_rename_mapping_file'
  id: sample_rename_mapping_file_tags
- doc: Rename sample IDs on-the-fly at runtime using the provided mapping file
  id: sample_rename_mapping_file
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--sample_rename_mapping_file", inputs['sample_rename_mapping_file_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'sample_rename_mapping_file'
  id: sample_rename_mapping_file_tags
- doc: Emit a log entry (level INFO) containing the full list of sequence data files
    to be included in the analysis (including files inside .bam.list or .cram.list
    files).
  id: showFullBamList
  type: boolean?
  inputBinding:
    prefix: --showFullBamList
- doc: If provided, output BAM/CRAM files will be simplified to include just key reads
    for downstream variation discovery analyses (removing duplicates, PF-, non-primary
    reads), as well stripping all extended tags from the kept reads except the read
    group identifier
  id: simplifyBAM
  type: boolean?
  inputBinding:
    prefix: --simplifyBAM
- doc: If provided, output BAM/CRAM files will be simplified to include just key reads
    for downstream variation discovery analyses (removing duplicates, PF-, non-primary
    reads), as well stripping all extended tags from the kept reads except the read
    group identifier
  id: simplifyBAM
  type: boolean?
  inputBinding:
    prefix: --simplifyBAM
- doc: Just output sites without genotypes (i.e. only the first 8 columns of the VCF)
  id: sites_only
  type: boolean?
  inputBinding:
    prefix: --sites_only
- doc: Just output sites without genotypes (i.e. only the first 8 columns of the VCF)
  id: sites_only
  type: boolean?
  inputBinding:
    prefix: --sites_only
- doc: Use static quantized quality scores to a given number of levels (with -BQSR)
  id: static_quantized_quals
  type:
  - 'null'
  - type: array
    items: int
    inputBinding:
      valueFrom: $(null)
  - int
  inputBinding:
    valueFrom: $(generateArrayCmd("--static_quantized_quals"))
- doc: Use static quantized quality scores to a given number of levels (with -BQSR)
  id: static_quantized_quals
  type:
  - 'null'
  - type: array
    items: int
    inputBinding:
      valueFrom: $(null)
  - int
  inputBinding:
    valueFrom: $(generateArrayCmd("--static_quantized_quals"))
- doc: Tag to identify this GATK run as part of a group of runs
  id: tag
  type: string?
  inputBinding:
    prefix: --tag
- doc: Tag to identify this GATK run as part of a group of runs
  id: tag
  type: string?
  inputBinding:
    prefix: --tag
- doc: 'Enable unsafe operations: nothing will be checked at runtime'
  id: unsafe
  type:
  - 'null'
  - type: enum
    symbols:
    - ALLOW_N_CIGAR_READS
    - ALLOW_UNINDEXED_BAM
    - ALLOW_UNSET_BAM_SORT_ORDER
    - NO_READ_ORDER_VERIFICATION
    - ALLOW_SEQ_DICT_INCOMPATIBILITY
    - LENIENT_VCF_PROCESSING
    - ALL
  inputBinding:
    prefix: --unsafe
- doc: 'Enable unsafe operations: nothing will be checked at runtime'
  id: unsafe
  type:
  - 'null'
  - type: enum
    symbols:
    - ALLOW_N_CIGAR_READS
    - ALLOW_UNINDEXED_BAM
    - ALLOW_UNSET_BAM_SORT_ORDER
    - NO_READ_ORDER_VERIFICATION
    - ALLOW_SEQ_DICT_INCOMPATIBILITY
    - LENIENT_VCF_PROCESSING
    - ALL
  inputBinding:
    prefix: --unsafe
- doc: Use the base quality scores from the OQ tag
  id: useOriginalQualities
  type: boolean?
  inputBinding:
    prefix: --useOriginalQualities
- doc: Use the base quality scores from the OQ tag
  id: useOriginalQualities
  type: boolean?
  inputBinding:
    prefix: --useOriginalQualities
- doc: How strict should we be with validation
  id: validation_strictness
  type:
  - 'null'
  - type: enum
    symbols:
    - STRICT
    - LENIENT
    - SILENT
  inputBinding:
    prefix: --validation_strictness
- doc: How strict should we be with validation
  id: validation_strictness
  type:
  - 'null'
  - type: enum
    symbols:
    - STRICT
    - LENIENT
    - SILENT
  inputBinding:
    prefix: --validation_strictness
- doc: Parameter to pass to the VCF/BCF IndexCreator
  id: variant_index_parameter
  type: int?
  inputBinding:
    prefix: --variant_index_parameter
- doc: Parameter to pass to the VCF/BCF IndexCreator
  id: variant_index_parameter
  type: int?
  inputBinding:
    prefix: --variant_index_parameter
- doc: Type of IndexCreator to use for VCF/BCF indices
  id: variant_index_type
  type:
  - 'null'
  - type: enum
    symbols:
    - DYNAMIC_SEEK
    - DYNAMIC_SIZE
    - LINEAR
    - INTERVAL
  inputBinding:
    prefix: --variant_index_type
- doc: Type of IndexCreator to use for VCF/BCF indices
  id: variant_index_type
  type:
  - 'null'
  - type: enum
    symbols:
    - DYNAMIC_SEEK
    - DYNAMIC_SIZE
    - LINEAR
    - INTERVAL
  inputBinding:
    prefix: --variant_index_type
- doc: Output version information
  id: version
  type: boolean?
  inputBinding:
    prefix: --version
- doc: Output version information
  id: version
  type: boolean?
  inputBinding:
    prefix: --version
- doc: The name of the library to keep, filtering out all others
  id: library
  type: string?
  inputBinding:
    prefix: --library
- doc: The name of the library to keep, filtering out all others
  id: library
  type: string?
  inputBinding:
    prefix: --library
- doc: Filter out reads with no stored bases (i.e. '*' where the sequence should be),
    instead of failing with an error
  id: filter_bases_not_stored
  type: boolean?
  inputBinding:
    prefix: --filter_bases_not_stored
- doc: Filter out reads with no stored bases (i.e. '*' where the sequence should be),
    instead of failing with an error
  id: filter_bases_not_stored
  type: boolean?
  inputBinding:
    prefix: --filter_bases_not_stored
- doc: Filter out reads with mismatching numbers of bases and base qualities, instead
    of failing with an error
  id: filter_mismatching_base_and_quals
  type: boolean?
  inputBinding:
    prefix: --filter_mismatching_base_and_quals
- doc: Filter out reads with mismatching numbers of bases and base qualities, instead
    of failing with an error
  id: filter_mismatching_base_and_quals
  type: boolean?
  inputBinding:
    prefix: --filter_mismatching_base_and_quals
- doc: Filter out reads with CIGAR containing the N operator, instead of failing with
    an error
  id: filter_reads_with_N_cigar
  type: boolean?
  inputBinding:
    prefix: --filter_reads_with_N_cigar
- doc: Filter out reads with CIGAR containing the N operator, instead of failing with
    an error
  id: filter_reads_with_N_cigar
  type: boolean?
  inputBinding:
    prefix: --filter_reads_with_N_cigar
- doc: Minimum read mapping quality required to consider a read for calling
  id: min_mapping_quality_score
  type: int?
  inputBinding:
    prefix: --min_mapping_quality_score
- doc: Minimum mapping quality of reads to count towards depth
  id: minMappingQuality
  type: int?
  inputBinding:
    prefix: --minMappingQuality
- doc: Insert size cutoff
  id: maxInsertSize
  type: int?
  inputBinding:
    prefix: --maxInsertSize
- doc: Insert size cutoff
  id: maxInsertSize
  type: int?
  inputBinding:
    prefix: --maxInsertSize
- doc: Allow a read to be filtered out based on having only 1 soft-clipped block.
    By default, both ends must have a soft-clipped block, setting this flag requires
    only 1 soft-clipped block.
  id: do_not_require_softclips_both_ends
  type: boolean?
  inputBinding:
    prefix: --do_not_require_softclips_both_ends
- doc: Allow a read to be filtered out based on having only 1 soft-clipped block.
    By default, both ends must have a soft-clipped block, setting this flag requires
    only 1 soft-clipped block.
  id: do_not_require_softclips_both_ends
  type: boolean?
  inputBinding:
    prefix: --do_not_require_softclips_both_ends
- doc: Value for which reads with less than this number of aligned bases is considered
    too short
  id: filter_is_too_short_value
  type: int?
  inputBinding:
    prefix: --filter_is_too_short_value
- doc: Value for which reads with less than this number of aligned bases is considered
    too short
  id: filter_is_too_short_value
  type: int?
  inputBinding:
    prefix: --filter_is_too_short_value
- doc: Discard reads with RG:PL attribute containing this string
  id: PLFilterName
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--PLFilterName"))
- doc: Discard reads with RG:PL attribute containing this string
  id: PLFilterName
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--PLFilterName"))
- doc: Discard reads with length greater than the specified value
  id: maxReadLength
  type: int?
  inputBinding:
    prefix: --maxReadLength
- doc: Discard reads with length greater than the specified value
  id: maxReadLength
  type: int?
  inputBinding:
    prefix: --maxReadLength
- doc: Discard reads with length shorter than the specified value
  id: minReadLength
  type: int?
  inputBinding:
    prefix: --minReadLength
- doc: Discard reads with length shorter than the specified value
  id: minReadLength
  type: int?
  inputBinding:
    prefix: --minReadLength
- doc: Read name to whitelist
  id: readName
  type: string?
  inputBinding:
    prefix: --readName
- doc: Read name to whitelist
  id: readName
  type: string?
  inputBinding:
    prefix: --readName
- doc: Discard reads on the forward strand
  id: filterPositive
  type: boolean?
  inputBinding:
    prefix: --filterPositive
- doc: Discard reads on the forward strand
  id: filterPositive
  type: boolean?
  inputBinding:
    prefix: --filterPositive
- doc: Default read mapping quality to assign to all reads
  id: default_mapping_quality
  type: int?
  inputBinding:
    prefix: --default_mapping_quality
- doc: Default read mapping quality to assign to all reads
  id: default_mapping_quality
  type: int?
  inputBinding:
    prefix: --default_mapping_quality
- doc: Original mapping quality
  id: reassign_mapping_quality_from
  type: int?
  inputBinding:
    prefix: --reassign_mapping_quality_from
- doc: Original mapping quality
  id: reassign_mapping_quality_from
  type: int?
  inputBinding:
    prefix: --reassign_mapping_quality_from
- doc: Desired mapping quality
  id: reassign_mapping_quality_to
  type: int?
  inputBinding:
    prefix: --reassign_mapping_quality_to
- doc: Desired mapping quality
  id: reassign_mapping_quality_to
  type: int?
  inputBinding:
    prefix: --reassign_mapping_quality_to
- doc: The name of the sample(s) to keep, filtering out all others
  id: sample_to_keep
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--sample_to_keep"))
- doc: The name of the sample(s) to keep, filtering out all others
  id: sample_to_keep
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--sample_to_keep"))
- doc: The name of the read group to keep, filtering out all others
  id: read_group_to_keep
  type: string?
  inputBinding:
    prefix: --read_group_to_keep
- doc: The name of the read group to keep, filtering out all others
  id: read_group_to_keep
  type: string?
  inputBinding:
    prefix: --read_group_to_keep
- doc: Calculate coverage statistics over this list of genes
  id: calculateCoverageOverGenes
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--calculateCoverageOverGenes", inputs['calculateCoverageOverGenes_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'calculateCoverageOverGenes'
  id: calculateCoverageOverGenes_tags
- doc: Calculate coverage statistics over this list of genes
  id: calculateCoverageOverGenes
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--calculateCoverageOverGenes", inputs['calculateCoverageOverGenes_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'calculateCoverageOverGenes'
  id: calculateCoverageOverGenes_tags
- doc: How should overlapping reads from the same fragment be handled?
  id: countType
  type:
  - 'null'
  - type: enum
    symbols:
    - COUNT_READS
    - COUNT_FRAGMENTS
    - COUNT_FRAGMENTS_REQUIRE_SAME_BASE
  inputBinding:
    prefix: --countType
- doc: Ignore sites consisting only of deletions
  id: ignoreDeletionSites
  type: boolean?
  inputBinding:
    prefix: --ignoreDeletionSites
- doc: Include information on deletions
  id: includeDeletions
  type: boolean?
  inputBinding:
    prefix: --includeDeletions
- doc: Include information on deletions
  id: includeDeletions
  type: boolean?
  inputBinding:
    prefix: --includeDeletions
- doc: Include sites where the reference is N
  id: includeRefNSites
  type: boolean?
  inputBinding:
    prefix: --includeRefNSites
- doc: Maximum quality of bases to count towards depth
  id: maxBaseQuality
  type: int?
  inputBinding:
    prefix: --maxBaseQuality
- doc: Maximum mapping quality of reads to count towards depth
  id: maxMappingQuality
  type: int?
  inputBinding:
    prefix: --maxMappingQuality
- doc: Minimum quality of bases to count towards depth
  id: minBaseQuality
  type: int?
  inputBinding:
    prefix: --minBaseQuality
- doc: Minimum quality of bases to count towards depth
  id: minBaseQuality
  type: int?
  inputBinding:
    prefix: --minBaseQuality
- doc: Minimum mapping quality of reads to count towards depth
  id: minMappingQuality
  type: int?
  inputBinding:
    prefix: --minMappingQuality
- doc: Number of bins to use for granular binning
  id: nBins
  type: int?
  inputBinding:
    prefix: --nBins
- doc: Do not output depth of coverage at each base
  id: omitDepthOutputAtEachBase
  type: boolean?
  inputBinding:
    prefix: --omitDepthOutputAtEachBase
- doc: Do not output depth of coverage at each base
  id: omitDepthOutputAtEachBase
  type: boolean?
  inputBinding:
    prefix: --omitDepthOutputAtEachBase
- doc: Do not calculate per-interval statistics
  id: omitIntervalStatistics
  type: boolean?
  inputBinding:
    prefix: --omitIntervalStatistics
- doc: Do not calculate per-interval statistics
  id: omitIntervalStatistics
  type: boolean?
  inputBinding:
    prefix: --omitIntervalStatistics
- doc: Do not calculate per-sample per-depth counts of loci
  id: omitLocusTable
  type: boolean?
  inputBinding:
    prefix: --omitLocusTable
- doc: Do not calculate per-sample per-depth counts of loci
  id: omitLocusTable
  type: boolean?
  inputBinding:
    prefix: --omitLocusTable
- doc: Do not output the summary files per-sample
  id: omitPerSampleStats
  type: boolean?
  inputBinding:
    prefix: --omitPerSampleStats
- doc: Do not output the summary files per-sample
  id: omitPerSampleStats
  type: boolean?
  inputBinding:
    prefix: --omitPerSampleStats
- doc: An output file created by the walker.  Will overwrite contents if file exists
  id: out
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--out", inputs['out_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'out'
  id: out_tags
- doc: An output file created by the walker.  Will overwrite contents if file exists
  id: out
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--out", inputs['out_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'out'
  id: out_tags
- doc: The format of the output file
  id: outputFormat
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--outputFormat", inputs['outputFormat_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'outputFormat'
  id: outputFormat_tags
- doc: Partition type for depth of coverage
  id: partitionType
  type:
  - 'null'
  - type: array
    items:
      type: enum
      symbols:
      - readgroup
      - sample
      - library
      - platform
      - center
      - sample_by_platform
      - sample_by_center
      - sample_by_platform_by_center
    inputBinding:
      valueFrom: $(null)
  - type: enum
    symbols:
    - readgroup
    - sample
    - library
    - platform
    - center
    - sample_by_platform
    - sample_by_center
    - sample_by_platform_by_center
  inputBinding:
    valueFrom: $(generateArrayCmd("--partitionType"))
- doc: Partition type for depth of coverage
  id: partitionType
  type:
  - 'null'
  - type: array
    items:
      type: enum
      symbols:
      - readgroup
      - sample
      - library
      - platform
      - center
      - sample_by_platform
      - sample_by_center
      - sample_by_platform_by_center
    inputBinding:
      valueFrom: $(null)
  - type: enum
    symbols:
    - readgroup
    - sample
    - library
    - platform
    - center
    - sample_by_platform
    - sample_by_center
    - sample_by_platform_by_center
  inputBinding:
    valueFrom: $(generateArrayCmd("--partitionType"))
- doc: Add base counts to per-locus output
  id: printBaseCounts
  type: boolean?
  inputBinding:
    prefix: --printBaseCounts
- doc: Add base counts to per-locus output
  id: printBaseCounts
  type: boolean?
  inputBinding:
    prefix: --printBaseCounts
- doc: Print the bin values and exit immediately
  id: printBinEndpointsAndExit
  type: boolean?
  inputBinding:
    prefix: --printBinEndpointsAndExit
- doc: Starting (left endpoint) for granular binning
  id: start
  type: int?
  inputBinding:
    prefix: --start
- doc: Ending (right endpoint) for granular binning
  id: stop
  type: int?
  inputBinding:
    prefix: --stop
- doc: Coverage threshold (in percent) for summarizing statistics
  id: summaryCoverageThreshold
  type:
  - 'null'
  - type: array
    items: int
    inputBinding:
      valueFrom: $(null)
  - int
  inputBinding:
    valueFrom: $(generateArrayCmd("--summaryCoverageThreshold"))
- doc: Coverage threshold (in percent) for summarizing statistics
  id: summaryCoverageThreshold
  type:
  - 'null'
  - type: array
    items: int
    inputBinding:
      valueFrom: $(null)
  - int
  inputBinding:
    valueFrom: $(generateArrayCmd("--summaryCoverageThreshold"))
outputs:
- id: library_summaryOutput
  doc: The library_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_summary')
- id: library_statisticsOutput
  doc: The library_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_statistics')
- id: library_interval_summaryOutput
  doc: The library_interval_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_interval_summary')
- id: library_interval_statisticsOutput
  doc: The library_interval_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_interval_statistics')
- id: library_gene_summaryOutput
  doc: The library_gene_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_gene_summary')
- id: library_gene_statisticsOutput
  doc: The library_gene_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_gene_statistics')
- id: library_cumulative_coverage_countsOutput
  doc: The library_cumulative_coverage_counts generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_cumulative_coverage_counts')
- id: library_cumulative_coverage_proportionsOutput
  doc: The library_cumulative_coverage_proportions generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_cumulative_coverage_proportions')
- id: read_group_summaryOutput
  doc: The read_group_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_summary')
- id: read_group_statisticsOutput
  doc: The read_group_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_statistics')
- id: read_group_interval_summaryOutput
  doc: The read_group_interval_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_interval_summary')
- id: read_group_interval_statisticsOutput
  doc: The read_group_interval_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_interval_statistics')
- id: read_group_gene_summaryOutput
  doc: The read_group_gene_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_gene_summary')
- id: read_group_gene_statisticsOutput
  doc: The read_group_gene_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_gene_statistics')
- id: read_group_cumulative_coverage_countsOutput
  doc: The read_group_cumulative_coverage_counts generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_cumulative_coverage_counts')
- id: read_group_cumulative_coverage_proportionsOutput
  doc: The read_group_cumulative_coverage_proportions generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_cumulative_coverage_proportions')
- id: sample_summaryOutput
  doc: The sample_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_summary')
- id: sample_statisticsOutput
  doc: The sample_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_statistics')
- id: sample_interval_summaryOutput
  doc: The sample_interval_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_interval_summary')
- id: sample_interval_statisticsOutput
  doc: The sample_interval_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_interval_statistics')
- id: sample_gene_summaryOutput
  doc: The sample_gene_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_gene_summary')
- id: sample_gene_statisticsOutput
  doc: The sample_gene_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_gene_statistics')
- id: sample_cumulative_coverage_countsOutput
  doc: The sample_cumulative_coverage_counts generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_cumulative_coverage_counts')
- id: sample_cumulative_coverage_proportionsOutput
  doc: The sample_cumulative_coverage_proportions generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_cumulative_coverage_proportions')
- id: library_summaryOutput
  doc: The library_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_summary')
- id: library_statisticsOutput
  doc: The library_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_statistics')
- id: library_interval_summaryOutput
  doc: The library_interval_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_interval_summary')
- id: library_interval_statisticsOutput
  doc: The library_interval_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_interval_statistics')
- id: library_gene_summaryOutput
  doc: The library_gene_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_gene_summary')
- id: library_gene_statisticsOutput
  doc: The library_gene_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_gene_statistics')
- id: library_cumulative_coverage_countsOutput
  doc: The library_cumulative_coverage_counts generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_cumulative_coverage_counts')
- id: library_cumulative_coverage_proportionsOutput
  doc: The library_cumulative_coverage_proportions generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'library_cumulative_coverage_proportions')
- id: read_group_summaryOutput
  doc: The read_group_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_summary')
- id: read_group_statisticsOutput
  doc: The read_group_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_statistics')
- id: read_group_interval_summaryOutput
  doc: The read_group_interval_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_interval_summary')
- id: read_group_interval_statisticsOutput
  doc: The read_group_interval_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_interval_statistics')
- id: read_group_gene_summaryOutput
  doc: The read_group_gene_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_gene_summary')
- id: read_group_gene_statisticsOutput
  doc: The read_group_gene_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_gene_statistics')
- id: read_group_cumulative_coverage_countsOutput
  doc: The read_group_cumulative_coverage_counts generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_cumulative_coverage_counts')
- id: read_group_cumulative_coverage_proportionsOutput
  doc: The read_group_cumulative_coverage_proportions generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'read_group_cumulative_coverage_proportions')
- id: sample_summaryOutput
  doc: The sample_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_summary')
- id: sample_statisticsOutput
  doc: The sample_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_statistics')
- id: sample_interval_summaryOutput
  doc: The sample_interval_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_interval_summary')
- id: sample_interval_statisticsOutput
  doc: The sample_interval_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_interval_statistics')
- id: sample_gene_summaryOutput
  doc: The sample_gene_summary generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_gene_summary')
- id: sample_gene_statisticsOutput
  doc: The sample_gene_statistics generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_gene_statistics')
- id: sample_cumulative_coverage_countsOutput
  doc: The sample_cumulative_coverage_counts generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_cumulative_coverage_counts')
- id: sample_cumulative_coverage_proportionsOutput
  doc: The sample_cumulative_coverage_proportions generated file
  type: File?
  outputBinding:
    glob: $(inputs.out + 'sample_cumulative_coverage_proportions')
