id: CommandLineGATK
cwlVersion: v1.0
baseCommand:
- java
- -jar
- /usr/GenomeAnalysisTK.jar
- --analysis_type
- CommandLineGATK
class: CommandLineTool
doc: |2-


   <h3>Info for end users</h3>

   <p>This is a list of options and parameters that are generally available to all tools in the GATK.</p>

   <p>There may be a few restrictions, which are indicated in individual argument descriptions. For example the -BQSR
   argument is only meant to be used with a subset of tools, and the -pedigree argument will only be effectively used
   by a subset of tools as well. Some arguments conflict with others, and some conversely are dependent on others. This
   is all indicated in the detailed argument descriptions, so be sure to read those in their entirety rather than just
   skimming the one-line summary in the table.</p>

   <h3>Info for developers</h3>

   <p>This class is the GATK engine itself, which manages map/reduce data access and runs walkers.</p>

   <p>We run command line GATK programs using this class. It gets the command line args, parses them, and hands the
   gatk all the parsed out information. Pretty much anything dealing with the underlying system should go here;
   the GATK engine should deal with any data related information.</p>
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
outputs: []
