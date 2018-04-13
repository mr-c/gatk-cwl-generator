id: GenotypeConcordance
cwlVersion: v1.0
baseCommand:
- java
- -jar
- /usr/GenomeAnalysisTK.jar
- --analysis_type
- GenotypeConcordance
class: CommandLineTool
doc: |2-


   <p>
    This tool takes in two callsets (vcfs) and tabulates the number of sites which overlap and share alleles,
    and for each sample, the genotype-by-genotype counts (e.g. the number of sites at which a sample was
    called homozygous-reference in the EVAL callset, but homozygous-variant in the COMP callset). It outputs these
    counts as well as convenient proportions (such as the proportion of het calls in the EVAL which were called REF in
    the COMP) and metrics (such as NRD and NRS).
    </p>

    <h3>Input</h3>
    <p>
    Genotype concordance requires two callsets (as it does a comparison): an EVAL and a COMP callset, specified via
    the -eval and -comp arguments. Typically, the EVAL callset is an experimental set you want to evaluate, while the
    COMP callset is a previously existing set used as a standard for comparison (taken to represent "truth").
    </p>
    <p>
    (Optional) Jexl expressions for genotype-level filtering of EVAL or COMP genotypes, specified via the -gfe and
    -cfe arguments, respectively.
    </p>

    <h3>Output</h3>
    <p>
    Genotype Concordance writes a GATK report to the specified file (via -o), consisting of multiple tables of counts
    and proportions. These tables are constructed on a per-sample basis, and include counts of EVAL vs COMP genotype
    states.
    </p>
    <h4>Tables</h4>
    <p>
    Headers for the (non-moltenized -- see below) GenotypeConcordance counts and proportions tables give the genotype of
    the EVAL callset followed by the genotype of the COMP callset. For example the value corresponding to HOM_REF_HET
    reflects variants called HOM_REF in the EVAL callset and HET in the COMP callset. Variants for which the alternate
    alleles between the EVAL and COMP sample did not match are excluded from genotype comparisons and given in the
    "Mismatching_Alleles" field.
    </p>
    <p>
    It may be informative to reshape rows of the GenotypeConcordance counts and proportions tables into separate row-major tables
    where the columns indicate the COMP genotype and the rows indicate the EVAL genotype for easy comparison between the
    two callsets. This can be done with the gsa.reshape.concordance.table function in the gsalib R library.
    In Excel this can be accomplished using the OFFSET function.
    </p>
    <ul>
        <li><i>GenotypeConcordance_CompProportions</i>: Gives the proportion of variants in each category normalized to the total number of called genotypes in the COMP callset</li>
        <li><i>GenotypeConcordance_Counts</i>: Gives the counts for number of genotypes in each category</li>
        <li><i>GenotypeConcordance_EvalProportions</i>: Gives the proportion of genotypes in each category normalized to the total number of called genotypes in the EVAL callset</li>
        <li><i>GenotypeConcordance_Summary</i>: Summary statistics for the sum of all samples and each sample individually. See below for definitions.</li>
        <li><i>SiteConcordance_Summary</i>: Gives comparison counts of called genotypes and their alleles between the two callsets. See below for definitions.</li>
    </ul>
    </p>

    <h4>Term and metrics definitions</h4>
   <p>
   <ul>
            <li><i>GenotypeConcordance_CompProportions</i>, <i>GenotypeConcordance_Counts</i>, and <i>GenotypeConcordance_EvalProportions</i></li>
            <ul>
                 <li>NO_CALL: reported genotype is ./., indicating not enough data to call</li>
                 <li>HET: heterozygous</li>
                 <li>HOM_REF: homozygous reference</li>
                 <li>HOM_VAR: homozygous variant</li>
                 <li>UNAVAILABLE: variant is not called in this callset</li>
                 <li>MIXED: something like ./1</li>
            </ul>
            <li><i>GenotypeConcordance_Summary</i></li>
            <ul>
                 <li>Non-Reference_Sensitivity (NRS): sensitivity of the EVAL calls to polymorphic calls in the COMP set, calculated by (# true positive)/(# true polymorphic)</li>
                 <li>Non-Reference_Discrepancy (NRD): genotype discordance excluding concordant reference sites, calculated by (# discordant sites)/(total excluding # HOM_REF_HOM_REF) = 1.0-(# HOM_VAR_HOM_VAR + # HET_HET)/(total excluding # HOM_REF_HOM_REF)</li>
                 <li>Overall_Genotype_Concordance: overall concordance calculated by (# concordant genotypes)/(# genotypes)</li>
            </ul>
            <li><i>SiteConcordance_Summary</i></li>
            <ul>
                 <li>ALLELES_MATCH: counts of calls at the same site where the alleles match</li>
                 <li>ALLELES_DO_NOT_MATCH: counts of calls at the same location with different alleles, such as the EVAL set calling a 'G' alternate allele, and the comp set calling a 'T' alternate allele</li>
                 <li>EVAL_SUBSET_TRUTH: (multi-alleleic sites only) ALT alleles for EVAL are a subset of ALT alleles for COMP. See also below.</li>
                 <li>EVAL_SUPERSET_TRUTH: (multi-allelic sites only) ALT alleles for COMP are a subset of ALT alleles for EVAL. See also below.</li>
                 <li>EVAL_ONLY: counts of sites present only in the EVAL set, not in the COMP set</li>
                 <li>TRUTH_ONLY: counts of sites present only in the COMP set, not in the EVAL set</li>
            </ul>
   </ul>
   </p>

   <h4>Site-level allelic concordance</h4>

   <p>
    For strictly bi-allelic VCFs, only the ALLELES_MATCH, EVAL_ONLY, TRUTH_ONLY fields will be populated,
    but where multi-allelic sites are involved counts for EVAL_SUBSET_TRUTH and EVAL_SUPERSET_TRUTH will be generated.
   </p>
   <p>
    For example, in the following situation
    <pre>
      eval:  ref - A   alt - C
      comp:  ref - A   alt - C,T
    </pre>
    then the site is tabulated as EVAL_SUBSET_TRUTH. Were the situation reversed, it would be EVAL_SUPERSET_TRUTH.
    However, in the case where EVAL has both C and T alternate alleles, both must be observed in the genotypes
    (that is, there must be at least one of (0/1,1/1) and at least one of (0/2,1/2,2/2) in the genotype field). If
    one of the alleles has no observations in the genotype fields of the EVAL, the site-level concordance is
    tabulated as though that allele were not present in the record.
   </p>

    <h4>Monomorphic Records</h4>
    <p>
    A site which has an alternate allele, but which is monomorphic in samples, is treated as not having been
    discovered, and will be recorded in the TRUTH_ONLY column (if a record exists in the COMP set), or not at all
    (if no record exists in the COMP set).
    </p>
    <p>
    That is, in the situation
    <pre>
     eval:  ref - A   alt - C   genotypes - 0/0  0/0  0/0 ... 0/0
     comp:  ref - A   alt - C   ...         0/0  0/0  ...
    </pre>
    is equivalent to
    <pre>
     eval:  ref - A   alt - .   genotypes - 0/0  0/0  0/0 ... 0/0
     comp:  ref - A   alt - C   ...         0/0  0/0  ...
    </pre>
    </p>
    <p>
    When a record is present in the COMP set the *genotypes* for the monomorphic site will still be used to evaluate
    per-sample genotype concordance counts.
   </p>

    <h4>Filtered Records</h4>
    Filtered records are treated as though they were not present in the VCF, unless -ignoreSiteFilters is provided,
    in which case all records are used. There is currently no way to assess concordance metrics on filtered sites
    exclusively. SelectVariants can be used to extract filtered sites, and VariantFiltration used to un-filter them.

   <h4>Moltenized tables</h4>

   <p>These tables may be optionally moltenized via the -moltenize argument. That is, the standard table

    <pre>
    Sample   NO_CALL_HOM_REF  NO_CALL_HET  NO_CALL_HOM_VAR   (...)
    NA12878       0.003        0.001            0.000        (...)
    NA12891       0.005        0.000            0.000        (...)
    </pre>

    would instead be displayed

   <pre>
    NA12878  NO_CALL_HOM_REF   0.003
    NA12878  NO_CALL_HET       0.001
    NA12878  NO_CALL_HOM_VAR   0.000
    NA12891  NO_CALL_HOM_REF   0.005
    NA12891  NO_CALL_HET       0.000
    NA12891  NO_CALL_HOM_VAR   0.000
    (...)
    </pre>

   <h3>Usage example</h3>
   <pre>
   java -jar GenomeAnalysisTK.jar \
     -T GenotypeConcordance \
     -R reference.fasta \
     -eval test_set.vcf \
     -comp truth_set.vcf \
     -o output.grp
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
- doc: Minimum read mapping quality required to consider a read for calling
  id: min_mapping_quality_score
  type: int?
  inputBinding:
    prefix: --min_mapping_quality_score
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
- doc: The variants and genotypes to compare against
  id: comp
  type: File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--comp", inputs['comp_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'comp'
  id: comp_tags
- doc: The variants and genotypes to compare against
  id: comp
  type: File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--comp", inputs['comp_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'comp'
  id: comp_tags
- doc: The variants and genotypes to evaluate
  id: eval
  type: File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--eval", inputs['eval_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'eval'
  id: eval_tags
- doc: The variants and genotypes to evaluate
  id: eval
  type: File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--eval", inputs['eval_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'eval'
  id: eval_tags
- doc: One or more criteria to use to set COMP genotypes to no-call. These genotype-level
    filters are only applied to the COMP rod.
  id: genotypeFilterExpressionComp
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--genotypeFilterExpressionComp"))
- doc: One or more criteria to use to set COMP genotypes to no-call. These genotype-level
    filters are only applied to the COMP rod.
  id: genotypeFilterExpressionComp
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--genotypeFilterExpressionComp"))
- doc: One or more criteria to use to set EVAL genotypes to no-call. These genotype-level
    filters are only applied to the EVAL rod.
  id: genotypeFilterExpressionEval
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--genotypeFilterExpressionEval"))
- doc: One or more criteria to use to set EVAL genotypes to no-call. These genotype-level
    filters are only applied to the EVAL rod.
  id: genotypeFilterExpressionEval
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--genotypeFilterExpressionEval"))
- doc: Filters will be ignored
  id: ignoreFilters
  type: boolean?
  inputBinding:
    prefix: --ignoreFilters
- doc: Molten rather than tabular output
  id: moltenize
  type: boolean?
  inputBinding:
    prefix: --moltenize
- doc: Molten rather than tabular output
  id: moltenize
  type: boolean?
  inputBinding:
    prefix: --moltenize
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
- doc: File to output the discordant sites and genotypes.
  id: printInterestingSites
  type: string?
  inputBinding:
    prefix: --printInterestingSites
- doc: File to output the discordant sites and genotypes.
  id: printInterestingSites
  type: string?
  inputBinding:
    prefix: --printInterestingSites
outputs: []
