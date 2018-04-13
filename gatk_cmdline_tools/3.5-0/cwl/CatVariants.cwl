id: CatVariants
cwlVersion: v1.0
baseCommand:
- java
- -jar
- /usr/GenomeAnalysisTK.jar
- --analysis_type
- CatVariants
class: CommandLineTool
doc: |2-


   <p>
   The main purpose of this tool is to speed up the gather function when using scatter-gather parallelization.
   This tool concatenates the scattered output VCF files. It assumes that:
   <ul>
       <li>All the input VCFs (or BCFs) contain the same samples in the same order.</li>
       <li>The variants in each input file are from non-overlapping (scattered) intervals.</li>
   </ul>
   </p>
   <p>When the input files are already sorted based on the intervals start positions, use -assumeSorted.</p>

   <h3>Input</h3>
   <p>
   Two or more variant sets to combine. They should be of non-overlapping genome intervals and with the same
   samples (sorted in the same order). If the files are ordered according to the appearance of intervals in the ref
   genome, then one can use the -assumeSorted flag.
   </p>

   <h3>Output</h3>
   <p>
   A combined VCF or BCF. The output file should have the same extension as the input(s).
   <\p>

   <h3>Important note</h3>
   <p>This is a command-line utility that bypasses the GATK engine. As a result, the command-line you must use to
   invoke it is a little different from other GATK tools (see example below), and it does not accept any of the
   classic "CommandLineGATK" arguments.</p>

   <h3>Usage example</h3>
   <pre>
   java -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
      -R reference.fasta \
      -V input1.vcf \
      -V input2.vcf \
      -out output.vcf \
      -assumeSorted
   </pre>

   <h3>Caveat</h3>
   <p>Currently the tool is more efficient when working with VCFs than with BCFs.</p>
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
- doc: genome reference file <name>.fasta
  id: reference
  type: File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--reference", inputs['reference_tags']))
  secondaryFiles:
  - .fai
  - ^.dict
- doc: genome reference file <name>.fasta
  id: reference
  type: File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--reference", inputs['reference_tags']))
  secondaryFiles:
  - .fai
  - ^.dict
- doc: assumeSorted should be true if the input files are already sorted (based on
    the position of the variants)
  id: assumeSorted
  type: boolean?
  inputBinding:
    prefix: --assumeSorted
- doc: assumeSorted should be true if the input files are already sorted (based on
    the position of the variants)
  id: assumeSorted
  type: boolean?
  inputBinding:
    prefix: --assumeSorted
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
- doc: output file
  id: outputFile
  type: File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--outputFile", inputs['outputFile_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'outputFile'
  id: outputFile_tags
- doc: output file
  id: outputFile
  type: File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--outputFile", inputs['outputFile_tags']))
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'outputFile'
  id: outputFile_tags
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'reference'
  id: reference_tags
- type:
  - 'null'
  - string
  - string[]
  doc: A argument to set the tags of 'reference'
  id: reference_tags
- doc: Input VCF file/s
  id: variant
  type:
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--variant", inputs['variant_tags']))
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'variant'
  id: variant_tags
- doc: Input VCF file/s
  id: variant
  type:
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--variant", inputs['variant_tags']))
- type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
  doc: A argument to set the tags of 'variant'
  id: variant_tags
- doc: the parameter (bin width or features per bin) to pass to the VCF/BCF IndexCreator
  id: variant_index_parameter
  type: int?
  inputBinding:
    prefix: --variant_index_parameter
- doc: which type of IndexCreator to use for VCF/BCF indices
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
