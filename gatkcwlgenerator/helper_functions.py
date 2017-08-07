#!/bin/python
"""
Collection of helper functions for cwl_generator.py and json2cwl.py
"""


def get_input_json(argument):
    """
    Returns the cwl syntax for expressing the given argument

    :param argument: The cwl argument, as specified in the json file
    :param inputs: The inputs object, to be written to with the correct information

    :returns: CWL object to describe the given argument
    """

    cwl_desc = {
        'doc': argument['summary'],
        'id': argument['name'].strip('-'),
    }

    prefix = argument['name']
    # Patch the incorrect description given by GATK for both --input_file
    if prefix == '--input_file' or prefix == '--input':
        argument['name'] = 'File'

    typ, is_array_type = get_CWL_type(argument)
    if not is_array_type:
        cwl_desc["inputBinding"] = {
            "prefix": argument["name"]
        }

    cwl_desc["type"] = typ

    if argument['defaultValue'] != "NA" and argument['defaultValue'] != "None":
        default_helper(cwl_desc, argument)

    if argument['name'] == '--reference_sequence':
        cwl_desc['secondaryFiles'] = ['.fai', '^.dict']
    elif 'requires' in argument['fulltext'] and 'files' in argument['fulltext']:
        cwl_desc['secondaryFiles'] = "$(self.location+'.'+self.basename.split('.').splice(-1)[0].replace('m','i'))"

    return cwl_desc


# You cannot get the enumeration information for an enumeration in a nested type, so they are hard coded here
enum_types = {
    # Example: https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_tools_walkers_variantutils_ValidateVariants.php
    "validationtype": ["ALL", "REF", "IDS", "ALLELES", "CHR_COUNTS"],
    # Example: https://software.broadinstitute.org/gatk/gatkdocs/3.7-0/org_broadinstitute_gatk_tools_walkers_cancer_contamination_ContEst.php#--lane_level_contamination
    "contaminationruntype": ['META', 'SAMPLE', 'READGROUP'],  # default is META
    # Example: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.php#--partitionType
    "partition": ["readgroup", "sample", "library", "platform", "center",
                  "sample_by_platform", "sample_by_center", "sample_by_platform_by_center"],
    "type": ['INDEL', 'SNP', 'MIXED', 'MNP', 'SYMBOLIC', 'NO_VARIATION']
}


def basic_GATK_type_to_CWL(argument, typ):
    """
    Gets the correct CWL type for an argument, given an argument's GATK type excluding List[...]

    typ may not be the same as the arguments type if you are overriding the type

    :param argument: The cwl argument, as specified in the json file
    :param typ: The GATK type given

    :returns: (cwl_type, is_array_type)
    """

    is_array_type = False

    if typ in ('long', 'double', 'int', 'string', 'float', 'boolean', 'bool'):
        cwl_type = typ
    elif typ == 'file':
        cwl_type = 'File'
    elif typ in ('byte', 'integer'):
        cwl_type = 'int'
    elif typ == 'set':
        # The caller function will turn this into an array of sets
        # Example of this -goodSM
        is_array_type = True
        cwl_type = 'string'
    elif argument['options']:
        # Check for enumerated types, and if they exist, ignore the specified type name
        cwl_type = {
            'type': 'enum',
            'symbols': [x['name'] for x in argument['options']]
        }
    elif typ in enum_types.keys():
        cwl_type = {
            "type": "enum",
            "symbols": enum_types[typ]
        }
    elif 'intervalbinding' in typ:
        # Type intervalbinding can be a list of strings or a list of files
        is_array_type = True
        cwl_type = 'File'
    elif "rodbinding" in typ:
         # Possible options: https://gist.github.com/ThomasHickman/b4a0552231f4963520927812ad29eac8
        cwl_type = 'File'
    else:
        print('\033[93m'+ 'WARNING: Unable to assign to a CWL type, defaulting to string\n Argument: {}   Type: {}'.format(
            argument['name'][2:], typ))
        
        cwl_type = "string"

    return (cwl_type, is_array_type)


def typcash(argument, typ, defVal):
    if typ == 'int':
        return int(defVal)
    elif typ == 'boolean':
        return bool(defVal)
    elif typ == 'string':
        return defVal
    elif typ == 'enum':
        return defVal
    elif typ == 'long':
        return long(defVal)
    elif typ == 'double':
        return float(defVal)
    else:
        raise ValueError('Failed to cash default value to assigned type. \n Argument: {}    Type: {}   DefaultValue: {}'.format(
            argument['name'], typ, defVal))


def get_CWL_type(argument):
    """
    Returns the CWL type of an argument, indicating if it is an array type

    :param argument: The cwl argument, as specified in the json file
    :param cwl_desc: The inputs object, to be written to with the correct information

    :returns: (type_json, is_array_type)
    If is_array_type is true, the caller shouldn't output inputBinding's prefix
    """

    is_array_type = False

    prefix = argument['name']
    # Patch the incorrect description given by GATK for both --input_file and the type intervalbinding
    if prefix == '--input_file' or prefix == '--input':
        typ = 'File'
    else:
        outer_type = argument['type'].lower()

        if 'list[' in outer_type or 'set[' in outer_type:
            inner_type = outer_type[outer_type.index('[') + 1 : -1]
            is_array_type = True
        elif '[]' in outer_type:
            inner_type = outer_type.strip('[]')
            is_array_type = True
        else:
            inner_type = outer_type

        typ, is_array_type2 = basic_GATK_type_to_CWL(argument, inner_type)

        if is_array_type2:
            is_array_type = True

        if is_array_type:
            # This needs to be done instead of adding [], as you can do correct
            # inputBinding prefixes and it works with object types
            typ = {
                "type": "array",
                "items": typ,
                "inputBinding": {
                    "prefix": prefix
                }
            }

        if argument['required'] == 'no':
            if isinstance(typ, list):
                typ.insert(0, 'null')
            else:
                typ = ['null', typ]

    return (typ, is_array_type)


def default_helper(inpt, argument):
    """
    Sets the default value as in its original type
    """
    if False:#cmd_line_args.dont_generate_default:
        return

    typ = inpt['type']
    defVal = argument['defaultValue']

    try:
        if isinstance(typ, list):
            typ = typ[1]
        if isinstance(typ, dict):
            if typ['type'] == 'array':
                item_type = typ['items']
                typ = typ['type']
            else:  # if type == 'enum'
                typ = typ['type']

    except:
        raise ValueError('Failed to identify the type of the input argument \n Input argument: {}    Input type: {}'.format(
            argument['name'], typ))

    if defVal == '[]':
        inpt['default'] = []
    else:
        if '[]' in typ and typ != '[]':
            typ = typ.strip('[]')
            inpt['default'] = [typcash(argument, typ, val)
                               for val in defVal[1:-1].replace(' ', '').split(',')]
        elif typ == "array":
            inpt['default'] = [typcash(argument, item_type, val)
                               for val in defVal[1:-1].replace(' ', '').split(',')]
        else:
            inpt['default'] = typcash(argument, typ, defVal)


def is_output_argument(argument):
    """
    Returns whether this argument's type indicates it's an output argument
    """
    return any(x in argument["type"] for x in ('PrintStream', 'Writer'))


def get_output_json(argument):
    """
    Modifies the `outputs` parameter with the cwl syntax for expressing a given output argument

    The default for output files is guessed based on the file type e.g. GATKSamFileWriter generates <NAME>.bam

    :param argument Object: The cwl argument, as specified in the json file
    :param outputs Object: The outputs object, to be written to with the correct information

    :returns: (input_json, output_json)
    """

    def helper(argument, globval, typ):
        return {
            'id': argument['name'],
            'type': typ,
            'outputBinding': {
                'glob': globval
            }
        }

    prefix = argument['name']
    name = prefix.strip('-')

    if argument['type'] == "GATKSAMFileWriter":
        output_path = '{}.bam'.format(name)
    elif argument['type'] == "PrintStream":
        output_path = '{}.txt'.format(name)
    elif argument['type'] == 'VariantContextWriter':
        output_path = '{}.vcf'.format(name)
    else:
        output_path = '{}.txt'.format(name)
        print("The input argument '{}' with input type '{}' will create an output file to the following output path: {} ".format(
            argument['name'], argument['type'], output_path))
        # when not an required argument,
    argument['type'] = 'string'  # option as an input to specify filepath

    if argument['defaultValue'] == "NA":
        """
        if the argument is not required and doesn't have a default value
        it is an optional input to which file is generated to
        if specified, the outputbinding should be the specified value
        """
        input_json = get_input_json(argument)                                                # input
        output_json = helper(argument, '$(inputs.{})'.format(
            name), ['null', 'File'])  # optional outputbinding
    else:
        """
        if the argument is not required but has a default value
        reset the default value to sth that isn't standard output
        it is an optional input to which file is generated to
        if specified, the outputbinding should be the specified value
        else default
        """
        argument['defaultValue'] = output_path                                       # reset default
        # input
        input_json = get_input_json(argument)
        output_json = helper(argument, '$(inputs.{})'.format(
            name), 'File')          # always an output

    return (input_json, output_json)