from os import sys, path
# Use fix from https://stackoverflow.com/a/19190695 to import from the base directory
sys.path.append(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))

import gatkcwlgenerator as cwl_gen

import subprocess
import os
from os import path
from multiprocessing import Process
import tempfile
import pytest

import requests

"""
Runs the specified command and reports it as an AssertionError if it fails (can override this with
expect_failure)
"""
def run_command(command, fail_message=None, expect_failure=False):
    process = subprocess.Popen(command.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    exitcode = process.returncode

    if exitcode != 0 and not expect_failure:
        raise AssertionError("{}\"{}\" fails with exit code {}\nstdout:\n{}\nstderr:\n{}".format(
            "" if fail_message is None else fail_message + ": ",
            command,
            exitcode,
            stdout,
            stderr
        ))

    return CommandOutput(stdout, stderr, exitcode)

class CommandOutput():
    def __init__(self, stdout, stderr, exitcode):
        self.stdout = stdout
        self.stderr = stderr
        self.exitcode = exitcode


base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
os.chdir(base_dir) # Need to be in the base directory for the cwl-runner to pick up the correct files

default_args = """
reference_sequence:
   class: File
   #path: /path/to/fasta/ref/file
   path: {0}/cwl-example-data/chr22_cwl_test.fa 
refIndex:
   class: File
   #path: /path/to/index/file
   path: {0}/cwl-example-data/chr22_cwl_test.fa.fai 
refDict:
   class: File
   #path: /path/to/dict/file
   path: {0}/cwl-example-data/chr22_cwl_test.fa.dict 
input_file: #must be BAM or CRAM
   class: File
   #path: /path/to/input/file 
   path: {0}/cwl-example-data/chr22_cwl_test.cram
out: out.gvcf.gz""".format(base_dir)

def run_haplotype_caller(extra_info, version, interval=1, filetext=None, expect_failure=False):
    return run_tool("HaplotypeCaller", extra_info, version, interval, filetext, expect_failure)

"""
Runs the haplotype_caller tool with the specified data
"""
def run_tool(toolname, extra_info, version, interval=1, filetext=None, expect_failure=False):
    if filetext is None:
        extra_info += "\nintervals: [chr22:10591400-{}]".format(10591400 + interval)
        filetext = "analysis_type: {}\n".format(toolname) + default_args + "\n" + extra_info

    with tempfile.NamedTemporaryFile() as f:
        f.write(filetext)
        f.flush()
        return run_command("cwl-runner gatk_cmdline_tools/{}/cwl/{}.cwl {}".format(
            version,
            toolname,
            f.name
        ), expect_failure=expect_failure)

# Unit tests

supported_versions = ["3.5", "3.6", "3.7", "3.8"]

@pytest.mark.parametrize("version", supported_versions)
class TestGenerateCWL:
    def test_get_json_links(self, version):
        json_links = cwl_gen.get_json_links(version)
        for link_type, links in json_links.__dict__.items():
            assert links, "There are no links of type '{}' in gatk version {}".format(
                    link_type,
                    version
                )

    def test_no_arguments_in_annotator(self, version):
        # If arguments are in annotator modules, we probably need to add them to the CWL file
        for url in cwl_gen.get_json_links(version).annotator_urls:
            ann_json = requests.get(url).json()
            assert not ann_json["arguments"]

# Integration tests

@pytest.mark.parametrize("version", supported_versions + ["4.beta-latest"])
def test_runs(version):
    run_command("python2 gatkcwlgenerator/main.py -v {} --dev".format(version))

@pytest.mark.parametrize("version", supported_versions)
class TestGeneratedCWLFiles:
    def get_base_cwl_path(self, version):
        return path.join(base_dir, "gatk_cmdline_tools/{}/cwl".format(version))
    
    def test_are_cwl_files_valid(self, version):
        exceptions = []
        for cwl_file in os.listdir(self.get_base_cwl_path(version)):
            try:
                print("Validated " + cwl_file)
                run_command("cwl-runner --validate " + path.join(self.get_base_cwl_path(version), cwl_file))
            except AssertionError as e:
                print(e)
                exceptions.append(e)

        if exceptions:
            raise AssertionError("Not all cwl files are valid:\n" + "\n\n".join(exceptions))

    def test_haplotype_caller(self, version):
        run_command("cwl-runner gatk_cmdline_tools/{}/cwl/HaplotypeCaller.cwl examples/HaplotypeCaller_inputs.yml".format(version))

    # Test if the haplotype caller accepts all the correct types

    def test_boolean_type(self, version):
        assert "ThreadEfficiencyMonitor" in run_haplotype_caller("monitorThreadEfficiency: True", version).stderr

    def test_integers_type(self, version):
        assert "42 data thread" in run_haplotype_caller("num_threads: 42", version, expect_failure=True).stderr

    def test_string_type(self, version):
        assert "Specified name does not exist in input bam files" in \
            run_haplotype_caller("sample_name: invalid_sample_name", version, expect_failure=True).stderr

    def test_file_type(self, version):
        BQSR_arg = """
BQSR:
   class: File
   path: {0}/cwl-example-data/chr22_cwl_test.fa
""".format(base_dir)
        assert "Bad input: The GATK report has an unknown/unsupported version in the header" in \
            run_haplotype_caller(BQSR_arg, version, expect_failure=True).stderr

    def test_enum_type(self, version):
        assert "Strictness is LENIENT" in run_haplotype_caller("validation_strictness: LENIENT", version).stderr

    def test_list_type(self, version):
        run_with_larger_intervals = run_haplotype_caller(extra_info="", version=version,
            filetext="analysis_type: HaplotypeCaller\n" + default_args + "\nintervals: [chr22:10591400-10591500, chr22:10591500-10591645]")

        assert "Processing 246 bp from intervals" in run_with_larger_intervals.stderr