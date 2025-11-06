######################################
# wrapper for rule: trimmomatic
######################################
import subprocess
import warnings
from os.path import dirname
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: trimmomatic \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

if snakemake.params.trim_adapters and not snakemake.params.qc_trim_reads:
    warnings.warn("Performing trimming qc as adapter removal was set to True")
    snakemake.params.qc_trim_reads = True

is_paired = len(snakemake.input.fastq) == 2

flag = ""
if snakemake.params.qc_trim_reads:
    if snakemake.params.trim_adapters:
        flag = " ILLUMINACLIP:adapters/" + snakemake.params.adapters + ".fa:2:30:10"
    else:
        flag = ""

if is_paired:
    trimmomatic_call = "trimmomatic PE " + snakemake.input.fastq[0] + " " + snakemake.input.fastq[1] + " " + snakemake.output.paired[0] + " " + str(snakemake.output.paired[0]).replace(".fastq", "_unpaired_R1.fastq") + " " + snakemake.output.paired[1] + " " + str(snakemake.output.paired[1]).replace(".fastq", "_unpaired_R2.fastq")
else:
    trimmomatic_call = "trimmomatic SE " + snakemake.input.fastq[0] + " " + snakemake.output.paired[0]

if snakemake.params.trim_adapters or snakemake.params.qc_trim_reads:

    command = trimmomatic_call + " " + flag + " LEADING:" + str(snakemake.params.leading) + " TRAILING:" + str(snakemake.params.trailing) \
        + " SLIDINGWINDOW:" + str(snakemake.params.slidingwindow) + " MINLEN:" + str(snakemake.params.minlen) \
        + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

else:

    command = "mkdir -p " + dirname(snakemake.output.paired[0]) + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = "cp " + snakemake.input.fastq[0] + " " + snakemake.output.paired[0] + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    if is_paired:

        command = "cp " + snakemake.input.fastq[1] + " " + snakemake.output.paired[1] + " >> " + log_filename + " 2>&1"
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)