######################################
# wrapper for rule: raw_fastq_qc
######################################
import subprocess
from os.path import dirname
from os.path import basename
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, "wt")
f.write("\n##\n## RULE: raw_fastq_qc\n")
f.close()

version = str(subprocess.Popen("conda list ", shell = True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, "at")
f.write("\n##\n## CONDA: "+version+"\n")
f.close()

command = "mkdir -p " + dirname(snakemake.output.html) + " >> " + log_filename + " 2>&1"
f = open(log_filename, "at")
f.write("\n##\n## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "fastqc -o " + dirname(snakemake.output.html) + " --noextract --format fastq --nogroup --threads " \
    + str(snakemake.threads) + " " + snakemake.input.raw_fastq + " >> " + log_filename + " 2>&1"
f = open(log_filename, "at")
f.write("\n##\n## COMMAND: "+command+"\n")
f.close()
shell(command)