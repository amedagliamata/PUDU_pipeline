######################################
# wrapper for rule: raw_nanoplot_qc
######################################
import subprocess
from os.path import dirname
from os.path import basename
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: raw_nanoplot_qc \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "NanoPlot --fastq " + snakemake.input.fastq + " -o " + snakemake.output.out + " >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND:\n"+command+"\n")
f.close()
shell(command)