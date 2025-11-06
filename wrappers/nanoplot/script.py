######################################
# wrapper for rule: nanoplot
######################################
import subprocess
from os.path import dirname
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: nanoplot \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "NanoPlot --fastq " + snakemake.input.fastq + " --no_static --outdir " + snakemake.output.out + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + snakemake.output.out + "/NanoStats.txt " + snakemake.output.out + "/" + snakemake.params.sample + "_NanoStats.txt >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)