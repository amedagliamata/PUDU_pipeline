######################################
# wrapper for rule: bracken_analysis
######################################
import subprocess
import os
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)
f = open(log_filename, 'wt')
f.write("\n##\n## RULE: bracken_analysis \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "bracken -d " + os.path.dirname(snakemake.input.database) + " -i " + snakemake.input.kraken_report + " -o " + snakemake.output.ab_report + " -l G -r " + str(snakemake.params.read_length) #+ "  >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
