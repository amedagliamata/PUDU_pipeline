#########################################
# wrapper for rule: centrifuger_nanopore
#########################################

import subprocess
from os.path import dirname
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)
f = open(log_filename, 'wt')
f.write("\n##\n## RULE: centrifuger_nanopore \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p " + dirname(snakemake.output.classified) + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if snakemake.params.is_paired:
    flag = " -1 " + snakemake.input.processed[0] + " -2 " + snakemake.input.processed[1]
else:
    flag = " -u " + snakemake.input.processed[0]

command = "centrifuger -x " + snakemake.params.centrifuge_db \
    + flag \
    + " --min-hitlen " + str(snakemake.params.min_hitlen) + " -t 10 -k " + str(snakemake.params.max_hits) \
    + " > " + snakemake.output.classified
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "centrifuger-kreport -x " + snakemake.params.centrifuge_db \
    + " " + snakemake.output.classified + " > " + snakemake.output.kreport #+ " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
