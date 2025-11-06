#############################################################
# wrapper for rule: emu_to_phyloseq
#############################################################
import os
import time
import resource
import platform
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'wt')
f.write("\n##\n## RULE: emu_to_phyloseq \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/emu_to_phyloseq.R "+\
            str(snakemake.input.experiment_design)+ " " +\
            str(snakemake.output.rds)+ " " +\
            " ".join(snakemake.input.emu_files) + " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)