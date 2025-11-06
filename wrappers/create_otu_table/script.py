#######################################
# wrapper for rule: otu_table_creation
#######################################

import subprocess
import os
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)
f = open(log_filename, 'wt')
f.write("\n##\n## RULE: otu_table_creation \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "python3 " + os.path.abspath(os.path.dirname(__file__)) + "/kraken2otu.py --inputfolder " + str(os.path.dirname(snakemake.input.kraken_reports[0])) + " --level " + snakemake.params.lv \
    + " --extension _kraken2_report.txt >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
