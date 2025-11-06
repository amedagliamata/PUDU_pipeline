######################################
# wrapper for rule: krona_plots
######################################

import subprocess
import os
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)
f = open(log_filename, 'wt')
f.write("\n##\n## RULE: krona_plots \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "ktImportTaxonomy -t 2 -m 5 " + snakemake.input.kraken_report + " -o " + str(snakemake.output.krona_graph) + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
