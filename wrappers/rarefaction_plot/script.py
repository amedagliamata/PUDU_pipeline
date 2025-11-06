######################################
# wrapper for rule: rarefaction_plot
######################################
import subprocess
import os
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)
f = open(log_filename, 'wt')
f.write("\n##\n## RULE: rarefaction_plot \n##\n")
f.close()

command = "python3 " +os.path.abspath(os.path.dirname(__file__))+ "/rarefaction_curve.py --input "+ snakemake.input.kraken_report + " --taxlevel " + snakemake.params.tax_lvl + " --output " + snakemake.output.rare_plot + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("COMMAND:"+command+"\n")
f.close()
shell(command)
