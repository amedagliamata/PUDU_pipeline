######################################
# wrapper for rule: kraken2
######################################
import subprocess
import os
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)
f = open(log_filename, 'wt')
f.write("\n##\n## RULE: kraken2 \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

if snakemake.params.is_paired:
    flag = " --paired " + snakemake.input.paired[0] + " " + snakemake.input.paired[1]
else:
    flag = snakemake.input.paired[0]


command = "kraken2 --db " + os.path.dirname(snakemake.input.database) + " --threads " + str(snakemake.threads) + " --report " + snakemake.output.kraken_report \
    + " --output " + snakemake.output.classified + " --confidence " + str(snakemake.params.conf_interval) + " --minimum-base-quality " + str(snakemake.params.min_base_qual) \
    + " " + flag + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
