######################################
# wrapper for rule: kraken2_reference
######################################
import subprocess
import os
import resource
import platform
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)
f = open(log_filename, 'wt')
f.write("\n##\n## RULE: kraken2_reference \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()


command = "kraken2-build --special " + os.path.dirname(snakemake.output.database[0]).replace("DB", "") + " -db " + os.path.dirname(snakemake.output.database[0]) + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "kraken2-build --build -db " + os.path.dirname(snakemake.output.database[0]) #+ " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

command = "bracken-build -d " + os.path.dirname(snakemake.output.database[0]) + " -t " + str(snakemake.threads) + " -k " + str(snakemake.params.kmer_length) + " -l " + str(snakemake.params.read_length) # + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "touch " + snakemake.output.reference + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("# COMMAND: " + command + "\n")
f.close()
shell(command)
