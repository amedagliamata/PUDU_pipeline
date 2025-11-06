######################################
# wrapper for rule: emu_analysis
######################################
import subprocess
import os
from snakemake.shell import shell
shell.executable("/bin/bash")


tax_lvl = snakemake.params.tax_lvl

tax_map = {
    "S": "species",
    "G": "genus",
    "F": "family",
    "O": "order",
    "C": "class",
    "P": "phylum",
    "D": "domain"
}

new_tax_lvl = tax_map.get(tax_lvl)  # 

log_filename = str(snakemake.log)
f = open(log_filename, 'wt')
f.write("\n##\n## RULE: bracken_analysis \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "emu abundance " + snakemake.input.processed +  " --keep-counts  --db " + os.path.dirname(snakemake.input.db) + " --threads " + str(snakemake.threads)  + " --output-dir " + os.path.dirname(snakemake.output.result) + " --output-basename " + os.path.basename(snakemake.input.processed).replace("_R1.fastq.gz","") + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "emu collapse-taxonomy " + snakemake.output.result + " " + new_tax_lvl +  " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
