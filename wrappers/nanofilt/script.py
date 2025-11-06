######################################
# wrapper for rule: nanofilt
######################################
import subprocess
from os.path import dirname
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: merge \n##\n")
f.close()

command = "mkdir -p " + dirname(snakemake.output.processed) + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "zcat " + str(snakemake.input.fastq) + " | NanoFilt --quality " + str(snakemake.params.quality) \
    + " --length " + str(snakemake.params.length) + " --maxlength " + str(snakemake.params.maxlength) \
    + " --minGC " + str(snakemake.params.mingc) + " --maxGC " + str(snakemake.params.maxgc) \
    + " --headcrop " + str(snakemake.params.headcrop) + " --tailcrop " + str(snakemake.params.tailcrop) + " > " + str(snakemake.output.processed).replace(".gz", "") + " 2>> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "gzip " + str(snakemake.output.processed).replace(".gz", "") + " 2>> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND:"+command+"\n")
f.close()
shell(command)