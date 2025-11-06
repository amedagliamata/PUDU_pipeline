#############################################################
# wrapper for rule: dada2_analysis
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'wt')
f.write("\n##\n## RULE: dada2_analysis \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/dada2_computation.R "+\
            str(snakemake.input.experiment_design)+ " " +\
            str(snakemake.params.for_trunc)+ " " +\
            str(snakemake.params.rev_trunc)+ " " +\
            str(snakemake.params.multithread)+ " " +\
            str(snakemake.input.train_set)+ " " +\
            str(snakemake.input.sp_assign)+ " "+\
            str(snakemake.params.for_error)+ " "+\
            str(snakemake.params.rev_error)+ " "+\
            " ".join(snakemake.input.fastq) + " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)