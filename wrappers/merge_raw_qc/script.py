######################################
# wrapper for rule: merge_raw_qc
######################################
import subprocess
from os.path import dirname
from os.path import basename
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: merge_raw_qc \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

if snakemake.params.long_reads:
	folder = "./qc_reports/*/raw_nanoplot/*"
else:
	folder =  "./qc_reports/*/raw_fastqc/*"

command = "multiqc -f -n " + snakemake.output.report + " " + folder + " --cl-config \"{{read_count_multiplier: 0.001, read_count_prefix: 'K', read_count_desc: 'thousands' }}\" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND:\n"+command+"\n")
f.close()
shell(command)
