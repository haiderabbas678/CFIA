import subprocess


a="ERR4627396.sra_1.fastq"
b="ERR4627396.sra_2.fastq"
d="ERR4627396.skesa.fa"


# fastp_cmd = ["fastp","-i", "{}".format(a),"-I", "{}".format(b),"-o", "ERR4627396.sra_1.fq", "-O", "ERR4627396.sra_2.fq"] 

# fastp = subprocess.run(fastp_cmd, stdout=subprocess.PIPE) 
c="ERR4627396.sra_1.fq"
e="ERR4627396.sra_2.fq"

skesa_cmd = ["skesa", "--fastq", "{},{}".format(c,e), ">", "{}".format(d)]
print(skesa_cmd)
# skesa =  subprocess.run(skesa_cmd, stdout=subprocess.PIPE)