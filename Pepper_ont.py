import subprocess
import os
import argparse
import multiprocessing
from pathlib import Path


# **Note run these commands in terminal: $ sudo groupadd docker, sudo usermod -aG docker $USER,
# newgrp docker to run the docker command. Create a separate conda environment to run PEPPER
class ont_reads:
    def __init__(self, args):
        self.query_folder = args.query_folder
        self.reference_folder = args.reference_folder
        query_file_ls = sorted(os.listdir(self.query_folder))
        reference_folder = self.reference_folder
        query_folder = self.query_folder
        Methods_Calling.out_put(query_file_ls, reference_folder, query_folder)


class Methods_Calling:
    @staticmethod
    def folder_create(output_folder):
        # Create output folder
        Path(output_folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def bam_file(ref_fold, base_name, bam, q_fold):
        # Generates the necessary bam files by piping samtools and minimap2
        # Location of the query fastq file
        query_location = str(q_fold) + '/' + base_name

        # Creating the bam file using minimap2 and samtools
        # ref_fold (reference folder), bam(bam file) and q_fold is the query folder
        minimap2_cmd = ['minimap2', '-x',
                        'map-ont', '-a', ref_fold, query_location]
        samtools_cmd = ['samtools', 'sort',
                        '-o', bam, '-@', '2']
        minimap2 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE)
        samtools = subprocess.Popen(samtools_cmd, stdin=minimap2.stdout, stdout=subprocess.PIPE)
        minimap2.stdout.close()
        samtools.communicate()

    @staticmethod
    def bai_file(bam, bai, ref_fold):
        # Indexes the bam file using samtools and moves files to the ref folder
        bai_cmd = ['samtools', 'index', bam, bai]
        mv_file = ['mv', bam, bai, ref_fold]
        subprocess.run(bai_cmd)

    @staticmethod
    def gvcf(output_folder, output_prefix, cpu, bam, ref_fold, ref_name):
        threads = str(cpu)
        # Pull docker
        Pull_cmd = ['docker', 'pull', 'kishwars/pepper_deepvariant:r0.7']
        # Generate the gvcf and vcf files
        gvcf_cmd = ['docker', 'run',
                    '-v', '{}:/ref'.format(ref_fold),
                    '-v', '{}:/output'.format(output_folder),
                    'kishwars/pepper_deepvariant:r0.7',
                    'run_pepper_margin_deepvariant', 'call_variant',
                    '-b', "ref/{}".format(os.path.basename(bam)),
                    '-f', "ref/{}".format(os.path.basename(ref_name)),
                    '-o', "/output",
                    '-p', output_prefix, '-t', threads,
                    '--ont_r9_guppy5_sup', '--gvcf']
        # subprocess.run(Pull_cmd)
        subprocess.run(gvcf_cmd)

    # Report and File generator
    @staticmethod
    def Generate(q_list, ref_fold, q_fold):
        base_name = (os.path.basename(q_list).split(".")[0]).split("_")[0]
        # q_list is the list of all the query fastq files
        # q_fold in the directory of the query folder
        # ref_fold is the directory for the folder that contains the reference genome

        # Set up input data
        # The reference file and query files were downloaded from ncbi either
        # or using the sratoolkit
        directory = os.getcwd()

        # Name of the query file in the query list
        query_file = base_name + ".fastq"

        # Name of the reference file in the reference folder, modify this to you liking
        ref_name = directory + "/ref/ref.fasta"

        # Name of the bam file
        bam = directory + "/ref/" + base_name + "_pepper.bam"

        # Name of bam.bai file
        bai = directory + "/ref/" + base_name + "_pepper.bam.bai"

        # Set the number of CPU's
        p = 1
        cpu = 2

        # Set the output prefix
        output_prefix = base_name + '_PEPPER_Margin_DeepVariant'
        # Set up output folder
        output_folder = directory + "/" + "output"

        # Creates the output folder
        Methods_Calling.folder_create(output_folder)
        
        # Creates the bam file for each query file in the list
        Methods_Calling.bam_file(ref_fold, query_file, bam, q_fold)
       
        # Indexes each bam file that was created in the previous step
        Methods_Calling.bai_file(bam, bai, ref_fold)

        # Produces the GVCF and VCF files using variant calling
        Methods_Calling.gvcf(output_folder, output_prefix, cpu, bam, ref_fold, ref_name)

    @staticmethod
    def out_put(query_file_ls, reference_folder, query_folder):
        for i in range(len(query_file_ls)):
            p = multiprocessing.Process(target=Methods_Calling.Generate, args=(query_file_ls[i],
                                        reference_folder, query_folder))
            p.start()


if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="VCF and GVCF file generation")
    parse.add_argument("query_folder", help="folder contains fastq mycobacterium query files")
    parse.add_argument("reference_folder", help="folder contains fasta mycobacterium reference file")
    arg = parse.parse_args()
    ont_reads(arg)
