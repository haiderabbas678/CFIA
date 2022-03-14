import subprocess
import os
import argparse
import multiprocessing


class short_reads:
    def __init__(self, args):
        self.fold = args.fold
        a = sorted(os.listdir(self.fold))
        a.remove("genome_assembly.py")
        b = MethodsCalling.Sorting(a)
        MethodsCalling.out_put(b)


class MethodsCalling:
    @staticmethod
    def Sorting(los):
        lis = []
        length = len(los) // 2
        for x in range(length):
            lis.append(los[:2])
            los = los[2:]
        return lis

    # Trimming Function:
    #    Responsible for removing adapter sequences, filtering by quality, and read pruning
    @staticmethod
    def fastp(r1, r2, t_r1, t_r2):
        fastp_cmd = ["fastp",
                     "--thread", "3",
                     "-i", "{}".format(r1),
                     "-I", "{}".format(r2),
                     "-o", t_r1,
                     "-O", t_r2]
        return subprocess.run(fastp_cmd)

    # Genome Assembly Function:
    #    Assembles the genome use the trimmed forward and reverse sequences
    @staticmethod
    def skesa(t_r1, t_r2, sk):
        skesa_cmd = ["skesa",
                     "--fastq", "{},{}".format(t_r1, t_r2),
                     "--contigs_out", sk,
                     "--cores", "2"]
        return subprocess.run(skesa_cmd)

    # Genome indexing Function
    #    Indexes the genome making it easier to find sequences of interest, mismatches, indels etc
    @staticmethod
    def bwa_index(sk):
        bwaindex_cmd = ["bwa", "index", sk]
        return subprocess.run(bwaindex_cmd)

    # Function for generating a sam file
    #    Maps the reads back to the reference genome to generate the sam file
    @staticmethod
    def bwamem(sk, r1, r2, sam):
        bwamem_cmd = ["bwa",
                      "mem",
                      "-t", "3",
                      "-M", sk, "{}".format(r1), "{}".format(r2),
                      "-o", sam]
        return subprocess.run(bwamem_cmd)

    # Function for converting a sam file into a bam file
    @staticmethod
    def samtool(sam, bam):
        samtools_cmd = ["samtools",
                        "sort", sam,
                        "-o", bam,
                        "-@", "3"]
        return subprocess.run(samtools_cmd, stdout=subprocess.PIPE)

    # Performs a qc check on the genome assembly
    @staticmethod
    def qualimap(bam, qm):
        qualimap_cmd = ["qualimap",
                        "bamqc",
                        "-bam", bam,
                        "-outfile", qm,
                        "-nt", "3"]
        return subprocess.run(qualimap_cmd)

    # Report and File generator
    @staticmethod
    def Generate(l):
        r1 = l[0]
        r2 = l[1]
        base_name = (os.path.basename(r1).split(".")[0]).split("_")[0]
        t1 = base_name + "_t1.fa"
        t2 = base_name + "_t2.fa"
        skesa = base_name + ".skesa.fa"
        sam = base_name + ".sam"
        bam = base_name + ".bam"
        qm = base_name + "_qualimap"
        a = MethodsCalling.fastp(r1, r2, t1, t2)
        b = MethodsCalling.skesa(t1, t2, skesa)
        c = MethodsCalling.bwa_index(skesa)
        d = MethodsCalling.bwamem(skesa, r1, r2, sam)
        e = MethodsCalling.samtool(sam, bam)
        f = MethodsCalling.qualimap(bam, qm)


    @staticmethod
    def out_put(l):
        for i in range(len(l)):
            p = multiprocessing.Process(target=MethodsCalling.Generate, args=(l[i],))
            p.start()


if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="Genome assembly and QC")
    parse.add_argument("fold", help="folder contains fastq files")
    arg = parse.parse_args()
    short_reads(arg)
