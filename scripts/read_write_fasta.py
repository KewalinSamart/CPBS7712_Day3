import re
import sys

def read_fastaQuery(query_path = './data/toy_query.fasta'):
    # read in query fasta file 
    f_query = open(query_path,'r')
    lines_query = f_query.readlines()
    query_seq = lines_query[1]

    return query_seq

def read_fastaReads(reads_path = './data/toy_READS.fasta'):
     # read in READS fasta file (database sequences) 
    f_reads = open(reads_path,'r')
    lines_reads = f_reads.readlines()

    hre = re.compile('>(\S+)')
    lre = re.compile('^(\S+)$')
    reads_dict = {}

    for line in lines_reads:
        outh = hre.search(line)
        if outh:
                id=outh.group(1)
        else:
                outl=lre.search(line)
                if(id in reads_dict.keys()):
                        reads_dict[id] += outl.group(1)
                else:
                        reads_dict[id]  =outl.group(1)
    return reads_dict 

def write_contig(output_dir, contig):
    # write output contig sequence in fasta format
    if isinstance(output_dir, str):
        output_file = open(output_dir + "ALLELES.fasta","w")
        output_file.write(">" + "finalcontig" + "\n" + contig + "\n")
        output_file.close()
    else:
        sys.exit("output_dir must be type of string")
