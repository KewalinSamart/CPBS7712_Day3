'''
Perform alignment step using BLAST algorithm to filter potential 
candidate sequences in the database for assemply step (BLAST preprocessing)
'''

import argparse, time
from read_write_fasta import *
from compare_sequences import *
from binary_search_tree import *

parser = argparse.ArgumentParser(description='Reconstruct sequence from reads and query for substring PartI: BLAST preprocessing and exact local alignment')
parser.add_argument('--QUERYfasta', metavar='QUERYfasta', type=str, default='./data/toy_query.fasta', help='input file containing query sequene')
parser.add_argument('--READSfasta', metavar='READSfasta', type=str, default='./data/toy_READS.fasta', help='input file containing database sequences i.e. READS')
parser.add_argument('--kmer', metavar='kmer', type=int, default=3, help='length of node string')
parser.add_argument('--outputFile', metavar='outputFile', type=str, default='./outputs/stats_output.txt', help='name of file to write intermediate' 
                    'statistics output to')
args = parser.parse_args()

def main(QUERYfasta, READSfasta, kmer, outputFile):
    # read in query sequence
    query_sequence = read_fastaQuery(QUERYfasta)
    # read in reads and format sequence database
    db_sequences = read_fastaReads(READSfasta)

    # compute numerical representation (key number unique to each unique kmer)
    kmer_keynum = format_database(db_sequences,k=kmer)
 
    output_file = open(outputFile,"w")
    output_file.write("Query sequence: " + query_sequence + "\n")
    output_file.write("Exploratory alignment results: " + "\n")
  
    for db_key, db_seq in db_sequences.items():
        output_file.write("-----------------------------" + "\n")
        output_file.write(db_key + ": " + db_seq + "\n")
        output_file.write("Length of database sequence: " + str(len(db_seq)) + "\n")
        
        start_seedsearch = time.time()
        seeding_dict = match_search(query_sequence, db_seq)
        stop_seedsearch = time.time()
        output_file.write("Time spent searching for seeds: " + str(stop_seedsearch-start_seedsearch) + " seconds" + "\n")
       
        if isinstance(seeding_dict,str): 
            output_file.write("No seeds found"  + "\n")
        else:
            for position, keynum in seeding_dict.items():
                output_file.write("kmer and numeric key: " + kmer_keynum[keynum] + "," + str(keynum)  + "\n")
                output_file.write("Seed position: " + str(position[0]) + "," + str(position[1])  + "\n")
        
        # explore data with Smith-Waterman local alignment
        start_waterman = time.time()
        percent_identity, max_score, align1, align2, symbol = smith_waterman(query_sequence, db_seq)
        stop_waterman = time.time()
        output_file.write("Time Smith-Waterman: " + str(stop_waterman-start_waterman) + " seconds" + "\n")
       
        percent_identity = round(percent_identity,2)
        output_file.write("percent matched: " + str(percent_identity) + "\n")
        output_file.write("Smith-Waterman max score: " + str(max_score) + "\n")
        output_file.write(align1 + "\n")
        output_file.write(symbol + "\n")
        output_file.write(align2 + "\n")

    output_file.close()

    print(outputFile, "was successfully written")

main(args.QUERYfasta, args.READSfasta, args.kmer, args.outputFile)