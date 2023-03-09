from read_write_fasta import *

full_seq = read_fastaQuery()
short_seq = full_seq[20:30]
output_file = open("/Users/kewalinsamart/Documents/GitHub/CPBS7712_ScratchRepo/toy_inputs/" + "toy_query.fasta","w")
output_file.write(">" + "queryseq" + "\n" + short_seq + "\n")
output_file.close()

output_file = open("/Users/kewalinsamart/Documents/GitHub/CPBS7712_ScratchRepo/toy_inputs/" + "toy_READS.fasta","w")
output_file.write(">" + "dbseq1" + "\n" + full_seq[:22] + "\n")
output_file.write(">" + "dbseq2" + "\n" + full_seq[20:35] + "\n")
output_file.write(">" + "dbseq3" + "\n" + full_seq[30:50] + "\n")
output_file.write(">" + "dbseq4" + "\n" + full_seq[45:58] + "\n")
output_file.write(">" + "dbseq5" + "\n" + full_seq[50:70] + "\n")
#output_file.close()