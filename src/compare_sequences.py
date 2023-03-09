import numpy as np
import itertools
from read_write_fasta import *
from binary_search_tree import * 

def get_kmers(input_sequence, k):
    '''
    This function creates a dictionary contataining  
    keys: kmers length of k;
    values: corresponding position on the input_sequence
    '''
    num_kmers = len(input_sequence) - k + 1
    # loop over the kmer start positions
    kmer_list = list()
    for i in range(num_kmers):
        # slice the string to get the kmer
        kmer = input_sequence[i:i+k]
        kmer_list.append(kmer)

    kmer_wpos = dict()
    for kmer in np.unique(kmer_list):
        kmer_indices  = [index+1 for (index, item) in enumerate(kmer_list) if item == kmer]
        kmer_wpos[kmer] = list(kmer_indices)    # look for index of element matching the kmer
    return kmer_wpos

# function to compute number key for each kmer
def compute_keynum(kmer):
    '''
    This function computes numerical representation i.e. numeric kmer (keynum) of a given kmer
    '''
    kmer_char = list(reversed([*kmer]))
    kmer_len = len(kmer)
    key_val = 0
    char_val = 0
    for i in reversed(range(kmer_len)):
        char = kmer_char[i]
        if char == "A":
            char_val = 0
        elif char == "C":
            char_val = 1
        elif char == "G":
            char_val = 2
        elif char == "T":
            char_val = 3
        key_val = key_val + char_val*(4**i) 
    return key_val

def format_database(reads_dict = read_fastaReads(), k = 3):
    '''
    This function computes key for each kmer in the database and store in a dictionary where
    keys: key numbers i.e. numeric kmers; 
    values: kmers
    '''
    db_seqs_list = list(reads_dict.values())
    kmer_keynum = dict()
    # get unque k-mer list
    unique_kmer_list = list()
    for db_seq in db_seqs_list:
        kmer_wpos = get_kmers(input_sequence = db_seq, k = k) # index to get kmer_list
        kmer_list = list(kmer_wpos.keys())
        unique_kmer_list = set(list(unique_kmer_list) + kmer_list)
    for kmer in unique_kmer_list:
        keynum = compute_keynum(kmer)
        kmer_keynum[keynum] = kmer
    return kmer_keynum 

def compute_kmer_within_HSSP1(kmer): 
    '''
    This function modifies one position in the input kmer by replacing one position
    at a time by other bases to get a list of modified kmers with score falling within HSSP threshold
    which would be later used for match search 
    '''
    # needs to be rewritten so it's universal for any k
    ## scoring scheme: match = 1; mismatch = -1, gap insertion = -1
    # modify one position in kmer
    HSSP_kmers = [kmer]
    for pos in range(len(kmer)):
        base_list = ["A","T","C","G"]
        # remaining bases
        base_list.remove(kmer[pos])
        # modify the pos position with the remaining bases
        new_kmer0 = kmer.replace(kmer[pos], base_list[0])
        new_kmer1 = kmer.replace(kmer[pos], base_list[1])
        new_kmer2 = kmer.replace(kmer[pos], base_list[2])
        HSSP_kmers.extend([new_kmer0, new_kmer1, new_kmer2])
    return HSSP_kmers 

def build_tree(elements, pos_list):
    '''
    This function build a binary search tree where 
    nodes: elements and their corresponding positions on the sequence being searched on
    '''
    root = BinarySearchTreeNode(elements[0], pos_list[0])
    for i in range(1,len(elements)):
        root.add_node(elements[i], pos_list[i])
    return root

def match_search(query_sequence, db_sequence, k = 3):
    '''
    This function builds a binary search tree from db_sequence
    and perform match searches on kmers from query_sequence to create a seeding dictionary where
    keys: seed positions found in db_sequence;
    values: key numbers (numeric kmers)
    '''
    seeding_dict = dict()
    search_kmers = list(get_kmers(query_sequence, k).keys())
    search_list = [compute_keynum(kmer) for kmer in search_kmers]
    kmer_wpos =  get_kmers(db_sequence, k)
    kmer_list = list(kmer_wpos.keys())
    keys = [compute_keynum(kmer) for kmer in kmer_list]
    pos_list = list(kmer_wpos.values())
    key_tree = build_tree(keys, pos_list)
    for search_index in range(len(search_list)):
        search_key = search_list[search_index]
        found_pos = key_tree.search(search_key)

        if found_pos:
            # record the found search key, kmer position and also its position on the db_sequence
            for y_pos in found_pos:
                x_pos = search_index
                seeding_pos = tuple([x_pos, y_pos])
                seeding_dict[seeding_pos] = search_key
    if bool(seeding_dict) == False:
        seeding_dict = "No seeeds found"
    return seeding_dict

def smith_waterman(query_sequence, db_sequence, 
                   match_score = 1, 
                   mismatch_score = -1, 
                   gap_penalty = -1):
    '''
    This function performs local alignment between a database and query sequences
    To be completed: extends seedings to perform local alignment between a database and query sequences
    '''
    align_mat = np.zeros((len(query_sequence) + 1, len(db_sequence) + 1), np.int)
    traceback_mat = np.zeros((len(query_sequence) + 1, len(db_sequence) + 1), np.int)
    max_score = 0
    for i, j in itertools.product(range(1, align_mat.shape[0]), range(1, align_mat.shape[1])):
        left_to_right = align_mat[i, j - 1] + gap_penalty
        top_to_bottom = align_mat[i - 1, j] + gap_penalty
        diagonal = align_mat[i - 1, j - 1] + (match_score if query_sequence[i - 1] == db_sequence[j - 1] else + mismatch_score)
        align_mat[i, j] = max(diagonal, top_to_bottom, left_to_right, 0)
        
        if align_mat[i,j] == 0: traceback_mat[i,j] = 0
        if align_mat[i,j] == left_to_right: 
            traceback_mat[i,j] = 1
        if align_mat[i,j] == top_to_bottom: 
            traceback_mat[i,j] = 2
        if align_mat[i,j] == diagonal: 
            traceback_mat[i,j] = 3
        # if we have a short list of kmers falling within HSPP score 
        # then we can get the index and iterate through the list to extend the seeds
        if align_mat[i][j] >= max_score:
            max_i = i
            max_j = j
            max_score = align_mat[i][j]

    align1, align2 = '', ''

    i,j = max_i,max_j # this i, j would correspond to each seed!!!
    # Traceback
    # Step 3. Traceback. Starting at the highest score in the scoring matrix 
    # H and ending at a matrix cell that has a score of 0, traceback based on the source of 
    # each score recursively to generate the best local alignment.
    while traceback_mat[i][j] != 0:
        if traceback_mat[i][j] == 3:
            base1 = query_sequence[i-1]
            base2 = db_sequence[j-1]
            i -= 1
            j -= 1
        elif traceback_mat[i][j] == 2:
            base1 = '-'
            base2 = db_sequence[j-1]
            j -= 1
        elif traceback_mat[i][j] == 1:
            base1 = query_sequence[i-1]
            base2 = '-'
            i -= 1
        align1 += base1
        align2 += base2

    align1 = align1[::-1]
    align2 = align2[::-1]
    symbol = ''
    identity_score = 0
        
    for i in range(len(align1)):
        base1 = align1[i]
        base2 = align2[i]
        if base1 == base2:                
            symbol += base1
            identity_score += 1
        elif base1 != base2 and base1 != '-' and base2 != '-': 
            symbol += ' '
        elif base1 == '-' or base2 == '-':          
            symbol += ' '
    percent_identity = identity_score / len(query_sequence) * 100
    return percent_identity, max_score, align1, align2, symbol
    

