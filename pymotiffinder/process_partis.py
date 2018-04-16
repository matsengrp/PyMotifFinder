import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# this should eventually include insertions as well as mutations, for
# now it just compares the naive sequence to the indel reversed
# sequences
def process_partis(partis_file, max_mutation_rate=1, use_indel_seqs=True):
    """Reads in mutations from a partis file

    Keyword arguments:
    partis_file -- A string giving the name of a partis file.

    Returns: A data frame with one row per mutation. Each row contains
    the mutated sequence, naive sequence, sequence id, the location of
    the mutation, the germline base, and the mutated base.

    """
    annotations = pd.read_csv(partis_file)
    # a list to store the rows of the data frame in
    mutation_rows = []
    # step through the partis data frame row by row
    for row in range(annotations.shape[0]):
        naive_seq = annotations["naive_seq"][row]
        # if no annotation, skip the sequence
        if pd.isnull(naive_seq):
            continue
        # if the mutation rate is higher than we want, skip the sequence
        if annotations["mut_freqs"][row] >= max_mutation_rate:
            continue
        # if partis called an indel and we don't want to use sequences
        # with indels, skip the sequence
        if not use_indel_seqs and annotations["indelfos"][row] != "[[]]":
            continue
        # we call mutations from the indel reversed sequences
        indel_reversed_seq = annotations["indel_reversed_seqs"][row]
        # no indel reversed sequence means that there were no indels
        # to reverse and we should use the input sequence
        if pd.isnull(indel_reversed_seq):
            mutated_seq = annotations["input_seqs"][row]
        else:
            mutated_seq = indel_reversed_seq
        # search along the sequence for mutations (mismatches between
        # naive and indel reversed sequence) and make one row per
        # mutation
        for i in range(len(mutated_seq)):
            if(mutated_seq[i] != naive_seq[i]):
                mutation_rows.append({
                        "mutated_seq": mutated_seq,
                        "naive_seq": naive_seq,
                        "mutated_seq_id": annotations["unique_ids"][row],
                        "mutation_index": i,
                        "gl_base": naive_seq[i],
                        "mutated_base": mutated_seq[i]
                        })
    mutation_df = pd.DataFrame(mutation_rows)
    return(mutation_df)


def process_partis_poly(partis_file, max_spacing, max_mutation_rate=1, use_indel_seqs=True):
    """Gets pairs of neighboring mutations for PolyMotifFinder.

    Keyword arguments:
    partis_file -- A string giving the name of a partis file.
    max_spacing -- The maximum distance allowed between mutations.

    Returns: A data frame. Each row corresponds to a pair of mutations
    within max_spacing of each other and contains the mutated
    sequence, naive sequence, sequence id, mutation indices, germline
    bases, and mutated bases.

    """
    # first get all the single mutations using process_partis
    mutation_df = process_partis(partis_file, max_mutation_rate, use_indel_seqs)
    # a list to hold the rows of the data frame
    poly_mutation_rows = []
    # loop through the sequences
    sequence_ids = set(mutation_df["mutated_seq_id"]) 
    for seq_id in sequence_ids:
        # just the mutations corresponding to the current sequence
        mutation_df_subset = mutation_df[mutation_df.mutated_seq_id == seq_id]
        # get all the pairs of mutations within max_spacing of each other
        poly_mutations = get_pairs(list(mutation_df_subset["mutation_index"]), \
                max_spacing)
        mutated_seq = list(mutation_df_subset["mutated_seq"])[0]
        naive_seq = list(mutation_df_subset["naive_seq"])[0]
        # for each pair of mutations, create a row describing those mutations
        for pm in poly_mutations:
            poly_mutation_rows.append({
                    "mutated_seq": mutated_seq,
                    "naive_seq": naive_seq,
                    "mutated_seq_id": seq_id,
                    "mutation_index": pm,
                    "gl_base": ''.join(naive_seq[i] for i in pm),
                    "mutated_base": ''.join(mutated_seq[i] for i in pm)
                    })
    return(pd.DataFrame(poly_mutation_rows))


def get_pairs(mut_indices, max_spacing):
    """Finds pairs of mutations within max_spacing of each other

    Keyword arguments:
    mut_indices -- A list containing mutation indices.
    max_spacing -- The maximum distance between pairs of mutations.

    Returns: A list containing all pairs of mutations within
    max_spacing of each other.

    """
    # a list to hold the pairs
    pairs = []
    mut_indices_len = len(mut_indices)
    # loop through all the pairs of mutations, put them in the list if
    # they are within max_spacing of each other.
    for i in range(0, mut_indices_len - 1):
        for j in range(i+1, mut_indices_len):
            if np.abs(mut_indices[j] - mut_indices[i]) <= max_spacing:
                pairs.append((mut_indices[i], mut_indices[j]))
            else:
                continue
    return(pairs)
