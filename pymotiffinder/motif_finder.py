import pandas as pd
import numpy as np
from process_partis import process_partis
from process_partis import process_partis_poly
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC


def motif_finder(partis_file, reference_fasta, k,
                 kmer_dict=None,
                 reverse_complement=False,
                 max_mutation_rate=1,
                 use_indel_seqs=True,
                 return_dict=False,
                 unique_mutations=False):
    """Matches mutations to potential gene conversion donors.

    Keyword arguments:
    partis_file -- A file containing partis output describing the
    germline and mature sequences.
    reference_fasta -- A fasta file containing the sequences in the
    donor gene set.
    k -- kmer size.
    reverse_complement -- If True, also look for gene conversion
    templates on the reverse complement.
    max_mutation_rate -- Remove any sequences from the partis file
    that have a mutation rate above max_mutation_rate.
    use_indel_seqs -- If False, remove any seqeunces with indels from
    partis_file.

    Returns: A data frame with the query name, mutation index, and the
    number of alignments in the reference set explaining that mutation.

    """
    mutations = process_partis(partis_file,
                                max_mutation_rate=max_mutation_rate,
                                use_indel_seqs=use_indel_seqs)
    n_mutations = get_n_mutations(mutations, unique=unique_mutations)
    if kmer_dict is None:
        kmer_dict = make_kmer_dict_from_fasta(reference_fasta,
                                              k,
                                              reverse_complement=reverse_complement)
    imf_out = indexed_motif_finder(mutations, kmer_dict, k)
    hits, total_mutations = templated_number(imf_out, dale_method=False)
    if not return_dict:
        return(imf_out,  hits / n_mutations)

    return (imf_out, hits / n_mutations, kmer_dict)

def poly_motif_finder(partis_file, reference_fasta, k,
                      reverse_complement=False,
                      max_mutation_rate=1,
                      use_indel_seqs=True,
                      kmer_dict=None,
                      dale_method=False):
    """Matches mutations to potential gene conversion donors.

    Keyword arguments:
    partis_file -- A file containing partis output describing the
    germline and mature sequences.
    reference_fasta -- A fasta file containing the sequences in the
    donor gene set.
    k -- kmer size.
    reverse_complement -- If True, also look for gene conversion
    templates on the reverse complement.
    max_mutation_rate -- Remove any sequences from the partis file
    that have a mutation rate above max_mutation_rate.
    use_indel_seqs -- If False, remove any seqeunces with indels from
    partis_file.

    Returns: A data frame with the query name, mutation index, and the
    number of alignments in the reference set explaining that mutation.

    """
    ## create the kmer dictionary if needed
    if kmer_dict is None:
        kmer_dict = make_kmer_dict_from_fasta(reference_fasta,
                                              k,
                                              reverse_complement=reverse_complement)
    ## get the groups of mutations to be explained
    poly_mutations = process_partis_poly(partis_file, max_spacing=k-1,
                               max_mutation_rate=max_mutation_rate,
                               use_indel_seqs=use_indel_seqs)
    ## data frame containing sets of mutations and information about their templates
    imf_out = indexed_motif_finder(poly_mutations, kmer_dict, k)
    if dale_method:
        hits, coverage_denom = templated_number(imf_out, dale_method=True)
        return (imf_out, hits / coverage_denom)
    else:
        hits, _ = templated_number(imf_out, dale_method=False)
        ## get the number of mutations
        single_mutations = process_partis(partis_file,
                                            max_mutation_rate=max_mutation_rate,
                                            use_indel_seqs=use_indel_seqs)
        n_mutations = single_mutations.shape[0]
        return (imf_out, hits / n_mutations)


def seed_starts(idx, seed_len, seq_len):
    """Finds starting positions for windows around mutations

    Keyword arguments:
    idx -- Either a single number (for a single mutation) or a tuple or list
    containing the indices of multiple mutations.
    seed_len -- The length of the window.
    seq_len -- The length of the sequence the mutations occur in.

    Returns: A pair with the lowest and highest indices a window containing
    the mutations can start at.

    """
    idx_hi = np.max(idx)
    idx_lo = np.min(idx)
    min_start = max(0, idx_hi - seed_len + 1)
    max_start = min(seq_len - seed_len, idx_lo)
    return((min_start, max_start))


def make_kmer_dict_from_fasta(fasta_file, k, reverse_complement=False):
    """Make a kmer dictionary from a fasta file

    Keyword arguments:

    fasta_file -- A file containing the sequences of potential
    conversion donors.
    k -- k-mer size.
    reverse_complement -- If True, also look for gene conversion
    templates on the reverse complement.

    Returns: A dictionary keyed by k-mers, mapping to sets of (name,
    location) pairs describing where the k-mer occured.

    """
    refs = [r for r in SeqIO.parse(fasta_file, "fasta", alphabet=IUPAC.unambiguous_dna)]
    kmer_dict = make_kmer_dictionary(refs, k, reverse_complement=reverse_complement)
    return(kmer_dict)


def make_kmer_dictionary(references, k, reverse_complement=False):

    """Make a dictionary mapping kmers to sequences they appear in

    Keyword arguments:
    references -- A list of SeqRecord objects, with names and
    sequences, sequences in IUPAC.unambiguous_dna.
    k -- The size of the k-mer.
    reverse_complement -- If True, also look for gene conversion
    templates on the reverse complement.

    Returns: A dictionary keyed by k-mers, mapping to sets of (name,
    sequence, location) pairs describing where the k-mer occurred.

    """
    d = {}
    if reverse_complement:
        references_rc = list()
        for ref in references:
            ref_rc = SeqRecord(name=ref.name + "_rc", seq=ref.reverse_complement().seq)
            references_rc.append(ref_rc)
        references = references + references_rc
    for ref in references:
        seq = str(ref.seq)
        seq_len = len(seq)
        for start in range(0, seq_len - k + 1):
            kmer = seq[start:(start + k)]
            if kmer in d.keys():
                d[kmer].add((ref.name, ref.seq, start))
            else:
                d[kmer] = set([(ref.name, ref.seq, start)])
    return d


def n_alignments_per_mutation(mutations, kmer_dict, k):
    """ Find the number of unique alignments in the reference set for each
    mutation

    Keyword arguments:
    mutations -- A data frame created by process_partis
    kmer_dict -- A kmer dictionary created by make_kmer_dictionary
    k -- The k used in the kmer dictionary

    Returns: A data frame with the query name, mutation index, and the
    number of alignments in the reference set explaining that mutation.
    """

    # a data frame describing the matches for each mutation in the references
    imf = indexed_motif_finder(mutations, kmer_dict, k)
    # dictionary that will hold the mutations and how many templates they have
    count_dict = {}
    # Each row in imf describes a mutation. If there is no template
    # for that mutation, the reference_alignment column is
    # np.nan. Otherwise, each row gives the location of one of the
    # templates for that mutation.
    for index, row in imf.iterrows():
        query = row["query_name"]
        query_index = row["query_mutation_index"]
        if np.isnan(row["reference_alignment"]):
            increment = 0
        else:
            increment = 1
        # mutations are described as a pair with the query name and
        # the location of the mutation
        if (query, query_index) in count_dict.keys():
            count_dict[(query, query_index)] += increment
        else:
            count_dict[(query, query_index)] = increment
    # build a DataFrame with query names, indices, and number of alignments
    rows = []
    for (query, query_index) in count_dict.keys():
        rows.append({
                "query_mutation_index": query_index,
                "query_name": query,
                "n_alignments": count_dict[(query, query_index)]
                })
    return(pd.DataFrame(rows))


def indexed_motif_finder(mutations, kmer_dict, k):
    """Find matches around a set of mutations.

    Keyword arguments:
    mutations -- A data frame containing the mutated sequences and
    mutation indices.
    kmer_dict -- A dictionary indexed by k-mers giving the sequences
    they appear in.
    k -- The k used in the kmer dictionary

    Returns: A pandas DataFrame containing the query sequence, the
    indices of the mutation(s) in the query, the name and sequence of
    the reference with a match, and a reference alignment. This is the
    position in the reference that matches the mutation in the query
    (if we are loking at single mutations) or the position in the
    reference matching the left-most mutation in a set of mutations in
    the query.

    """
    row_list = []
    for index, row in mutations.iterrows():
        sequence_list = list(row["mutated_seq"])
        q = "".join(sequence_list)
        q_id = row["mutated_seq_id"]
        seq_len = len(q)
        mut_idx = row["mutation_index"]
        # search through all the windows around the mutation and check whether
        # they occur in the references
        found_match = False
        (min_start, max_start) = seed_starts(mut_idx, k, seq_len)
        for start in range(min_start, max_start + 1):
            # create the seed
            seed = q[start:(start + k)]
            mut_offset = np.min(mut_idx) - start
            if seed in kmer_dict:
                for (ref_name, ref_seq, ref_idx) in kmer_dict[seed]:
                    row_list.append({
                            "query_sequence": q,
                            "query_name": q_id,
                            "query_mutation_index": mut_idx,
                            "naive_sequence": row["naive_seq"],
                            "mutated_base": row["mutated_base"],
                            "reference_name": ref_name,
                            "reference_sequence": str(ref_seq),
                            "reference_alignment": ref_idx + mut_offset
                            })
                found_match = True
        # if there wasn't a match, we still put the sequence in DataFrame,
        # with np.nan as the value for reference_alignment
        if not found_match:
            row_list.append({
                    "query_sequence": q,
                    "query_name": q_id,
                    "query_mutation_index": mut_idx,
                    "naive_sequence": row["naive_seq"],
                    "mutated_base": row["mutated_base"],
                    "reference_name": "",
                    "reference_sequence": "",
                    "reference_alignment": np.nan
                    })
    return pd.DataFrame(row_list).drop_duplicates().reset_index(drop=True)


def extend_matches(df):
    """Extends matches from indexed_motif_finder. Adds a columns for left-most 
    and right-most match indices and the match extent"""

    row_count = df.shape[0]
    df["match_extent"] = pd.Series([0] * row_count, index=df.index)
    df["query_left_idx"] = pd.Series([0] * row_count, index=df.index)
    df["query_right_idx"] = pd.Series([0] * row_count, index=df.index)
    for row in range(0, row_count):
        if df.loc[row, "reference_sequence"] == "":
            continue
        left = 0
        right = 0
        query_idx = int(np.min(df.loc[row, "query_mutation_index"]))
        ref_idx = int(df.loc[row, "reference_alignment"])
        query_seq = df.loc[row, "query_sequence"]
        ref_seq = df.loc[row, "reference_sequence"]
        while True:
            if query_idx - left - 1 < 0:
                break
            elif ref_idx - left - 1 < 0:
                break
            elif ref_seq[ref_idx - left - 1] == query_seq[query_idx - left - 1]:
                left += 1
            else:
                break
        while True:
            if ref_idx + right + 1 >= len(ref_seq):
                break
            elif query_idx + right + 1 >= len(query_seq):
                break
            elif ref_seq[ref_idx + right + 1] == \
                    query_seq[query_idx + right + 1]:
                right += 1
            else:
                break
        df.loc[row, "match_extent"] = left + right + 1
        df.loc[row, "query_left_idx"] = query_idx - left
        df.loc[row, "query_right_idx"] = query_idx + right

def templated_number(df, dale_method=False):
    """Number of templated mutations from PolyMotifFinder.

    Keyword arguments:
    df -- A pandas DataFrame created by indexed_motif_finder, with query_mutation_index giving a set of indices corresponding to mutations that are explained or not by templated mutagenesis.
    Returns: The number of mutations that have a template.
    If dale_method=True, this is the number of unique mutations, where a unique mutation is one that corresponds to the same base and the same germline sequence.
    If dale_method=False, this is the total number of mutations, so if we see the same mutation in more than one mutated sequence and it has a template both times, we count it twice.
    """
    already_seen_set = set()
    ## keyed by mutated sequence, mutation index pairs
    ## values are a pair, first element the sequence, index element of the mutation matrix, second element the sequence, index element of the scoring matrix
    scoring_and_mutation_matrices = {}
    for (index, row) in df.iterrows():
        templated = not np.isnan(row["reference_alignment"])
        seen_before = (row["naive_sequence"], row["query_mutation_index"], row["mutated_base"]) in already_seen_set
        set_scoring_and_mutation_dict(row["query_sequence"], row["query_mutation_index"], templated, seen_before, scoring_and_mutation_matrices)
        # if we're trying to count unique mutations, add the pair to the set of already seen mutation pairs
        if dale_method:
            already_seen_set.add((row["naive_sequence"], row["query_mutation_index"], row["mutated_base"]))

    hit_number = sum([value[0] and (not value[1]) for value in scoring_and_mutation_matrices.values()])
    total_mutations = sum((not value[1]) for value in scoring_and_mutation_matrices.values())
    return float(hit_number), float(total_mutations)

def set_scoring_and_mutation_dict(mutated_seq, mutation_index, is_templated, was_seen_before, scoring_and_mutation_dict):
    if type(mutation_index) is not tuple:
        mutation_index = [mutation_index]
    for m in mutation_index:
        if (mutated_seq, m) in scoring_and_mutation_dict:
            old_value = scoring_and_mutation_dict[(mutated_seq, m)]
        else:
            old_value = [False,False]
        scoring_and_mutation_dict[(mutated_seq, m)] = [old_value[0] or is_templated, old_value[1] or was_seen_before]


def get_n_mutations(mutation_df, unique=False):
    """Counts the number of mutations

    Keyword arguments:
    mutation_df: The output from process_partis
    unique: If True, counts the number of "unique" mutations, so if we see the same mutation away from the naive sequence in two different mutated sequences, it only counts once. Otherwise, counts the total number of mutations.

    Returns: The number of mutations
    """
    ## we have one row per mutation, so if we're not collapsing on the naive sequence the number of mutations is the number of rows in mutation_df
    if not unique:
        return mutation_df.shape[0]
    ## otherwise we make a set to store the mutations
    mutation_set = set()
    for index, row in mutation_df.iterrows():
        naive_seq = row["naive_seq"]
        mutation_index = row["mutation_index"]
        mutation_identity = row["mutated_base"]
        mutation_set.add((naive_seq, mutation_index, mutation_identity))
    return len(mutation_set)



def likelihood_given_gcv(partis_file, kmer_dict, k, max_mutation_rate, use_indel_seqs):
    """Finds the likelihood of mutations conditional on being due to gcv

    Keyword arguments:
    partis_file -- A partis csv describing the mutations.
    kmer_dict -- A kmer dictionary describing the references.
    k -- The minimum match length for gcv tracts.

    Returns: A data frame giving the probability of seeing each
    observed mutation. Mutations are described by the name of the
    query sequence and the position of the mutation in that query
    sequence.

    """
    bases = ["A", "C", "G", "T"]
    mut_df = process_partis(partis_file, max_mutation_rate=max_mutation_rate, use_indel_seqs=use_indel_seqs)
    # make a data frame containing all the mutations we didn"t see
    unobs_mut_rows = []
    for index, row in mut_df.iterrows():
        for b in bases:
            if b not in set([row["gl_base"], row["mutated_base"]]):
                r = row.copy()
                unseen_seq = list(r["mutated_seq"])
                unseen_seq[r["mutation_index"]] = b
                r["mutated_seq"] = "".join(unseen_seq)
                r["mutated_base"] = b
                unobs_mut_rows.append(r)
    unobs_mut_df = pd.DataFrame(unobs_mut_rows)

    # run motif finder on the observed and unobserved mutations
    motifs_obs = n_alignments_per_mutation(mut_df, kmer_dict, k)
    motifs_unobs = n_alignments_per_mutation(unobs_mut_df, kmer_dict, k)
    obs_and_unobs = pd.merge(motifs_obs, motifs_unobs,
                             how="outer",
                             on=["query_mutation_index", "query_name"],
                             validate="one_to_one")

    # get the probabilities of seeing the observed
    # mutations. n_alignments_x is the number of alignments for
    # observed mutations because motifs_obs was in the first position
    # in pd.merge
    def get_prob(row):
        n_obs = row["n_alignments_x"]
        n_unobs = row["n_alignments_y"]
        if n_obs + n_unobs == 0:
            return(np.nan)
        return n_obs / (n_obs + n_unobs + 0.)

    obs_and_unobs["prob"] = obs_and_unobs.apply(get_prob, axis=1)
    return(obs_and_unobs)


def per_base_alignments(partis_file, kmer_dict, k, max_mutation_rate, use_indel_seqs):
    """Finds the number of templates for each potential base at each mutated site.

    Keyword arguments:
    partis_file -- A partis csv describing the mutations.
    kmer_dict -- A kmer dictionary describing the references.
    k -- The minimum match length for gcv tracts.

    Returns: A data frame giving the probability of seeing each
    observed mutation. Mutations are described by the name of the
    query sequence and the position of the mutation in that query
    sequence.

    """
    bases = ["A", "C", "G", "T"]
    mut_df = process_partis(partis_file, max_mutation_rate=max_mutation_rate, use_indel_seqs=use_indel_seqs)
    output_rows = []
    # make a data frame containing all the mutations we didn"t see
    for index, row in mut_df.iterrows():
        output_row = row.copy()
        for b in bases:
            r = row.copy()
            unseen_seq = list(r["naive_seq"])
            unseen_seq[r["mutation_index"]] = b
            r["mutated_seq"] = "".join(unseen_seq)
            r["mutated_base"] = b
            output_row[b] = n_alignments_per_mutation(pd.DataFrame([r]), kmer_dict, k)["n_alignments"].item()
        output_rows.append(output_row)

    return(pd.DataFrame(output_rows))


def templates_per_base(partis_file, kmer_dict, k):
    """

    Keyword arguments:
    partis_file -- A partis csv describing the mutations.
    kmer_dict -- A kmer dictionary describing the references.
    k -- The minimum match length for gcv tracts.

    Returns: A data frame where each row contains the number of
    templates for a putative mutation. The mutation is described in
    terms of its germline base, the target base, the id of the mutated
    sequence, the location within that sequence of the mutation, and
    the mutation that actually occurred.

    """
    bases = ["A", "C", "G", "T"]
    mut_df = process_partis(partis_file)
    output = []
    # make a data frame containing all the mutations we didn"t see
    for index, row in mut_df.iterrows():
        for b in bases:
            r = row.copy()
            r["mutated_seq"] = make_mutated_sequence(r["naive_seq"], r["mutation_index"], b)
            n = n_alignments_per_mutation(r, kmer_dict, k)
            output_row = {"gl_base": r["gl_base"],
                          "template_base": b,
                          "mutated_seq_id": r["mutated_seq_id"],
                          "mutation_index": r["mutation_index"],
                          "true_mutation": r["mutated_base"]}
            output.append(output_row)
    return(pd.DataFrame(output))


def make_mutated_sequence(naive, index, base):
    """

    Keyword arguments:
    naive -- A string giving the naive sequence.
    index -- The position in the naive sequence to be changed.
    base -- The target base.

    Returns: A string with the mutated sequence.

    """
    seq_as_list = list(naive)
    seq_as_list[index] = base
    return("".join(seq_as_list))
