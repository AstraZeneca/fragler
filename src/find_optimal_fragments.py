#!/usr/bin/env python3

__author__ = "Petr Volkov"

import sys
import pickle
import itertools
from collections import OrderedDict
from typing import Dict, List


def read_fasta(fasta_file):
    """Read a fasta file"""
    seqs = []
    ids = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                ids.append(line.rstrip().lstrip('>'))
            else:
                seqs.append(line.rstrip().upper())

    return (seqs, ids)


def complement(seq):
    map = {"A": "T", "C": "G", "T": "A", "G": "C"}
    return "".join([map[nuc] for nuc in seq])


def reverse(seq):
    seq = seq[::-1]
    return seq


def findall(match, string, start=0, end=-1):
    """Find a match in a string"""
    if end == -1:
        end = len(string)

    res = string.find(match, start, end)
    if res == -1:
        return []
    else:
        return [res] + findall(match, string, res + 1, end)


def fpsp(seq, kmers, last_cut=0, shortest_fragment=300, longest_fragment=1000):
    """Main recursion function"""
    # recursion stops if last cut is closer to the end than shorter fragment
    if last_cut > len(seq) - shortest_fragment:
        return [{}]

    # we are going to look for next cut in the following region:
    start = last_cut + shortest_fragment
    end = min(last_cut + longest_fragment, len(seq) - shortest_fragment)

    # for each of the kmers, find where the cut can be made
    results = [{}]
    for kmer in kmers:
        matches = findall(kmer, seq, start, end)

        # for each match, make a copy of kmer list without the current k-mer
        new_kmers = kmers.copy()
        new_kmers.remove(kmer)

        # for every match, execute a recursion
        possible_cuts = []
        for match in matches:
            # fpsp will return a list of dictionaries
            next_cuts = fpsp(seq, new_kmers, match)
            current_cut = {match: kmer}
            possible_cuts += [{**current_cut, **cut} for cut in next_cuts]

        results += possible_cuts
    return results


def is_redundant_cache(cache):
    if cache == {}:
        return True

    if max(get_fragment_frequenies(cache)) == 1:
        return True
    return False


def rank_cuts_for_input_string(cache, input_string):
    """Create a dictionary with pickled cut as a key
    and a rank of the cut as a value
    """
    cut_ranks = {}
    for fragment, cut_set in cache.items():
        if input_string in cut_set:
            for cut in cut_set[input_string]:

                fragment_importance = len(cut_set)
                if cut in cut_ranks:
                    cut_ranks[cut] -= fragment_importance
                else:
                    cut_ranks[cut] = len(pickle.loads(cut)) - \
                        fragment_importance

    return cut_ranks


def cache_to_cut_dictionary(cache,
                            reorder_cuts_by_rank=True,
                            load_pickle=False):
    """
    Convert a cache to a disctionary with input strings as keys
    and a list of (unique) cuts from the cache that correspond to this
    input string as values.
    """
    cut_dict = OrderedDict()
    for cut_set in cache.values():
        for input_string, cuts in cut_set.items():
            if input_string in cut_dict.keys():
                cut_dict[input_string].update(set(cuts))
            else:
                cut_dict[input_string] = set(cuts)

    for input_string, cuts in cut_dict.items():
        cut_ranks = rank_cuts_for_input_string(cache, input_string)
        cut_dict[input_string] = sorted(
            cuts,
            key=lambda cut: cut_ranks[cut]
        )

    for input_string, cuts in cut_dict.items():
        if load_pickle:
            cut_dict[input_string] = [pickle.loads(cut) for cut in list(cuts)]
        else:
            cut_dict[input_string] = list(cuts)

    for input_string, cuts in cut_dict.items():
        cut_dict[input_string] = cuts[:30]

    # reorder cut dict by the order of input strings
    for key in sorted(cut_dict.keys()):
        cut_dict[key] = cut_dict.pop(key)
    return cut_dict


def cache_diagnostics(cache):
    print()
    print()
    print()
    print("-----------------------------")
    print("-----------------------------")
    print("-----------------------------")
    print("PRINTING CACHE")
    print("-----------------------------")
    print("-----------------------------")
    print("-----------------------------")
    print("Number of input strings in cache: %d" %
          len(get_unique_input_strings(cache)))
    print("-----------------------------")

    print("Number of fragments in cache: %d" %
          len(get_cache_fragments(cache)))
    print("First string: %s" % get_cache_fragments(cache)[0])

    # print number of occurences of all strings in cache
    print("Fragment frequencies")
    print(get_fragment_frequenies(cache))
    print()

    print("-----------------------------")
    print("Number of cuts per fragment")
    [print(get_number_of_cuts_per_fragment(cache, fragment), end=' ')
     for fragment in get_cache_fragments(cache)]

    print()
    print("-----------------------------")
    print("Cache complexity")
    cache_complexity = number_of_combinations_to_test(cache)
    print(cache_complexity)


def is_a_perfect_fragment(cache, fragment):
    # check if a fragment comes from all srings in cache
    return set(cache[fragment].keys()) == set(get_unique_input_strings(cache))


def find_cut_in_cache(cache, cut):
    """Bolean;
    checks if a cut (pickled) is in cache (pickeled)
    """
    for fragment, cut_set in cache.items():
        for input_string, cuts in cut_set.items():
            if cut in cuts:
                return True
    return False


def number_of_combinations_to_test(cache):
    product = 1
    for input_string in get_unique_input_strings(cache):
        product *= len(set(get_cuts_for_input_string(cache, input_string)))
    return product


def get_number_of_cuts_per_fragment(cache, fragment):
    return sum([len(cuts) for cuts in cache[fragment].values()])


def get_cuts_for_input_string(cache, input_index):
    """Return a list of cuts for an input string
    (characterized by a number, so 0 - first input string) from the cache"""
    res = []
    for cut_set in cache.values():
        if input_index in cut_set.keys():
            res += cut_set[input_index]
    return res


def get_first_cut(cache):
    """Only for diagnostic purposes"""
    return cache[list(cache.keys())[0]][0][0]


def collapse_cache(cache):
    """For each of the strings in cache, do the following:
    Convert each protein index into a list of pickle strings instead of cuts
    """
    input_cuts = {}

    for fragment, cutinfo in cache.items():
        for input_index, cuts in cutinfo.items():
            if input_index in input_cuts.keys():
                input_cuts[input_index].update({fragment: cuts})
            else:
                input_cuts[input_index] = {fragment: cuts}

    # for each string
    for index, fragments in input_cuts.items():
        # for each fragment, find it has overlapping cuts with other fragments
        # (cause if not - we cak take any) OBS - not any, the shortest
        for fragment_i, cuts_i in fragments.items():
            res = []
            pickle_i = [pickle.dumps(cut) for cut in cuts_i]

            for fragment_j, cuts_j in fragments.items():
                pickle_j = [pickle.dumps(cut) for cut in cuts_j]
                if fragment_i != fragment_j:
                    overlap = set(pickle_i).intersection(set(pickle_j))
                    if overlap:
                        res += list(overlap)

            if res == []:
                res = [cuts_i[0]]
            input_cuts[index][fragment_i] = res
    return input_cuts


def cuts_to_native(cache):
    """
    Convert pickle cuts in a cache object
    into native python dict-based strings
    """
    for fragment, cut_set in cache.items():
        for input_string, cuts in cut_set.items():
            cache[fragment][input_string] = [pickle.loads(cut) for cut in cuts]

    return cache


def thinner_cache(cache):
    """We will find fragments that makes sence to keep for furher analysis"""

    # step one - remove potential fragments that appear only in 1 string
    fragments_to_keep = {k: v for k, v in cache_sizes(cache).items() if v > 1}

    # the key idea of the algorithm - keep only the shortest of
    # the nexted fragments
    filtered_fragments = filter_nested_fragments(
        list(fragments_to_keep.keys())
    )
    filtered_fragments += filter_nested_fragments(
        list(fragments_to_keep.keys()),
        keep_largest=False
    )
    fragments_to_keep = {k: v for k, v in fragments_to_keep.items()
                         if k in filtered_fragments}

    return filter_fragments(cache, fragments_to_keep)


def filter_fragments_below_frequency(cache, frequency=1):
    """Remove from the cache fragments that occur
    with this frequency or less
    """
    fragments_to_keep = {k: v for k, v in cache_sizes(cache).items()
                         if v > frequency}
    return {k: v for k, v in cache.items() if k in fragments_to_keep.keys()}


def filter_unitary_fragments(cache):
    """
    Takes a cache and returs a cache without fragments that
    only occure once (come from only 1 input string)
    """
    return filter_fragments_below_frequency(cache, 1)


def filter_fragments(cache, fragments_to_keep):
    """
    Filter fragments from cache
    """
    return {k: v for k, v in cache.items()
            if k in fragments_to_keep}


def filter_nested_fragments(strs, keep_largest=True):
    """
    Remove fragments which contain other fragments
    so that only the longest nested string remains.
    Can't find a better way to do it now - refactor.
    """
    strs = list(strs)
    strs.sort()

    inds_to_remove = []
    for i, pattern in enumerate(strs):
        for j, target in enumerate(strs):
            if i != j and target.find(pattern) != -1:
                if keep_largest:
                    inds_to_remove.append(i)
                else:
                    inds_to_remove.append(j)
    return [str for i, str in enumerate(strs) if i not in inds_to_remove]


def find_potential_split_positions(seq,
                                   kmers,
                                   shortest_fragment=300,
                                   longest_fragment=1000):
    """
    Wrapper function to call the recurive algorithm of finding
    the potential cuts performs some additional checks on the
    fragments returned by the recursion
    """
    cuts = fpsp(seq,
                kmers,
                shortest_fragment=shortest_fragment,
                longest_fragment=longest_fragment)

    # recursion returns cuts that have the last piece
    # longer than the longest fragment
    # or shorter that the shortest fragment
    # so we filter the recursion cutsults
    if cuts[0] == {}:
        cuts.pop(0)

    cuts = [cut for cut in cuts
            if max(cut.keys()) > len(seq) - longest_fragment
            and max(cut.keys()) < len(seq) - shortest_fragment]

    # if the string is short enough, add a dummy fragment to the end:
    if (len(seq) < longest_fragment):
        cuts += [{}]

    return cuts


def nfragments_per_cut_combination(seqs: List[str],
                                   cuts: List[Dict[int, str]]):
    """ This is an objective funciton that takes a list of sequences
    and a list of cuts (one cut for each string)
    and returns a number of fragments """
    if len(seqs) != len(cuts):
        raise("Number of cuts in a combination of cuts should equal the \
               number of seqs to which the cuts are proposed")

    fragments = []

    for i in range(0, len(seqs)):
        a = [0] + list(cuts[i].keys())
        b = list(cuts[i].keys()) + [len(seqs[i])]
        fragments += [[y - x, x, y, seqs[i]] for x, y in zip(a, b)]

    fragments = sorted(fragments)
    identical_fragments = 0

    for i in range(0, len(fragments) - 1):
        if fragments[i][0] == fragments[i + 1][0]:
            seq1 = fragments[i][3]
            seq2 = fragments[i + 1][3]

            fragment1 = seq1[(fragments[i][1]): (fragments[i][2])]
            fragment2 = seq2[fragments[i+1][1]: fragments[i+1][2]]

            if fragment1 == fragment2:
                identical_fragments += 1
    return len(fragments) - identical_fragments


def cut_to_fragment_coordinates(seq: str, cut: Dict[int, str]):
    """ get fragment seuqneces for a single cut """
    starts = [0] + list(cut.keys())
    ends = list(cut.keys()) + [len(seq)]
    return list(zip(starts, ends))


def coordinate_pair_to_seq(seq, coord):
    """Takes a touple of coordinates and
    returns a string between this coordinates
    """
    return seq[coord[0]:coord[1]]


def fragment_coordinates_to_fragment_seqs(seq, coords):
    """Takes a list of fragment coordinates
    and returns a list of fragment strings
    """
    return [coordinate_pair_to_seq(seq, coord) for coord in coords]


def cut_to_fragment_seqs(seq: str, cut: Dict[int, str]):
    """Takes a sequence and a cut
    and returns a list of fragment sequences
    """
    coords = cut_to_fragment_coordinates(seq, cut)
    return fragment_coordinates_to_fragment_seqs(seq, coords)


def cuts_to_fragment_seqs(seq: str, cuts: List[Dict[int, str]]):
    """Takes a sequence and list of cuts (usually all potential cuts),
    and returns a list (dict?) of fragmets
    """
    cut_fragments = [cut_to_fragment_seqs(seq, cut) for cut in cuts]
    return cut_fragments


def fragment_seqs_to_fragment_cache(fragment_seqs_per_cut,
                                    cuts,
                                    add_left_overhangs=True):
    """Takes a list of cuts and corresponding fragment sequences
    and builds a dictionary for the fragment sequences, where:

    keys - fragment sequences
    values - lists of cuts where this sequences come from
    """
    fragment_cache = OrderedDict()
    for cut, fragment_seqs in zip(cuts, fragment_seqs_per_cut):
        for ind, fragment_seq in enumerate(fragment_seqs):

            # add overhang sequence to the left side of the cut
            if add_left_overhangs and ind != len(fragment_seqs) - 1:
                fragment_seq += fragment_seqs[ind+1][:4]

            if fragment_seq in fragment_cache.keys():
                fragment_cache[fragment_seq] += [cut]
            else:
                fragment_cache[fragment_seq] = [cut]
    return fragment_cache


def input_seqs_to_fragment_cache(seqs,
                                 kmers,
                                 shortest_fragment=300,
                                 longest_fragment=1000):
    """
    Takes a set of input sequences and returns a sequence cache, same as in
    fragment_seqs_to_fragment_cache, but including the number of input string
    """

    input_cache = {}

    for index, seq in enumerate(seqs):
        cuts = find_potential_split_positions(seq,
                                              kmers,
                                              shortest_fragment,
                                              longest_fragment)
        fragment_seqs = cuts_to_fragment_seqs(seq, cuts)
        fragment_cache = fragment_seqs_to_fragment_cache(fragment_seqs, cuts)
        for fragment, cuts in fragment_cache.items():
            if fragment in input_cache.keys():
                input_cache[fragment][index] = cuts
            else:
                input_cache[fragment] = {}
                input_cache[fragment][index] = cuts

    return input_cache


def cuts_to_pickle(cache):
    '''Convert native python dict-based cuts in a
    cache object into pickle strings
    '''
    for fragment, cut_set in cache.items():
        for input_string, cuts in cut_set.items():
            cache[fragment][input_string] = [pickle.dumps(cut) for cut in cuts]

    return cache


def cache_sizes(cache):
    """This function converts cache in a simple dictionary of sizes;
    sizes in this case means the number of strings in which
    a potential fragment can appear"""
    return {k: len(v) for k, v in cache.items()}


def get_cache_fragments(cache):
    return list(cache.keys())


def get_fragment_frequency(cache, fragment):
    return len(cache[fragment].keys())


def get_unique_input_strings(cache):
    res = set()
    for cut_set in cache.values():
        res.update(cut_set.keys())

    return list(res)


def iterate_through_cut_dictionary(cut_dict):
    return itertools.product(*cut_dict.values())


def get_fragment_frequenies(cache):
    return [get_fragment_frequency(cache, fragment)
            for fragment in get_cache_fragments(cache)]


def get_fragment_lengths(cache):
    return [len(fragment) for fragment in get_cache_fragments(cache)]


def order_cache(cache):
    """ This function makes sure the cache is ordered in the order of
    1 - fragment rank (decending)
    2 - fragmnet length (descending)
    """
    ordered_cache = OrderedDict()

    lengths = get_fragment_lengths(cache)
    frequencies = get_fragment_frequenies(cache)
    fragments = get_cache_fragments(cache)

    sort_order = list(zip(lengths, frequencies, fragments))
    sort_order = sorted(sort_order, key=lambda x: x[0], reverse=True)
    sort_order = sorted(sort_order, key=lambda x: x[1], reverse=True)

    for (l, q, f) in sort_order:
        ordered_cache[f] = cache[f]

    return ordered_cache


def merge_cut_sets(first, second):
    """ Take two cut sets and merge them.
    If an input string is in both, and overlap is nonempty -
        take an overlap + cuts shorter than overlap shortest cut
    If an input string is in both, and overlap is emtpy -
        take cuts from the first + cuts shorter than in the first
    If an input string is in one but not another -
        take cuts from the existing one
    """
    new_cut_set = OrderedDict()
    overlapping_input_strings = set(
        first.keys()).intersection(set(second.keys())
                                   )

    if overlapping_input_strings == set():
        new_cut_set = {**first, **second}
        return new_cut_set

    for s in overlapping_input_strings:
        overlap = set(first[s]).intersection(set(second[s]))
        if overlap != set():
            min_length_to_keep = min([len(pickle.loads(cut))
                                      for cut in overlap])
            keep = set(filter(
                lambda x: len(pickle.loads(x)) < min_length_to_keep, first[s]))
            new_cut_set[s] = overlap
            new_cut_set[s].update(keep)
        else:
            min_length_to_keep = min([len(pickle.loads(cut))
                                      for cut in first[s]])
            keep = set(filter(
                lambda x: len(pickle.loads(x)) < min_length_to_keep, second[s])
                      )

            # we always take the first if there is no overlap
            new_cut_set[s] = set(first[s])
            new_cut_set[s].update(keep)

    for s in set(first.keys() - second.keys()):
        new_cut_set[s] = first[s]

    for s in set(second.keys() - first.keys()):
        new_cut_set[s] = second[s]

    return new_cut_set


def the_algoritm(cache):
    """
    Algorithm to recursively merge cut sets
    """
    # start with the longest fragment of the highest rank
    candidate = list(cache.keys())[0]
    candidate_cut_set = cache[candidate]

    # and then merge cutsets one by one
    for fragment, cut_set in cache.items():
        candidate_cut_set = merge_cut_sets(candidate_cut_set, cut_set)

    return candidate_cut_set


def brute_force(cutdict):
    """
    Attempt at brute force - not used
    """
    i = 1
    min = 1000
    best_cut = []
    for cut in iterate_through_cut_dictionary(cutdict):
        i += 1
        if i % 100000 == 0:
            print("Iteration %d" % i)

        if i > 100000000:
            break

        res = nfragments_per_cut_combination(list(cutdict.keys()), list(cut))
        if res < min:
            best_cut = list(cut)
            min = res

    return list(cutdict.keys()), list(best_cut)


def print_cut(seqs, cuts, fasta, insert_sequence="", add_left_overhangs=True):
    fasta_dict = {}
    for seq, seq_id in zip(fasta[0], fasta[1]):
        fasta_dict[seq] = seq_id

    print("[Proposed fragmentation]")
    fragment_frequencies = OrderedDict()
    for seq, cut in zip(seqs, cuts):
        print(fasta_dict[seq], end=": ")

        starts = [0] + list(cut.keys())
        ends = list(cut.keys()) + [len(seq)]
        fragments = [seq[s:e] for s, e in zip(starts, ends)]

        for i, fragment in enumerate(fragments):
            # add overhang sequence to the left side of the cut
            if add_left_overhangs and i != len(fragments) - 1:
                fragment += fragments[i+1][:4]
            if i != len(fragments) - 1:
                fragment += reverse(complement(insert_sequence))
            if i != 0:
                fragment = insert_sequence + fragment[:4] + fragment[4:]
            if fragment in fragment_frequencies.keys():
                ind = list(fragment_frequencies.keys()).index(fragment) + 1
                fragment_frequencies[fragment] += 1
            else:
                fragment_frequencies[fragment] = 1
                ind = len(fragment_frequencies)
            print("[%d]" % ind, end=' ')
        print()

    print()
    print("[Fragment dictionary]")

    for i, fragment in enumerate(fragment_frequencies):
        print("%d: \"%s\"" % (i+1, fragment))


def run(fasta,
        kmers_set,
        shortest_fragment,
        longest_fragment,
        insert_sequence):
    seqs = fasta[0]
    cache = input_seqs_to_fragment_cache(
        seqs, kmers_set, shortest_fragment, longest_fragment
    )
    cache = cuts_to_pickle(cache)
    cache = order_cache(cache)

    candidate_cut_set = the_algoritm(cache)
    cut_dict = {}
    for ind, cut in candidate_cut_set.items():
        cut_dict[seqs[ind]] = [pickle.loads(c) for c in list(cut)]

    best_cut = brute_force(cut_dict)
    print_cut(*best_cut, fasta, insert_sequence, add_left_overhangs=True)


def main():
    """Main function"""
    fasta_file = sys.argv[1]
    kmers_set = sys.argv[2]
    shortest_fragment = int(sys.argv[3])
    longest_fragment = int(sys.argv[4])
    kmers_set = kmers_set.split(',')
    insert_sequence = sys.argv[5]

    fasta = read_fasta(fasta_file)
    run(fasta, kmers_set, shortest_fragment, longest_fragment, insert_sequence)


if __name__ == "__main__":
    main()
