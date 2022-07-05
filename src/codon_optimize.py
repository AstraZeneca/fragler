__author__ = "Dean Sumner"

import sys
import requests
import json
import re
import numpy as np


def read_fasta(file):
    """Function for reading in multiple protein sequences from fasta format.
    Args:
        file (str): string containing the path to fasta file.
    Returns:
        sequences (dict): codon optimized form of DNA with specified arguments.
    """
    seqs = {}
    for line in open(file).readlines():
        if line.startswith(">"):
            seqid = line[1:].strip()
            seqs[seqid] = ""
        else:
            seqline = line.strip()
            seqline = seqline.replace(" ", "")
            seqs[seqid] += seqline
    return seqs


def find_common_substrings(seqa, seqb, min_fragment_size):
    """
    Function which acts as the base operation for collon substring algorithm.
    Args:
        seqa (str): string containing the path to fasta file.
        seqa (str): string containing the path to fasta file.
        min_fragment_size (int): int specifying minimum size for each
                                 identified substring
    Returns:
        common (list): list of common sequences found
    """

    apos = 0
    bpos = 0
    alen = len(seqa)
    blen = len(seqb)
    common = []

    while apos < alen - min_fragment_size:
        aword = seqa[apos: apos + min_fragment_size]

        hit = seqb[bpos:].find(aword)
        if "[" not in aword and "]" not in aword and hit >= 0:

            bpos += hit
            i = 0
            while (
                apos + i < alen
                and bpos + i < blen
                and seqa[apos + i] != "["
                and seqa[apos + i] == seqb[bpos + i]
            ):
                i += 1
            i += -1
            common.append([apos, bpos, i + 1])
            apos += i
            bpos += i
        else:
            apos += 1
    return common


def get_common_seqs(seqs, passage, common_seqs_used, minLength=100):
    """Function which acts as the wrapper for substring algorithm.
    Args:
        seqs (dict): string containing the sequences to be searched.
        passage (int): int containing the current passage or epoch.
        minLength (int): int specifying minimum size for each
    identified substring.
    Returns:
        boolean (True): placeholder for next epoch.
        seqs (dict): dict containing the updated sequences with resepctive ID's
        commonSeqsUsed (list): list containing all of the common
    sequences identified.
    """
    new_seqs = seqs
    common_seqs = []
    for seqa in seqs:
        for seqb in seqs:
            if seqa != seqb:

                sequencea = seqs[seqa]
                sequenceb = seqs[seqb]
                commons = find_common_substrings(
                    sequencea,
                    sequenceb,
                    minLength)
                for common in commons:
                    if (
                        sequencea[common[0]: common[0] + common[2]] not in
                            common_seqs and common != []
                    ):
                        common_seqs.append(sequencea[
                            common[0]: common[0] + common[2]
                        ])

    if len(common_seqs) == 0:
        return (False, new_seqs, common_seqs_used)

    common_seqs.sort(key=len)
    new_seqs = {}
    ci = 0
    for cseq in common_seqs:
        for seqname, seq in seqs.items():
            sequence = seq
            if cseq in sequence:
                seqs[seqname] = sequence.replace(
                    cseq, "[" + str(passage) + "_" + str(ci) + "]"
                )
                common_seqs_used["[" +
                                 str(passage) +
                                 "_" + str(ci) +
                                 "]"] = cseq
        ci += 1
    return (True, seqs, common_seqs_used)


def get_lib(path):
    """
    Function which returns a pre-defined library saved to directory.
    Used for sequence back-translation.
    Args:
        path(str): string containing the path to file directory
                   containing mapping.
    Returns:
        dictionary (dict): dictionary containing the {protein : DNA mapping}.
    """

    with open(path) as file:
        keys = []
        vals = []
        for keyvalue in file:
            codon = (keyvalue.strip()).split(":")
            keys.append(codon[0])
            vals.append(codon[1])
        dictionary = dict(zip(keys, vals))
        return dictionary


def back_translate(sequence, library, singleSequence=False):
    """
    Function which maps a protein sequence to a generic DNA sequence.
    Args:
        sequence(str): seqence to be back-transalted.
        library(str): library imported using call to getLib() function.
        singleSequence(Boolean): Boolean indicating if function is applied
    to single sequence or set. It is recommended to apply to a single sequence
    iteratively if input is list or set.
    Returns:
        dna (str): string of converted protein input sequence to generic dna
                   output.
    """
    if type(sequence) != list:
        sequence = list(sequence)

    if not singleSequence:
        frag = sequence[0]
        dna = str()
    else:
        frag = sequence
        dna = str()

    for amino in frag:
        dna += library[amino]

    return dna


def translate(sequence, library, single_sequence=True):
    """
    Function which maps a generic DNA sequence to a protein sequence.
    Opposite of back-translate.
    Uses the getLib function output together with DNA sequence to be translated
    to original protein sequence.

    Args:
        sequence(str): seqence to be transalted.
        library(str): library imported using call to getLib() function.
        singleSequence(Boolean): Boolean indicating if function is applied to
    single sequence or set. It is recommended to apply to a single sequence
    iteratively if input is list or set.

    Returns:
        protein (str): string of converted dna input sequence
    to protein output.
    """
    if type(sequence) == list:
        sequence = "".join(sequence)
    frag = sequence
    protein = str()

    # split the dna sequence into codons
    frag = [frag[x: x + 3] for x in range(0, len(frag), 3)]
    for dna in frag:
        for key, values in library.items():
            if dna in values:
                protein += key
    return protein


def define_sets(seqs, common_seqs_used):
    """
    Function which concatentates each variable protein sequence with
    the constant defined regions from 'getCommonSeqs()' output.

    Args:
        sequence(str): protein sequence with constant regions ID substituted.
        comonSeqsUSed(list): output from getCommonSeqs() function.
    Returns:
        allSets(dict): concatenated dictionary containing both constant and
    variable regions.
    """

    v = []
    for key, value in seqs.items():
        no_brackets = re.sub(r"([\(\[]).*?([\)\]])", ",", seqs[key])
        nb = list(no_brackets.split(","))
        if "" in nb:
            nb.remove("")
        v.append(nb)
    v_seqs = sorted(
        np.unique([val for sublist in v for val in sublist]),
        key=len,
        reverse=True
    )

    # combine into a dictionary:
    i = 0
    variable_seqs = {}
    for variable in v_seqs:
        if variable != "":
            variable_seqs.update({"[" + "v_" + str(i) + "]": variable})
            i += 1

    # append variable regions to common regions:
    all_sets = {**common_seqs_used, **variable_seqs}
    return all_sets


# Functions to reconstruct the recipe back to original sequence.
def strip_word(word):
    x = ["[" + part.strip() + "]"
         for part in re.split(r"[\[\]]", word) if part.strip()]
    return x


def cook(word, all_sets):
    """
    Function which converst a single ID substituted 'recipe' into a fully
    'baked' sequence.
    Args:
        word(str): string containing the recipe for a given sequence.
    This will be split and used as a key to index the corresponding sequence of
    the provided dict(allSets).
        allSets(dict): dictionay containing all of the constant and variable
    regions stored by ID(key), depending on values can be
    allSetsOptimised (DNA)
        or allSets (protein).

    Returns:
        seqs(dict): dictionary of {sequenceName: allSets}.
    """

    stripped = strip_word(word)
    for strip in stripped:
        match = re.search(r"\[(.*?)\]", strip)
        ID = match[1]  # val in brackets
        for key, val in all_sets.items():
            if key == strip:
                word = re.sub(r"\[(%s)\]" % ID, val, word)
    return word


# reconstruct collection of recipes (all encoded strings)
def cookAll(recipe, all_sets):
    """Function which calls cook() to evey recipe in a dictionary.

    Args:

        recipe(dict): sdictionary, values are strings containing the recipe for
    a key indexed sequence. This will be split and used as a key to
    convert the corresponding recipe to the value prodived by allSets.

        allSets(dict): dictionay containing all of the constant and
    variable regions stored by ID(key), depending on values can be
    allSetsOptimised (DNA) or allSets (protein).

    Returns:
        cooked(dict): dictionary of {sequenceName: allSets} depending on values
    of allSets, can be protein or optimizedDNA.
    """

    cooked = {key: cook(val, all_sets) for (key, val) in recipe.items()}
    return cooked


def get_codon_lib(name="resources/codon_table.txt"):
    """Function which provieds a mapping dictionarye.g. {optimizedDNA: protein}
    while considering dna redundancy
    by using a translation table. Should be used as the libary parameter for
    translate() function.

    Args:
        name(str): string path to file containing transaltion table.
    Returns:
        mapping(dict): dictionary containing the mapping
    {amino-acid: list(dna alternative)}.
    """

    # map dna to protein using a codon translation table
    # (using IUPAC generic alphabet)
    with open(name) as file:
        mapping = {}
        vals = []
        for line in file:
            keyval = (line.strip()).split(":")
            key = keyval[0]
            val = keyval[1].split(",")

            for codon in val:
                vals.append(codon)
            # print(vals)

            mapping[key] = vals
            vals = []
        return mapping


def reverse_compliment(s):
    """Function which returns the reverse compliment of an input sequence 's'.

    Args:
        s(str): DNA sequence to be used for reverse compliment proceedure.

    Returns:
        complimentReverseSeq(str): string of the reverse compliment of
    the input.
    """

    # reverse the seq
    reverse_seq = s[::-1]

    # make reverse compliment funct for golden hinge API
    # optimization preventSeq.
    compliment_reverse_seq = ""
    for char in reverse_seq:
        if char == "A":
            compliment_reverse_seq += "T"
        elif char == "T":
            compliment_reverse_seq += "A"
        elif char == "C":
            compliment_reverse_seq += "G"
        elif char == "G":
            compliment_reverse_seq += "C"
    return compliment_reverse_seq


def contact_api(sequence, organism, prevent_seq, x_api_key):
    """Function which connects to ThermoAPI for codon optimization of a
    generic DNA seqeunce.

    Args:
        sequence (str): DNA generic back-translated protein seq
    of common/variable sequences.
        organism(str): Specify the organism of choice e.g
    'H.sapiens' or 'Esirichia.coli'.
        preventSeq(str): String containing sequence not to be
    introduced during call to API
        codon optimization.
        x_api_key(str): API token
    Returns:
        str or obj(uncomment as desired): codon optimized form of DNA
    with specified arguments.
    """

    start_zero_based = 0
    length = len(sequence)

    # params
    optimization_style = "FULLEXPRESSIONOPTIMIZATION"

    # write JSON
    data = {
        "sequence": sequence,
        "organism": organism,
        "orfs": [{"startZeroBased": start_zero_based, "length": length}],
        "optimizationStyle": optimization_style,
        "customMotifs": {
            "preventSeq": prevent_seq,
            "reversePreventSeq": reverse_compliment(prevent_seq),
        },
    }

    accept = "application/json"
    content = "application/json"
    headers = {"x-api-key": x_api_key,
               "accept": accept,
               "Content-Type": content}

    post_url = (
        "https://www.thermofisher.com/order/gene-design-services"
        "/api/optimization/v1/"
    )
    post = requests.post(post_url, data=json.dumps(data), headers=headers)
    print("job status: pending", file=sys.stderr)

    while post.json()["status"] != "done":
        post = requests.post(post_url, data=json.dumps(data), headers=headers)

    print("job status: done", file=sys.stderr)

    get_seq = post.json()["content"]["sequence"]

    return get_seq


def main():
    fasta_file = sys.argv[1]
    threshold = int(sys.argv[2])
    # full_length_threshold = int(sys.argv[3])
    organism = sys.argv[4]
    prevent_seq = sys.argv[5]
    api_token = sys.argv[6]

    seqs = read_fasta(fasta_file)

    # make a copy of seqs
    original = seqs.copy()

    # initialise vars
    cont = True
    passage = 0

    # run substring algorithm if sequences exceed minimum
    # threshold (i.e. > 1kbp)
    # minimal_cases = {k: v for (k, v) in seqs.items()
    #                  if len(v) < full_length_threshold}

    if len(seqs.values()) != 0:
        # Run function to get common seqs
        common_seqs_used = {}
        while cont:
            cont, seqs, common_seqs_used = get_common_seqs(
                seqs, passage, common_seqs_used, minLength=threshold
            )
            passage += 1
    # else:
    #     minimal_seqs = optimise_simple(
    # sequence=minimal_cases, organism=organism)

    lib = get_lib("resources/generic.fas")
    all_sets_back_translated = {}
    all_sets_optimized = {}

    # Define API vars

    # create recipe for any seqs over threshold
    if len(seqs.items()) > 0:
        all_sets = define_sets(seqs, common_seqs_used)

        # back transalte each seq in CommonSeqsUsed to generic DNA
        for key, value in all_sets.items():
            all_sets_back_translated.update(
                {key: back_translate(sequence=value,
                                     library=lib,
                                     singleSequence=True
                                     )}
            )

        print("Displaying all recipe items", file=sys.stderr)
        print(all_sets, file=sys.stderr)
        print("\n", file=sys.stderr)
        print("Converting allSets to optimal DNA:", file=sys.stderr)
        print("\n", file=sys.stderr)

        # codon optimize each of the generic DNA seqs in allSetsBackTranslate
        for key, value in all_sets_back_translated.items():
            all_sets_optimized.update(
                {
                    key: contact_api(
                        sequence=value,
                        organism=organism,
                        prevent_seq=prevent_seq,
                        x_api_key=api_token,
                    )
                }
            )

        # overwite common regions (var and constant) with ID's
        # to produce recipe
        recipe = seqs
        for key, value in recipe.items():
            for key_set, val_set in all_sets.items():
                if val_set in value:
                    value = value.replace(val_set, key_set)
                recipe[key] = value

        print("\n", file=sys.stderr)
        print("Displaying constructed recipe:", file=sys.stderr)
        print(recipe, file=sys.stderr)

        # apply function to map encoding back to original sequences
        cooked = cookAll(recipe, all_sets_optimized)

        # re-translate the optimised (cooked) DNA to protein on a single
        # sequence (test)

        # re-translate the optimised DNA to he original protein sequence using
        # codon table(see library)
        trans_opt_seqs = {}
        for key, value in cooked.items():
            trans_opt_seqs.update(
                {
                    key: translate(
                        sequence=value,
                        library=get_codon_lib(),
                        single_sequence=True
                    )
                }
            )

        # check to see if mapping from recipe -> DNA ->
        # original protein is correct
        print("\n", file=sys.stderr)
        print(
            (
                "Checking if all re-translated values equal "
                "their original sequences:"
            ),
            file=sys.stderr,
        )

        same_as_input = []

        for key, val in trans_opt_seqs.items():
            print(key, trans_opt_seqs[key] == original[key], file=sys.stderr)
            print(key, "optimized seq", trans_opt_seqs[key], file=sys.stderr)
            print(key, "original seq", original[key], file=sys.stderr)
            same_as_input.append(trans_opt_seqs[key] == original[key])

    print("\n", file=sys.stderr)
    print("(If any False, try running whole script again).", file=sys.stderr)

    if all(same_as_input):
        for key, value in cooked.items():
            print(">", key.strip(), sep="")
            print(value.strip(), sep="")


if __name__ == "__main__":
    main()
