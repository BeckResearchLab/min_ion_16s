#!/usr/bin/env python3
"""Assembles a set of fasta files and RDP output to build an amplicon data sheet"""

import os

from Bio.SeqIO.FastaIO import SimpleFastaParser
import click
import numpy as np
import pandas as pd

CLASSIFIED_COLUMNS = [
    "id", "classified",
    "rank1name", "rank1", "rank1conf",
    "rank2name", "rank2", "rank2conf",
    "rank3name", "rank3", "rank3conf",
    "rank4name", "rank4", "rank4conf",
    "rank5name", "rank5", "rank5conf",
    "rank6name", "rank6", "rank6conf",
    "rank7name", "rank7", "rank7conf",
    "rank8name", "rank8", "rank8conf",
    "rank9name", "rank9", "rank9conf",
]
TAX = ['domain', 'phylum', 'class','subclass', 'order', 'suborder', 'family', 'genus']

@click.command()
@click.option("--base", "-b", multiple=True, help="base name of .fasta and .classified.txt files to unify, may be multiple")
@click.option("--cache/--no-cache", default=True, help="should tab separated value caches be used if present to speed up processing?")
@click.option("--sequence_join", default=False, help="should the samples be joined on sequences?")
@click.option("--uniq", default=False, help="should the list of unique sequences be written to unique_sequences.txt?")
@click.option("--taxa_join", default=True, help="should the samples be joined on taxaonomy?")
@click.option("--taxa-join-conf-cut", default=.85, type=float)
@click.option("--taxa-statistics", default=True, help="should a statistical analysis be performed of the taxa distributions")
def main(base, cache, sequence_join, uniq, taxa_join, taxa_join_conf_cut, taxa_statistics):
    """Driver for amplicon data sheet tool"""
    df_dict = {}
    for basefile in base:
        print(f"processing for {basefile}")
        fasta = basefile + ".fasta"
        classified = basefile + ".classified.txt"
        if not os.path.exists(fasta):
            raise FileNotFoundError(f"unable to find the file {fasta}")
        if not os.path.exists(classified):
            raise FileNotFoundError(f"unable to find the file {classified}")

        # if there is a cached main table, read it and skip building it
        main = basefile + ".main.tsv"
        if os.path.exists(main) and cache:
            print(f"* found cached main table, reading from {main}")
            sample_df = pd.read_csv(main, sep='\t', index_col=False)
            print(f"+ main table has {sample_df.shape[0]} rows and {sample_df.shape[1]} columns")
            # save for merger later
            df_dict[basefile] = sample_df
            continue

        print(f"* reading {classified}")
        class_tsv = basefile + "_classified.tsv"
        if not os.path.exists(class_tsv) or not cache:
            orig_class_df = pd.read_csv(classified, sep="\t", names=CLASSIFIED_COLUMNS)
            print(f"+ found {orig_class_df.shape[0]} rows and {orig_class_df.shape[1]} columns in original classified.txt")
            class_df = rebuild_classified_df(orig_class_df)
            print(f"+ writing cached version of reprocessed data to {class_tsv}")
            class_df.to_csv(class_tsv, sep='\t', index=False)
        else:
            print(f"+ found cached version of reprocessed classifications, reading from {class_tsv}")
            class_df = pd.read_csv(class_tsv, sep='\t', index_col=False)
        class_df.set_index("id", inplace=True)

        print(f"* reading {fasta}")
        fasta_tsv = fasta + ".tsv"
        if not os.path.exists(fasta_tsv) or not cache:
            ids = []
            seqs = []
            lens = []
            with open(fasta) as f_in:
                for idt, sequence in SimpleFastaParser(f_in):
                    ids.append(idt.split()[0])
                    sequence = "".join(sequence)
                    lens.append(len(sequence))
                    seqs.append(sequence)
            fasta_df = pd.DataFrame({"id": ids, "sequence_length": lens, "sequence": seqs})
            print(f"+ writing cached version {fasta_tsv}")
            fasta_df.to_csv(fasta_tsv, sep='\t', index=False)
        else:
            print(f"+ reading cached version {fasta_tsv}")
            fasta_df = pd.read_csv(fasta_tsv, sep='\t', index_col=False)
        fasta_df.set_index("id", inplace=True)
        print(f"+ derived table has {fasta_df.shape[0]} rows and {fasta_df.shape[1]} columns")

        sample_df = class_df.merge(fasta_df, how="outer", on="id", validate="1:1")
        # count the occurences of unique sequences
        print(f"+ found {len(sample_df['sequence'].unique())} unique sequences")
        # setup groupby like table to be joined back
        sample_counts_df = sample_df.value_counts("sequence").reset_index(name="count")
        # join back, note this will be interesting if the sequences are not unique
        sample_df.merge(sample_counts_df, how="inner", on="sequence", validate="1:1")
        print(f"> merged table has {sample_df.shape[0]} rows and {sample_df.shape[1]} columns")
        print(f"+ writing combined file to {main}")
        sample_df.reset_index(inplace=True)
        sample_df.to_csv(main, sep='\t', index=False)

        # save for merger later
        df_dict[basefile] = sample_df

    if len(df_dict) < 2:
        print("only one base file supplied; nothing more to do")
        return

    if sequence_join:
        sequence_join_dfs(df_dict, uniq)

    if taxa_join:
        taxa_join_dfs(df_dict, taxa_join_conf_cut, taxa_statistics)

    return


def taxa_join_dfs(df_dict: dict, taxa_join_conf_cut: float, taxa_statistics: bool):

    if taxa_statistics:
        #print(f"doing a statistical analysis of taxonomy distributions across samples")
        pass

    print(f"joining all samples on taxonomy")
    df_main = None
    last_key = None
    for df_key in df_dict:
        df = df_dict[df_key]
        print(f". merging table {df_key}")
        if df_main:
            # suffixes will not work if joins is odd or even?  last join will leave unsuffixed columns, dacb thought
            df_main = df_main.merge(df, how="outer", on=TAX, suffixes=["_"+last_key, "_"+df_key])
        else:
            df_main = df
        print(df_main.columns)
        last_key = df_key
    df_main.reset_index(inplace=True)
    print(f"+ final taxa joined table has {df_main.shape[0]} rows and {df_main.shape[1]} columns")
    print(f"> writing main taxa joined table to main.taxa.tsv")
    df_main.to_csv("main.taxa.tsv", sep='\t', index=False)

    return


def sequence_join_dfs(df_dict: dict, uniq: bool):
    print(f"finding unique sequences across all samples")
    seqs = pd.Series(dtype="object")
    seq_count = 0
    for df_key in df_dict:
        df = df_dict[df_key]
        seq_count = seq_count + df.shape[0]
        seqs = pd.concat([seqs, df["sequence"]], ignore_index=True)
    print(f"+ total sequences  {seq_count}")
    seqs.sort_values(inplace=True, ignore_index=True)
    print(f"+ unique sequences {len(seqs.unique())}")
    if uniq:
        print(f"+ writing unique sequences to unique_sequences.txt")
        seqs.to_csv("unique_sequences.txt", index=False)

    print(f"joining all samples on sequence")
    df_main = pd.DataFrame(seqs, columns=["sequence"])
    last_key=""
    for df_key in df_dict:
        df = df_dict[df_key]
        print(f". merging table {df_key}")
        # suffixes will not work if joins is odd or even?  last join will leave unsuffixed columns, dacb thought
        df_main = df_main.merge(df, how="outer", on="sequence", suffixes=["_"+last_key, "_"+df_key])
        print(df_main.columns)
        last_key = df_key
    df_main.reset_index(inplace=True)
    print(f"+ final sequence joined table has {df_main.shape[0]} rows and {df_main.shape[1]} columns")
    print(f"> writing main sequence joined table to main.sequences.tsv")
    df_main.to_csv("main.sequences.tsv", sep='\t', index=False)

    return


def rebuild_classified_df(class_df_in):
    print(f"+ reprocessing the classified dataframe for conformity")
    class_df_in.reset_index(inplace=True)
    new_df_rows = class_df_in.shape[0]

    columns = ["id", "classified"]
    for tax in TAX:
        columns.append(tax)
        columns.append(tax+"_conf")

    class_df_dict = {}
    for column in columns:
        class_df_dict[column] = np.empty(new_df_rows, dtype="object")

    for index, row in class_df_in.iterrows():
        class_df_dict["id"][index] = row["id"]
        class_df_dict["classified"][index] = row["classified"]
        # build up the tax dict from the row data
        tax_dict = {}
        for rank in range(1, 10):
            ranks = str(rank)
            tax_dict[row["rank"+ranks]] = (row["rank"+ranks+"name"], row["rank"+ranks+"conf"])
        # assign the tax from the dict
        for tax in TAX:
            name, conf = tax_dict.get(tax, (None, None))
            class_df_dict[tax][index] = name
            class_df_dict[tax+"_conf"][index] = conf
    # create the df from the dictionary of numpy arrays
    class_df = pd.DataFrame(class_df_dict)
    print(f"+ after reprocessing class dataframe had {class_df.shape[0]} rows and {class_df.shape[1]} columns")
    return class_df


if __name__ == "__main__":
    main()
    print("goodbye")
