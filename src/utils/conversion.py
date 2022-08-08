#!/usr/bin/env python

import pandas as pd
import re

def search(
    attr: str,
    entry: str
) -> str:
    """Checks if an attribute is present for a given .gtf entry

    Args:
        attr (str): .gtf attribute
        entry (str): .gtf attributes entry

    Returns:
        str: Attribute=value concatenation
    """

    #Check if attribute is present for given entry, and if it is, return it's value
    if re.search(f'{attr} "', entry):
        return "{}={}".format(attr, re.findall(rf"(?<={attr} \")[^\"]+", entry)[0])


def bed_to_gff3(
    input: str,
    output: str
) -> None:
    """Convert .bed file to .gff3 file (modeled after R's rtracklayer)

    Args:
        input (str): Input .bed file path
        output (str): Output .gff3 file path
    """

    #Reading the .bed file as a dataframe
    bed_df = pd.read_csv(input, sep="\t", header=None, comment='#')

    #Converting the .bed dataframe to a .gff3 dataframe
    gff_df = pd.DataFrame({
        "seqid": bed_df.iloc[:, 0],
        "source": ["."]*len(bed_df),
        "type": ["sequence_feature"]*len(bed_df),
        "start": bed_df.iloc[:, 1]+1,
        "end": bed_df.iloc[:, 2],
        "score": bed_df.iloc[:, 4],
        "strand": bed_df.iloc[:, 5],
        "phase": ["."]*len(bed_df),
        "attributes": "name="+bed_df.iloc[:, 3]
    })

    #Outputting the gff3 dataframe
    with open(output, "w") as f_out:
        f_out.write("##gff-version 3\n")
        f_out.write("##source-version ROSE\n")
        gff_df.to_csv(f_out, sep="\t", header=False, index=False, mode="a")


def gtf_to_gff3(
    input: str,
    output: str,
    full: bool = False
) -> None:
    """Convert .gtf file to .gff3 file (modeled after R's rtracklayer)

    Args:
        input (str): Input .bed file path
        output (str): Output .gff3 file path
        full (bool, optional): Boolean to fully convert the .gtf atttributes to .gff3 format. Defaults to False.
    """

    #Initialising variables
    gtf_attributes = ["gene_id", "db_xref", "gbkey", "gene", "gene_biotype", "transcript_id",
                      "model_evidence", "product", "exon_number", "protein_id", "anticodon",
                      "inference", "note", "exception", "transl_except", "pseudo", "partial"]

    #Reading the .gtf file as a dataframe
    df = pd.read_csv(input, sep="\t", header=None, comment='#')
    if full:
        df.iloc[:, 8] = [";".join([a for a in [search(attr, annot) for attr in gtf_attributes] if a is not None]) for annot in df.iloc[:, 8].values]
    else:
        df.iloc[:, 8] = [";".join([a for a in [search(attr, annot) for attr in ["gene_id", "transcript_id"]] if a is not None]) for annot in df.iloc[:, 8].values]

    #Outputting the gff3 dataframe
    with open(output, "w") as f_out:
        f_out.write("##gff-version 3\n")
        f_out.write("##source-version ROSE\n")
        df.to_csv(f_out, sep="\t", header=False, index=False, mode="a")
