#!/usr/bin/env python

import pandas as pd

from typing import Dict, List, Tuple


# def importRefseq(
#     refseqFile: str
# ) -> Tuple[pd.core.frame.DataFrame, Dict[str, List[int]]]:
#     """Open UCSC refseq file and create dictionary of gene entries

#     Args:
#         refseqFile (str): Path to the UCSC refseq file

#     Returns:
#         Tuple[pd.core.frame.DataFrame, Dict[str, List[int]]]: Tuple of UCSC refseq dataframe and dictionary gene entries
#     """

#     #Initialising variables
#     refseqDict = {}

#     #Reading the UCSC refseq file in as a dataframe
#     refseqTable = pd.read_csv(refseqFile, sep="\t")
#     for key, val in refseqTable["name"].to_dict().items():
#         if val in refseqDict:
#             refseqDict[val].append(key)
#         else:
#             refseqDict[val] = [key]

#     return refseqTable, refseqDict


def makeStartDict(
    annotFile: str
) -> pd.core.frame.DataFrame:
    """Create subsetted dataframe from UCSC refseq file

    Args:
        annotFile (str): Path to UCSC refseq annotation file

    Returns:
        pd.core.frame.DataFrame: Subsetted annotation dataframe
    """

    #Reading the annotation file in as a dataframe
    refseqTable = pd.read_csv(annotFile, sep="\t")

    #Remove duplicate "name" entries and extract specific data
    startDict = refseqTable[refseqTable.columns.intersection(["name", "strand", "chrom", "txStart", "txEnd", "name2"])].copy()
    startDict.drop_duplicates(subset=["name"], inplace=True)
    startDict.loc[startDict["strand"]=="-", ["txStart", "txEnd"]] = (startDict.loc[startDict["strand"] == "-", ("txEnd", "txStart")].values)
    startDict.rename({"name": "id", "strand": "sense", "chrom": "chr", "txStart": "start", "txEnd": "end", "name2": "name"}, axis=1, inplace=True)

    return startDict