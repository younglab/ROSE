#!/usr/bin/env python

'''
PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS,
AND RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS
APRIL 11, 2013
VERSION 0.1
CONTACT: youngcomputation@wi.mit.edu
'''

import argparse
# import ROSE_utils
import shutil

from pathlib import Path
from src.utils.annotation import makeStartDict
from src.utils.conversion import bed_to_gff3, gtf_to_gff3
from src.utils.file_helper import get_path, check_file, check_path
from typing import Any, Dict

def str2bool(
    v: str
) -> bool:
    """Convert string to boolean

    Args:
        v (str): boolean string

    Raises:
        argparse.ArgumentTypeError: String is not named "true" or "false"

    Returns:
        bool: Booleanised string
    """
    
    if v.lower() == "true":
        return True
    elif v.lower() == "false":
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def parseArgs() -> argparse.Namespace:
    """Parse arguments from CLI

    Returns:
        argparse.Namespace: Argparse space containing parsed arguments
    """

    parser = argparse.ArgumentParser(description='Stitch regions together to form enhancers, map read density to stitched regions and \
                                                  rank enhancers by read desnity to discover super-enhancers')

    #Required arguments
    parser.add_argument('-g', '--genome', type=str, help='Genome build (MM10, MM9, MM8, HG18, HG19, HG38)')
    parser.add_argument('-i', '--input', type=str, help='File (.bed, .gff or .gtf) containing binding sites to make enhancers')
    parser.add_argument('-o', '--output', type=str, help='Output directory name')
    parser.add_argument('-r', '--rankby',  type=str, help='File (.bam) to rank enhancer by')

    #Optional arguments
    parser.add_argument('-b', '--bams', nargs='*', help="Comma separated list of additional files (.bam) to map to")
    parser.add_argument('-c', '--control',  type=str, nargs='?', help="File (.bam) to rank enhancer by")
    parser.add_argument('-s', '--stitch', type=int, nargs='?', default=12500, help="Max linking distance for stitching")
    parser.add_argument('-t', '--tss', type=int, nargs='?', default=0, help="Distance from TSS to exclude (0 = no TSS exclusion)")
    parser.add_argument('-v', '--verbose', type=str2bool, nargs='?', const=True, default=False, help='Print verbose messages')


    #Printing arguments to the command line
    args = parser.parse_args()

    print("Called with args:")
    print(f"{args}\n")

    #Ensuring that argument files exist
    check_file(args.input)
    check_file(args.rankby)
    check_file(args.control)

    return args


def main():
    '''
    main run call
    '''

    #Parse arguments from the command line
    args = parseArgs()
    path = get_path()

    #Initialising variables
    # annotFile = genomeDict[upper(args.genome)]
    debug = False
    genomeDict = {
        "HG18": Path(path, "data", "annotation", "hg18_refseq.ucsc"),
        "HG19": Path(path, "data", "annotation", "hg19_refseq.ucsc"),
	    "HG38": Path(path, "data", "annotation", "hg38_refseq.ucsc"),
        "MM8":  Path(path, "data", "annotation", "mm8_refseq.ucsc"),
        "MM9":  Path(path, "data", "annotation", "mm9_refseq.ucsc"),
        "MM10": Path(path, "data", "annotation", "mm10_refseq.ucsc"),
        }
    stitchWindow = int(args.stitch)
    tssWindow = int(args.tss)

    # if tssWindow != 0:
    #     removeTSS = True
    # else:
    #     removeTSS = False
    if args.control:        
        bamFileList = [args.rankby, args.control]
    else:
        bamFileList = [args.rankby]
    # if options.bams:
    #     bamFileList += options.bams.split(',')
    #     bamFileLIst = ROSE_utils.uniquify(bamFileList)

    #Ensuring necessary output directories exist
    output = check_path(Path(path, args.output))
    gffFolder = check_path(Path(path, args.output, "gff"))
    mappedFolder = check_path(Path(path, args.output, "mappedGFF"))

    #Copying/creating the input .gff3 file
    if Path(args.input).suffix == ".bed":
        if args.verbose:
            print("Converting input .bed file to .gff3 format")
        inputGFFFile = str(Path(path, "output", "gff", Path(args.input).stem)) + ".gff3"
        # bed_to_gff3(args.input, inputGFFFile)

    elif Path(args.input).suffix == ".gtf":
        if args.verbose:
            print("Converting input .gtf file to .gff3 format")
        inputGFFFile = str(Path(path, "output", "gff", Path(args.input).stem)) + ".gff3"
        # gtf_to_gff3(args.input, inputGFFFile, full=False)

    elif Path(args.input).suffix == ".gff" or Path(args.input).suffix == ".gff3":
        if args.verbose:
            print("Copying input .gff file to new directory")
        inputGFFFile = str(Path(path, "output", "gff", Path(args.input).stem)) + ".gff3"
        # shutil.copyfile(args.input, inputGFFFile)
        
    else:
        raise ValueError("Input file must be a .bed, .gtf, .gff or gff3 file")

    #Setting the bound region file to define enhancers
    if args.verbose:
        print(f"Using {inputGFFFile} as the input .gff file\n")
    inputName = str(Path(inputGFFFile).stem)

    # annotFile = genomeDict[upper(genome)]

    #Making the start dict
    if args.verbose:
        print("Making the start dict")
    startDict = makeStartDict(check_file(str(genomeDict[args.genome.upper()])))


    # #LOADING IN THE BOUND REGION REFERENCE COLLECTION
    # print('LOADING IN GFF REGIONS')
    # referenceCollection = ROSE_utils.gffToLocusCollection(inputGFFFile)

    # #CHECKING INPUT REGIONS FOR FORMATTING
    # print('CHECKING INPUT TO MAKE SURE EACH REGION HAS A UNIQUE IDENTIFIER')
    # checkRefCollection(referenceCollection) #makes sure that all input regions have a unique ID

    # #NOW STITCH REGIONS
    # print('STITCHING REGIONS TOGETHER')
    # stitchedCollection,debugOutput = regionStitching(inputGFFFile,stitchWindow,tssWindow,annotFile,removeTSS)

    
    # #NOW MAKE A STITCHED COLLECTION GFF
    # print('MAKING GFF FROM STITCHED COLLECTION')
    # stitchedGFF=ROSE_utils.locusCollectionToGFF(stitchedCollection)
    
    # if not removeTSS:
    #     stitchedGFFFile = '%s%s_%sKB_STITCHED.gff' % (gffFolder,inputName,stitchWindow/1000)
    #     stitchedGFFName = '%s_%sKB_STITCHED' % (inputName,stitchWindow/1000)
    #     debugOutFile = '%s%s_%sKB_STITCHED.debug' % (gffFolder,inputName,stitchWindow/1000)
    # else:
    #     stitchedGFFFile = '%s%s_%sKB_STITCHED_TSS_DISTAL.gff' % (gffFolder,inputName,stitchWindow/1000)
    #     stitchedGFFName = '%s_%sKB_STITCHED_TSS_DISTAL' % (inputName,stitchWindow/1000)
    #     debugOutFile = '%s%s_%sKB_STITCHED_TSS_DISTAL.debug' % (gffFolder,inputName,stitchWindow/1000)

    # #WRITING DEBUG OUTPUT TO DISK
        
    # if debug:
    #     print('WRITING DEBUG OUTPUT TO DISK AS %s' % (debugOutFile))
    #     ROSE_utils.unParseTable(debugOutput,debugOutFile,'\t')

    # #WRITE THE GFF TO DISK
    # print('WRITING STITCHED GFF TO DISK AS %s' % (stitchedGFFFile))
    # ROSE_utils.unParseTable(stitchedGFF,stitchedGFFFile,'\t')



    # #SETTING UP THE OVERALL OUTPUT FILE
    # outputFile1 = outFolder + stitchedGFFName + '_ENHANCER_REGION_MAP.txt'

    # print('OUTPUT WILL BE WRITTEN TO  %s' % (outputFile1))
    
    # #MAPPING TO THE NON STITCHED (ORIGINAL GFF)
    # #MAPPING TO THE STITCHED GFF


    # # bin for bam mapping
    # nBin =1

    # #IMPORTANT
    # #CHANGE cmd1 and cmd2 TO PARALLELIZE OUTPUT FOR BATCH SUBMISSION
    # #e.g. if using LSF cmd1 = "bsub python bamToGFF.py -f 1 -e 200 -r -m %s -b %s -i %s -o %s" % (nBin,bamFile,stitchedGFFFile,mappedOut1)

    # for bamFile in bamFileList:

    #     bamFileName = bamFile.split('/')[-1]

    #     #MAPPING TO THE STITCHED GFF
    #     mappedOut1 ='%s%s_%s_MAPPED.gff' % (mappedFolder,stitchedGFFName,bamFileName)
    #     #WILL TRY TO RUN AS A BACKGROUND PROCESS. BATCH SUBMIT THIS LINE TO IMPROVE SPEED
    #     cmd1 = "python ROSE_bamToGFF.py -f 1 -e 200 -r -m %s -b %s -i %s -o %s &" % (nBin,bamFile,stitchedGFFFile,mappedOut1)
    #     print(cmd1)
    #     os.system(cmd1)

    #     #MAPPING TO THE ORIGINAL GFF
    #     mappedOut2 ='%s%s_%s_MAPPED.gff' % (mappedFolder,inputName,bamFileName)
    #     #WILL TRY TO RUN AS A BACKGROUND PROCESS. BATCH SUBMIT THIS LINE TO IMPROVE SPEED
    #     cmd2 = "python ROSE_bamToGFF.py -f 1 -e 200 -r -m %s -b %s -i %s -o %s &" % (nBin,bamFile,inputGFFFile,mappedOut2)
    #     print(cmd2)
    #     os.system(cmd2)
        

    
    # print('PAUSING TO MAP')
    # time.sleep(10)

    # #CHECK FOR MAPPING OUTPUT
    # outputDone = False
    # ticker = 0
    # print('WAITING FOR MAPPING TO COMPLETE. ELAPSED TIME (MIN):')
    # while not outputDone:

    #     '''
    #     check every 5 minutes for completed output
    #     '''
    #     outputDone = True
    #     if ticker%6 == 0:
    #         print(ticker*5)
    #     ticker +=1
    #     #CHANGE THIS PARAMETER TO ALLOW MORE TIME TO MAP
    #     if ticker == 144:
    #         print('ERROR: OPERATION TIME OUT. MAPPING OUTPUT NOT DETECTED')
    #         exit()
    #         break
    #     for bamFile in bamFileList:
            
    #         #GET THE MAPPED OUTPUT NAMES HERE FROM MAPPING OF EACH BAMFILE
    #         bamFileName = bamFile.split('/')[-1]
    #         mappedOut1 ='%s%s_%s_MAPPED.gff' % (mappedFolder,stitchedGFFName,bamFileName)

    #         try:
    #              mapFile = open(mappedOut1,'r')
    #              mapFile.close()
    #         except IOError:
    #             outputDone = False

    #         mappedOut2 ='%s%s_%s_MAPPED.gff' % (mappedFolder,inputName,bamFileName)
            
    #         try:
    #             mapFile = open(mappedOut2,'r')
    #             mapFile.close()
    #         except IOError:
    #             outputDone = False
    #     if outputDone == True:
    #         break
    #     time.sleep(300)
    # print('MAPPING TOOK %s MINUTES' % (ticker*5))

    # print('BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS')
    # #CALCULATE DENSITY BY REGION
    # mapCollection(stitchedCollection,referenceCollection,bamFileList,mappedFolder,outputFile1,refName = stitchedGFFName)


    # time.sleep(10)

    # print('CALLING AND PLOTTING SUPER-ENHANCERS')


    # if options.control:

    #     rankbyName = options.rankby.split('/')[-1]
    #     controlName = options.control.split('/')[-1]
    #     cmd = 'R --no-save %s %s %s %s < ROSE_callSuper.R' % (outFolder,outputFile1,inputName,controlName)

    # else:
    #     rankbyName = options.rankby.split('/')[-1]
    #     controlName = 'NONE'
    #     cmd = 'R --no-save %s %s %s %s < ROSE_callSuper.R' % (outFolder,outputFile1,inputName,controlName)
    # print(cmd)
    # os.system(cmd)



if __name__ == "__main__":
    main()