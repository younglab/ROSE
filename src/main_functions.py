def checkRefCollection(referenceCollection):

    '''
    makes sure the names of all loci in the reference collection are unique
    '''

    namesList = [locus.ID() for locus in referenceCollection.getLoci()]
    
    if len(namesList) != len(ROSE_utils.uniquify(namesList)):
        print("ERROR: REGIONS HAVE NON-UNIQUE IDENTIFIERS")
        print("THE SECOND COLUMN OF THE INPUT .GFF OR THE FOURTH COLUMN OF THE INPUT .BED MUST HAVE A UNIQUE IDENTIFIER FOR EACH REGION")
        sys.exit()
    else:
        print("REFERENCE COLLECTION PASSES QC")
        return


def regionStitching(inputGFF,stitchWindow,tssWindow,annotFile,removeTSS=True):
    print('PERFORMING REGION STITCHING')
    #first have to turn bound region file into a locus collection

    #need to make sure this names correctly... each region should have a unique name
    boundCollection = ROSE_utils.gffToLocusCollection(inputGFF)

    debugOutput = []
    #filter out all bound regions that overlap the TSS of an ACTIVE GENE
    if removeTSS:
        #first make a locus collection of TSS
        startDict = ROSE_utils.makeStartDict(annotFile)

        #now makeTSS loci for active genes
        removeTicker=0
        #this loop makes a locus centered around +/- tssWindow of transcribed genes
        #then adds it to the list tssLoci
        tssLoci = []
        for geneID in startDict.keys():
            tssLoci.append(ROSE_utils.makeTSSLocus(geneID,startDict,tssWindow,tssWindow))


        #this turns the tssLoci list into a LocusCollection
        #50 is the internal parameter for LocusCollection and doesn't really matter
        tssCollection = ROSE_utils.LocusCollection(tssLoci,50)

        #gives all the loci in boundCollection
        boundLoci = boundCollection.getLoci()

        #this loop will check if each bound region is contained by the TSS exclusion zone
        #this will drop out a lot of the promoter only regions that are tiny
        #typical exclusion window is around 2kb
        for locus in boundLoci:
            if len(tssCollection.getContainers(locus,'both'))>0:
                
                #if true, the bound locus overlaps an active gene
                boundCollection.remove(locus)
                debugOutput.append([locus.__str__(),locus.ID(),'CONTAINED'])
                removeTicker+=1
        print('REMOVED %s LOCI BECAUSE THEY WERE CONTAINED BY A TSS' % (removeTicker))

    #boundCollection is now all enriched region loci that don't overlap an active TSS
    stitchedCollection = boundCollection.stitchCollection(stitchWindow,'both')

    if removeTSS:
        #now replace any stitched region that overlap 2 distinct genes
        #with the original loci that were there
        fixedLoci = []
        tssLoci = []
        for geneID in startDict.keys():
            tssLoci.append(ROSE_utils.makeTSSLocus(geneID,startDict,50,50))


        #this turns the tssLoci list into a LocusCollection
        #50 is the internal parameter for LocusCollection and doesn't really matter
        tssCollection = ROSE_utils.LocusCollection(tssLoci,50)
        removeTicker = 0
        originalTicker = 0
        for stitchedLocus in stitchedCollection.getLoci():
            overlappingTSSLoci = tssCollection.getOverlap(stitchedLocus,'both')
            tssNames = [startDict[tssLocus.ID()]['name'] for tssLocus in overlappingTSSLoci]
            tssNames = ROSE_utils.uniquify(tssNames)
            if len(tssNames) > 2:
            
                #stitchedCollection.remove(stitchedLocus)
                originalLoci = boundCollection.getOverlap(stitchedLocus,'both')
                originalTicker+=len(originalLoci)
                fixedLoci+=originalLoci
                debugOutput.append([stitchedLocus.__str__(),stitchedLocus.ID(),'MULTIPLE_TSS'])
                removeTicker+=1
            else:
                fixedLoci.append(stitchedLocus)

        print('REMOVED %s STITCHED LOCI BECAUSE THEY OVERLAPPED MULTIPLE TSSs' % (removeTicker))
        print('ADDED BACK %s ORIGINAL LOCI' % (originalTicker))
        fixedCollection = ROSE_utils.LocusCollection(fixedLoci,50)
        return fixedCollection,debugOutput
    else:
        return stitchedCollection,debugOutput


def mapCollection(stitchedCollection,referenceCollection,bamFileList,mappedFolder,output,refName):


    '''
    makes a table of factor density in a stitched locus and ranks table by number of loci stitched together
    '''

    
    print('FORMATTING TABLE')
    loci = stitchedCollection.getLoci()

    locusTable = [['REGION_ID','CHROM','START','STOP','NUM_LOCI','CONSTITUENT_SIZE']]

    lociLenList = []

    #strip out any that are in chrY
    for locus in list(loci):
        if locus.chr() == 'chrY':
            loci.remove(locus)
    
    for locus in loci:
        #numLociList.append(int(stitchLocus.ID().split('_')[1]))
        lociLenList.append(locus.len())
        #numOrder = order(numLociList,decreasing=True)
    lenOrder = ROSE_utils.order(lociLenList,decreasing=True)
    ticker = 0
    for i in lenOrder:
        ticker+=1
        if ticker%1000 ==0:
            print(ticker)
        locus = loci[i]

        #First get the size of the enriched regions within the stitched locus
        refEnrichSize = 0
        refOverlappingLoci = referenceCollection.getOverlap(locus,'both')
        for refLocus in refOverlappingLoci:
            refEnrichSize+=refLocus.len()

        try:
            stitchCount = int(locus.ID().split('_')[0])
        except ValueError:
            stitchCount = 1
        
        locusTable.append([locus.ID(),locus.chr(),locus.start(),locus.end(),stitchCount,refEnrichSize])
        
            

    print('GETTING MAPPED DATA')
    for bamFile in bamFileList:
        
        bamFileName = bamFile.split('/')[-1]

        print('GETTING MAPPING DATA FOR  %s' % bamFile)
        #assumes standard convention for naming enriched region gffs
        
        #opening up the mapped GFF
        print('OPENING %s%s_%s_MAPPED.gff' % (mappedFolder,refName,bamFileName))

        mappedGFF =ROSE_utils.parseTable('%s%s_%s_MAPPED.gff' % (mappedFolder,refName,bamFileName),'\t')        

        signalDict = defaultdict(float)
        print('MAKING SIGNAL DICT FOR %s' % (bamFile))
        mappedLoci = []
        for line in mappedGFF[1:]:

            chrom = line[1].split('(')[0]
            start = int(line[1].split(':')[-1].split('-')[0])
            end = int(line[1].split(':')[-1].split('-')[1])
            mappedLoci.append(ROSE_utils.Locus(chrom,start,end,'.',line[0]))
            try:
                signalDict[line[0]] = float(line[2])*(abs(end-start))
            except ValueError:
                print('WARNING NO SIGNAL FOR LINE:')
                print(line)
                continue
                
                
        
        mappedCollection = ROSE_utils.LocusCollection(mappedLoci,500)
        locusTable[0].append(bamFileName)

        for i in range(1,len(locusTable)):
            signal=0.0
            line = locusTable[i]
            lineLocus = ROSE_utils.Locus(line[1],line[2],line[3],'.')
            overlappingRegions = mappedCollection.getOverlap(lineLocus,sense='both')
            for region in overlappingRegions:
                signal+= signalDict[region.ID()]
            locusTable[i].append(signal)

    ROSE_utils.unParseTable(locusTable,output,'\t')
