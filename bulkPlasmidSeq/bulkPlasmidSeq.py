import argparse
import os
import sys
import shutil
import subprocess

def main():
    '''
    See dependencies in accompanying text file.
    
    BulkPlasmidSeq usage examples (Not updated with Porechop implementation):
    
    For filtered sequences with known reference sequence:
    
        python bulkPlasmidSeq.py -i my_reads.fastq -r my_plasmid_genome.fa -o output_directory -t 4
  
        python bulkPlasmidSeq.py -i path/to/reads -r path/to/plasmids
        
    For unfiltered reads - Guppy filters to Q7 ~85% basecalling accuracy:
    
        python bulkPlasmidSeq.py -i my_reads.fastq -r my_plasmid_genome.fa -q 8 
        
        
        
    '''
    args = getArgs()
    
    reads, reference, outputDir = loadReads(args.input_reads, args.reference, args.output_dir)
    
    if args.porechop and reference is not None:
        
        sys.exit('Porechop has no option "reference". If you would like to bin with Porechop'
                 'and also polish with medaka, please only specify -r or -BC'
                'The same set of seqs will be used for both.')
    
    if args.porechop:
        runPorechop(reads, outputDir, args.barcodes, args.barcode_threshold,
                              args.threads, args.end_size, args.medaka, args.porechop_iterations)
    
    if args.medaka and args.porechop == False:
        runMedaka(reads, reference, outputDir, args.threads, args.porechop)
    

def getArgs():
    '''
    Parses arguments:
        - See code for available options
        - (i, r, o, t) are the standard arguments needed to process reads which have been filtered and have known reference
    
    '''
    
    ap = argparse.ArgumentParser()
    
    mainArgs = ap.add_argument_group('Main options')
    mainArgs.add_argument('-i', '--input_reads', required=True,
                          help= 'input reads (directory to reads or .fastq)')
    mainArgs.add_argument('-r', '--reference', required=False,
                          help = 'plasmid sequences (directory or .fa)')
    mainArgs.add_argument('-o', '--output_dir', required = False,
                          help = 'output directory for medaka')
    mainArgs.add_argument('-t', '--threads', required = False, default = 2,
                          help = 'number of threads for medaka, default 2')
    mainArgs.add_argument('-M', '--medaka', action = 'store_true', required = False,
                          help = 'Run medaka')
    mainArgs.add_argument('-P', '--porechop', action = 'store_true', required = False,
                          help = 'Bin reads using porechop')
    
    porechopArgs = ap.add_argument_group('Porechop Options')
    porechopArgs.add_argument('--barcode_threshold', required = False, default = 75)
    porechopArgs.add_argument('--end_size', required = False, default = 1000)
    porechopArgs.add_argument('-BC', '--barcodes', required = False, help = 'Make barcodes')
    porechopArgs.add_argument('-N', '--porechop_iterations', required = False, default = 1,
                              help = 'Number of rounds of Porechop binning')
    

    args = ap.parse_args()
                
    if args.medaka and args.porechop:
        
        print('Binning the reads with porechop first, then assembling consensus using medaka')
        
    if args.medaka == False and args.porechop == False:
        
        sys.exit('Need to specify -M to generate consensus without binning (Medaka),'
        'or -P to bin reads on barcodes (Porechop)')
       
    
   
    return args

def loadReads(inputFiles, referenceFiles, outputDir):
    '''
    You're gonna regret putting outoutDir here when you can't find it later........
    '''
    if os.path.isdir(inputFiles):
        
        print('Directory was given for input reads, concatenating these to output_reads.fastq')
        
        subprocess.run(['cat %s/*.fastq > output_reads.fastq' % inputFiles], shell = True)
        reads = 'output_reads.fastq'
                                          
    if os.path.isfile(inputFiles):
        #Checks if fastq format
        if str(inputFiles).lower().endswith('.fastq'):
            reads = inputFiles
        else:
            sys.exit('Error: Reads should be in .fastq format ')
    
    if referenceFiles is not None:
        if os.path.isdir(referenceFiles):

            print('Directory was given for plasmid reference, concatenating these to output_ref.fasta')

            subprocess.run(['cat %s/*.fa* > output_ref.fasta' % referenceFiles], shell = True)
            reference = 'output_ref.fasta'


        if os.path.isfile(referenceFiles):
            #Checks if reference is already assembled in to plasmid genome
            if str(referenceFiles).lower().endswith('.fa') or str(referenceFiles).lower().endswith('.fasta'):
                reference = referenceFiles
            else:
                sys.exit('Error: Reference sequence should be in .fasta or .fa format')
                
    else:
        reference = referenceFiles
            
    return reads, reference, outputDir

def runMedaka(reads, reference, outputDir, threads, porechop = False):
    #Make this more modular, here run medaka
    print('----------------------------------\n')
    print('Checking Medaka Args')
    print('----------------------------------\n')
    print('Input files: %s \nReference Files: %s \nOutput directory: %s \nThreads: %s \nPorechop: %s'
         % (reads, reference, outputDir, str(threads), str(porechop)))
    
    subprocess.run(['medaka_consensus -i %s -d %s -o %s -t %s -m r941_min_fast_g344'
                    % (reads, reference, outputDir, str(threads))], shell = True)
    
    processMedakaOutput(outputDir, porechop)
    
    return 


def filterReads():
    return

    
def processMedakaOutput(outputDir, porechop):
    '''
    Objective:
    Splits the consensus.fasta file into individual plasmid seqs. 
    
    '''
    if os.path.isdir(outputDir):
        filePath = outputDir+'/consensus.fasta'
        consensusFasta = file_object = open(filePath, 'rt')
        for line in consensusFasta:
            if line.startswith('>'):
                header = line.split(' ')[0]  
            else:
                seq = line
                out = open('%s.fasta' % header[1:-1], 'w')
                out.write('%s %s' % (header, seq))

    
    elif porechop == True:
        #This should automatically be set by -o when running medaka, if it's not a dir, idk what to do with that
        print('Porechop barcode binning was done prior to running medaka')
        
        
def runPorechop(reads, outputDir, barcodes, barcodeThreshold, threads, endSize, medaka, iterations):
    #print('Pretending to run Porechop here, processing existing directories')
    if not os.path.isdir(outputDir):
            os.makedirs(outputDir)
    detailsOut = open(os.path.join(outputDir, 'PorechopRunDetails.txt'), 'wt')
    
    detailsOut.write('Porechop run details \nInput files: %s \nBarcodes: %s \nOutput Directory:'
                     '%s \nBarcode Threshold: %s \nThreads: %s \nEnd Size: %s \nMedaka: %s \nInterations %s\n'
         % (reads, barcodes, outputDir, barcodeThreshold, str(threads), str(endSize), str(medaka), str(iterations)))
    
    detailsOut.close()
    
    for x in range(int(iterations)):
        
        baseRun = 'porechop -i %s -b %s -BC %s -t %s --barcode_threshold %s \
        --no_split --untrimmed -v 0 --end_size %s --discard_unassigned'
        
        if x == 0:  
            subprocess.run([baseRun % (reads, outputDir_0, barcodes, str(threads),
                                       str(barcodeThreshold), str(endSize))], shell = True)
            
        else:
            subprocess.run([baseRun % (reads, outputDir+'_'+str(x), barcodes, str(threads),
                                       str(barcodeThreshold), str(endSize))], shell = True)
    
    if medaka:
        
        runMedaka(reads, reference, barcodes, threads, porechop = True)
        
    processPorechopOutput(outputDir, iterations, barcodes)
    
    return

def processPorechopOutput(outputDir, iterations, barcodes):
    
    paths = [path for path in os.listdir('.') if path.startswith(outputDir)]
    
    if len(paths) != int(iterations):
        
        sys.exit('Number of directories is not the same as number of iterations, one or more'
                'rounds of Porechop likely failed, try again with new parameters (--end_size) recommended')
    
    else:
        
        #barcodes = [BC for BC in os.listdir(paths[0]) if BC.startswith('BC')]
        
        finalOutput = outputDir + '_Final'
        
        os.mkdir(finalOutput)
        
        finalFastq = open('allBinned.fastq', 'wt')
        
        
        for directory in paths:
            
            subprocess.run(['cat %s/*.fastq >> %s/allBinned.fastq' % (directory, finalOutput)], shell = True)
            #subprocess.run(['rm -r %s' % directory], shell = True)
        
        finalFastq.close()
        
        runMinimap2(outputDir, barcodes)
        
    return 

        
def runMinimap2(outputDir, reference):
    
    finalReads = outputDir + '_Final/allBinned.fastq'
    
    minimap2Run = 'minimap2 -ax map-ont %s %s | samtools sort -o %s'
    
    #samtoolsRun = 'samtools sort | samtools view -q > %s/filtered_sorted_binned_reads.sam'
    
    subprocess.run([minimap2Run % (reference, finalReads, outputDir+'_Final/sorted_reads.sam')], shell = True)
        
    return 
       
if __name__ == "__main__":

    main()