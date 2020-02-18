'''
Camille Mumm - Boyle Lab rotation 2020

Uses Medaka (https://github.com/nanoporetech/medaka) and Porechop (https://github.com/rrwick/Porechop)
to process reads from bulk plasmid sequencing. 
'''
import argparse
import os
import sys
import shutil
import subprocess
import time

from collections import defaultdict

def main():
    '''
    See dependencies in accompanying text file.
    
    BulkPlasmidSeq usage examples (Not updated with Porechop implementation):
    
    For alignment, polishing, and creating consensus without binning (Medaka):
    
        python bulkPlasmidSeq.py -i 'path/to/reads' -r 'path/to/plasmids' -o outputDirectory  -M
    
    For binning and alignment (Porechop):
    
        python bulkPlasmidSeq.py -i 'path/to/reads' -o outputDirectory -BC 'path/to/plasmids' -P
        
        python bulkPlasmidSeq.py -i 'path/to/reads' -o outputDirectory -BC 'path/to/plasmids' -P --end_size 1000 -N 3        
    '''
    args = getArgs()
    
    if None not in (args.input_reads, args.output_dir):
    
        reads, reference, outputDir = loadReads(args.input_reads, args.reference, args.output_dir)
        
    else:
        sys.exit('Missing one of the following arguments, --input(-i), --reference(-r), --output_dir(-o)')
    
    if args.nanofilt:
        '''
        Runs Nanofilt with specified filters, updates reads to use for downstream analysis
        If none of the defaults are changed, tell the user and tell them to change min_length to something different
        To ignore this check
        '''
        if args.max_length == 100000 and args.min_length == 0 and args.min_quality == 7:
            sys.exit('No arguments given for nanofilt, will not filter out anything from your reads unless'
                     'you are using failed reads from Guppy. If using failed reads use --min_length 1')
            
        reads = filterReads(reads, args.max_length, args.min_length, args.min_quality, args.porechop, args.medaka)
    
    if args.porechop:
        #Check the inputs
        '''
        Checks that all the necessary inputs are provided. Note that -BC is used for building custom barcodes
        and also later as a reference for minimap. This program checks that -r isn't used even though it would do the
        exact same thing.
        '''
        if None not in (outputDir, args.barcodes):
            runPorechop(reads, outputDir, args.barcodes, args.barcode_threshold,
                        args.threads, args.end_size, args.porechop_iterations, args.screenshot, args.igv)
            
        elif reference is not None:
            sys.exit('Porechop has no option "reference". Please use -BC which will be used for barcoding and reference')
            
        else:
            sys.exit('Porechop needs input reads (-i), output directory (-o), and barcodes (-BC)')
    
    if args.medaka and args.porechop == False:
        '''
        Running Medaka, check that Porechop is false. Check the necessary inputs. Checks for use of -r vs -BC. 
        These two args take the same thing! A file or directory full of plasmid, the distiction is used for user clarity.
        '''
        
        if None not in (reads, reference, outputDir):
            runMedaka(reads, reference, outputDir, args.threads, args.porechop, args.screenshot, args.igv)
            
        elif args.barcodes is not None:
            sys.exit('If you are looking for consensus building without binning (Medaka) use -r for reference, not -BC')
        
        else:
            sys.exit('Porechop needs input reads (-i), output directory (-o), and reference sequences (-r)')
    
    
def getArgs():
    '''
    Parses arguments:
        - See -h, --help
    
    '''
    ap = argparse.ArgumentParser()
    
    mainArgs = ap.add_argument_group('Main options')
    mainArgs.add_argument('-i', '--input_reads', required=True,
                          help= 'input reads (directory to reads or .fastq)')
    mainArgs.add_argument('-r', '--reference', required=False,
                          help = 'plasmid sequences (directory or .fa)')
    mainArgs.add_argument('-o', '--output_dir', required = False,
                          help = 'output directory')
    mainArgs.add_argument('-t', '--threads', required = False, default = 2,
                          help = 'number of threads, default 2')
    mainArgs.add_argument('-M', '--medaka', action = 'store_true', required = False,
                          help = 'Run medaka')
    mainArgs.add_argument('-P', '--porechop', action = 'store_true', required = False,
                          help = 'Bin reads using porechop')
    mainArgs.add_argument('-F', '--nanofilt', action = 'store_true', required = False, 
                         help = 'Filter reads with NanoFilt')
    
    porechopArgs = ap.add_argument_group('Porechop Options')
    porechopArgs.add_argument('--barcode_threshold', required = False, default = 75)
    porechopArgs.add_argument('--end_size', required = False, default = 1000)
    porechopArgs.add_argument('-BC', '--barcodes', required = False, help = 'Make barcodes')
    porechopArgs.add_argument('-N', '--porechop_iterations', required = False, default = 1,
                              help = 'Number of rounds of Porechop binning')
    
    nanofilt = ap.add_argument_group('Nanofilt arguments')
    nanofilt.add_argument('--max_length', required = False, default = 100000, help = 'Filtering reads by maximum length')
    nanofilt.add_argument('--min_length', required = False, default = 0,  help = 'Filtering reads by minimum length')
    nanofilt.add_argument('-q', '--min_quality', required = False, default = 7, help = 'Filter reads by quality score > N')
    
    igv = ap.add_argument_group('IGV arguments for screenshotting')
    igv.add_argument('--screenshot', required = False, action = 'store_true',
                     help = 'Take screenshots of IGV for each plasmid')
    igv.add_argument('--igv', required = False, default = 'igv_hidpi.sh', help = 'Path to igv.sh or igv_hidpi.sh')
    

    args = ap.parse_args()
        
    if args.medaka == False and args.porechop == False and args.nanofilt == False:
        
        sys.exit('Need to specify -M to generate consensus without binning (Medaka),'
        ' -P to bin reads on barcodes (Porechop), or -F for filtering (Nanofilt)')
        
    return args

def loadReads(inputFiles, referenceFiles, outputDir):
    '''
    This function checks to make the input, reference, and outputDir. If input or reference arguments are
    directories, this function concatenates the .fasta or fastq files to make the input reads easier to work with
    and the plasmids into a 'Plasmid Genome.' Checks that outputDir is not a file.
    '''
    if os.path.isdir(inputFiles):
        
        print('Directory was given for input reads, concatenating these to output_reads.fastq')
        
        subprocess.run(['cat %s/*.fastq > input_reads.fastq' % inputFiles], shell = True)
        reads = 'input_reads.fastq'
                                          
    if os.path.isfile(inputFiles):
        #Checks if fastq format
        if str(inputFiles).lower().endswith('.fastq'):
            reads = inputFiles
        else:
            sys.exit('Error: Reads should be in .fastq format ')
    
    if referenceFiles is not None:
        if os.path.isdir(referenceFiles):

            print('Directory was given for plasmid reference, concatenating these to output_ref.fasta')

            subprocess.run(['cat %s/*.fa* > plasmid_genome_ref.fasta' % referenceFiles], shell = True)
            reference = 'plasmid_genome_ref.fasta'


        if os.path.isfile(referenceFiles):
            #Checks if reference is already assembled in to plasmid genome
            if str(referenceFiles).lower().endswith('.fa') or str(referenceFiles).lower().endswith('.fasta'):
                reference = referenceFiles
            else:
                sys.exit('Error: Reference sequence should be in .fasta or .fa format')
                
    if os.path.isfile(outputDir):
        
        sys.exit('It looks like your output directory exists exists as a file. Please remove this file or change where' 
                 'you want the output to be written to')
                
    else:
        reference = referenceFiles
            
    return reads, reference, outputDir

def runMedaka(reads, reference, outputDir, threads, porechop, screenshot):
    '''
    Takes the input reads, reference, and output directory (after being processed through loadReads) and runs Medaka.
    Please see Dependencies.txt or https://github.com/nanoporetech/medaka for more information about installing medaka.
    Only medaka_consensus used in this pipeline.
    '''
    
    print('----------------------------------\n')
    print('Checking Medaka Args')
    print('----------------------------------\n')
    print('Input files: %s \nReference Files: %s \nOutput directory: %s \nThreads: %s \nPorechop: %s'
         % (reads, reference, outputDir, str(threads), str(porechop)))
    
    subprocess.run(['medaka_consensus -i %s -d %s -o %s -t %s -m r941_min_fast_g344'
                    % (reads, reference, outputDir, str(threads))], shell = True)
    
    processMedakaOutput(outputDir, porechop)
    
    return 


def filterReads(reads, maxLength, minLength, quality, porechop, medaka):
    '''
    Filters reads using NanoFilt. More options are available if you use NanoFilt itself please see
    https://github.com/wdecoster/nanofilt for more information. This program only allows for read length and q score 
    filtering. If filtering first, please remember to change your input path for further analysis. 
    '''
    
    subprocess.run(['NanoFilt -q %s -l %s --maxlength %s %s> filtered_output.fastq' %
                    (quality, minLength, maxLength, reads)], shell = True)
    
    new_reads = 'filtered_output.fastq'
    
    if porechop == False and medaka == False:
        
        sys.exit('Done: Filtered Reads written to "filtered_output.fastq" no Porechop or Medaka called')
    
    return new_reads

    
def processMedakaOutput(outputDir, porechop, screenshot, igv = None):
    '''
    Splits the consensus.fasta file that was generated by Medaka into individual plasmid seqs. 
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
        
    elif screenshot and igv:
        takeScreenshots(reference, outputDir, igv)
        
        
def runPorechop(reads, outputDir, barcodes, barcodeThreshold, threads, endSize, iterations, screenshot, igv):
    '''
    As of 20200206: Relies on CM build of rrwick's Porechop.
    This function sets up target directories for running Porechop with input reads, the sequence that will be used for
    custom barcoding as well as other Porechop specific parameters. Outputs that were written to various directories are
    processed with processPorechopOutput.
    '''
    
    if not os.path.isdir(outputDir):
            os.makedirs(outputDir)
    detailsOut = open(os.path.join(outputDir, 'PorechopRunDetails.txt'), 'wt')
    
    detailsOut.write('Porechop run details \n'
                     '_____________________________________\n'
                     'Input files: %s \nBarcodes: %s \nOutput Directory:'
                     '%s \nBarcode Threshold: %s \nThreads: %s \nEnd Size: %s \nInterations %s\n'
                     '_____________________________________\n'
         % (reads, barcodes, outputDir, barcodeThreshold, str(threads), str(endSize), str(iterations)))
    
    detailsOut.close()
    
    for x in range(int(iterations)):
        
        baseRun = 'porechop -i %s -b %s -BC %s -t %s --barcode_threshold %s \
        --no_split --untrimmed -v 0 --end_size %s --discard_unassigned'
        
#         if x == 0:
#             #First run is 
#             subprocess.run([baseRun % (reads, outputDir + '_0', barcodes, str(threads),
#                                        str(barcodeThreshold), str(endSize))], shell = True)
            
#         else:
        subprocess.run([baseRun % (reads, outputDir+'_'+str(x),
                                       barcodes, str(threads), str(barcodeThreshold), str(endSize))], shell = True)
            
    processPorechopOutput(outputDir, iterations, barcodes, screenshot, igv)
    
    return

def processPorechopOutput(outputDir, iterations, barcodes, screenshot, igv):
    
    #This causes problems if you name multiple output directories similar names. I recommend using the date
    #All the run info will be included in an output file, so a very descriptive name is not necessary
    paths = [path for path in os.listdir('.') if path.startswith(outputDir)]
        
    finalOutput = outputDir + '_Final'
        
    os.mkdir(finalOutput)
                
    runSummary = open(finalOutput+'/runSummary.txt', 'wt')
        
        
    for directory in range(len(paths)):
            
        subprocess.run(['cat %s/PorechopRunDetails.txt >> %s/runSummary.txt' % 
                            (paths[directory], finalOutput)], shell = True)
            
        if directory == 0:
            continue
                
        else:
            subprocess.run(['cat %s/*.fastq >> %s/allBinned.fastq' % (paths[directory], finalOutput)], shell = True)
            
    runSummary.close()
       
    runMinimap2(outputDir, barcodes, screenshot, igv)
        
    return 

        
def runMinimap2(outputDir, reference, screenshot, igv):
    
    finalReads = outputDir + '_Final/allBinned.fastq'
    
    minimap2Run = 'minimap2 -ax map-ont -L %s %s | samtools view -bq 1 | samtools sort -o %s'
    
    subprocess.run([minimap2Run % (reference, finalReads, outputDir+'_Final/filtered_sorted_reads.bam')], shell = True)
    
    subprocess.run(['samtools markdup -r %s/filtered_sorted_reads.bam %s/final_processed.bam' % 
                    (outputDir+'_Final', outputDir+'_Final')], shell = True)
    
    subprocess.run(['samtools index %s/final_processed.bam' % (outputDir + '_Final')], shell = True)
    
    #subprocess.run([samtools view -q 1 ])
    
    if screenshot:
        
        takeScreenshots(reference, outputDir, screenshot, igv)
        
    return 


def takeScreenshots(reference, outputDir, screenshot, igv):
    
    subprocess.run(['mkdir %s_Final/screenshots'% outputDir], shell = True)
    plasmidList = []
    with open(reference, 'rt') as f:#Determine what screenshots to take based on reference
        for line in f:
            if line.startswith('>'):
                plasmidList.append(line.strip('> \n'))
                
    stringForIGV = ''
    for plasmid in plasmidList:
        #####################Testing thing here
        stringForIGV += 'goto %s\nsnapshot\nsetSleepInterval 10\n' % (plasmid)
        #stringForIGV += 'goto %s\nsnapshot %s.png\n' % (plasmid, plasmid)
    
    igvBatch = open('%s_Final/igv_screenshot.bat'%outputDir, 'wt')
    
    #IGV is very VERY picky about formatting here
    template = 'new\nsnapshotDirectory %s\ngenome %s\nload %s\nmaxPanelHeight 400' \
    '\n%sexit'% (outputDir + '_Final/screenshots', reference, outputDir + '_Final/final_processed.bam', stringForIGV)
        
    igvBatch.write(template)
        
    igvBatch.close()
    
    subprocess.run(['%s -b %s' % (igv, outputDir + '_Final/igv_screenshot.bat')], shell = True)
    
       
if __name__ == "__main__":
    main()