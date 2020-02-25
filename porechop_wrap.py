import os
import sys
#import shutil
import subprocess

import bulkPlasmidSeq

def run(reads, reference, outputDir, args):
    '''
    Checks that all the necessary inputs are provided. Note that -BC is used for building custom barcodes
    and also later as a reference for minimap. 
    '''
    if None not in (reference, outputDir):
        try:
            runPorechop(reads, outputDir, reference, args.barcode_threshold,
                    args.threads, args.end_size, args.rounds, args.screenshot, args.igv)
    
        except IOError:
            sys.exit('IO error: check that you have activated medaka conda environment/ medaka_consensus available')
            
    else:
        sys.exit('Porechop needs input reads (-i), output directory (-o), and barcodes (-BC)')


def runPorechop(reads, outputDir, barcodes, barcodeThreshold, threads, endSize, iterations, screenshot, igv):
    '''
    Relies on crmumm fork of rrwick's Porechop.
    This function sets up target directories for running Porechop with input reads, the sequence that will be used for
    custom barcoding as well as other Porechop specific parameters. Outputs that were written to various directories are
    processed with processPorechopOutput.
    '''
    if not os.path.isdir(outputDir):
            os.makedirs(outputDir)
            
    detailsOut = open(os.path.join(outputDir, 'run_details.txt'), 'wt')
    
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

        subprocess.run([baseRun % (reads, outputDir + '/' + 'porechop_'+str(x),
                                       barcodes, str(threads), str(barcodeThreshold), str(endSize))], shell = True)
            
    processPorechopOutput(outputDir, iterations, barcodes, screenshot, igv)
    
    return


def processPorechopOutput(outputDir, iterations, barcodes, screenshot, igv):
    '''
    Looks in all the directories make by Porechop and concatinates the run information into runSummary
    and cats the runs together for minimap2. Zips the fastq files that were assigned to each bin.
    '''
    
    #This causes problems if you name multiple output directories similar names. I recommend using the date
    #All the run info will be included in an output file, so a very descriptive name is not necessary
    paths = [path for path in os.listdir(outputDir) if path.startswith('porechop_')]
                
    runSummary = open(outputDir+'/run_details.txt', 'a+')
        
    for directory in range(len(paths)):
        
        subprocess.run(['cat %s/run_details.txt >> %s/run_details.txt' % 
                            (outputDir + '/' + paths[directory], outputDir)], shell = True)

        subprocess.run(['cat %s/*.fastq >> %s/all_binned.fastq' %
                            (outputDir + '/' + paths[directory], outputDir)], shell = True)

        subprocess.run(['gzip %s/*.fastq' % (outputDir + '/' + paths[directory])], shell = True)
            
    runSummary.close()
       
    runMinimap2(outputDir, barcodes, screenshot, igv)
        
    return


def runMinimap2(outputDir, reference, screenshot, igv):
    '''
    Maps the filtered reads from Porechop processing. Then filters and sorts the alignments. 
    Calls takeScreenshots if args.screenshot == True.
    
    Args:
    reference - plasmid fasta
    outoutDir - output directory for reading and writing
    igv, screenshot - Path to igv.sh, screenshot must be true
    
    '''
    
    finalReads = outputDir + '/allBinned.fastq'
    
    minimap2Run = 'minimap2 -ax map-ont -L %s %s | samtools view -bq 1 | samtools sort -o %s'
    
    subprocess.run([minimap2Run %
                    (reference, outputDir+'/all_binned.fastq', outputDir+'/filtered_sorted_reads.bam')], shell = True)
    
    subprocess.run(['samtools markdup -r %s/filtered_sorted_reads.bam %s/final_processed.bam' % 
                    (outputDir, outputDir)], shell = True)
    
    subprocess.run(['samtools index %s/final_processed.bam' % (outputDir)], shell = True)
    
    if screenshot:
        
        bulkPlasmidSeq.takeScreenshots(reference, outputDir, igv)
        
    return 
