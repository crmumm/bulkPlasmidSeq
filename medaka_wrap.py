import os
import sys
import subprocess

import bulkPlasmidSeq 

def run(reads, reference, outputDir, args):
    '''
    Running Medaka, check that Porechop is false. Takes as input fastq reads or directory to reads,
    a fasta file with all the plasmid or directory pointing to the plasmid files. 
    '''
        
    if None not in (reads, reference, outputDir):
        
        runMedaka(reads, reference, outputDir, args.threads, args.model, args.screenshot, args.igv)
            
        
    else:
        sys.exit('Medaka needs input reads (-i), output directory (-o), and reference sequences (-r)')



def runMedaka(reads, reference, outputDir, threads, model, screenshot, igv):
    '''
    Takes the input reads, reference, and output directory (after being processed through loadReads) and runs Medaka.
    Please see Dependencies.txt or https://github.com/nanoporetech/medaka for more information about installing medaka.
    Only medaka_consensus used in this pipeline.
    '''
    
    print('----------------------------------\n')
    print('Running Medaka with the following arguments, more info in medaka_log.txt')
    print('----------------------------------\n')
    print('Input files: %s \nReference Files: %s \nOutput directory: %s \nThreads: %s \nModel: %s\n'
         % (reads, reference, outputDir, str(threads), model))
    
    with open('%s/medaka_log.txt' % outputDir, 'wt') as log:
        # Quietly run medaka, make a medaka_log
        subprocess.run(['medaka_consensus', '-i', str(reads), '-d', str(reference),
                        '-o', str(outputDir), '-t', str(threads), '-m', str(model)],
                       stdout = log, stderr = log, check = True)
    
    processMedakaOutput(outputDir, reference, screenshot, igv)
    
    return 


def processMedakaOutput(outputDir, reference, screenshot, igv = None):
    '''
    Splits the consensus.fasta file that was generated by Medaka into individual plasmid seqs.
    Processes calls_to_draft.bam to filter out supplemental alignments and make index of that. 
    '''
    if os.path.isdir(outputDir):
        filePath = outputDir+'/consensus.fasta'
        subprocess.run(['mkdir', '%s/consensus_sequences' % outputDir])
        consensusFasta = open(filePath, 'rt')
        for line in consensusFasta:
            if line.startswith('>'):
                header = line.split(' ')[0]  
            else:
                seq = line
                out = open('%s/consensus_sequences/%s.fasta' % (outputDir, header[1:]), 'wt')
                out.write('%s \n%s' % (header, seq))
    
    subprocess.run(['samtools view -bq 1 %s/calls_to_draft.bam > %s/final_processed.bam' %
                    (outputDir, outputDir)], shell = True)
    subprocess.run(['samtools index %s/final_processed.bam' % outputDir], shell = True)
    
                
    #Run take screenshots based on final_processed.bam    
    if screenshot:
        bulkPlasmidSeq.takeScreenshots(reference, outputDir, igv)     