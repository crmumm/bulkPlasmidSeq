import os
import sys
import subprocess
import bulkPlasmidSeq

from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline

def run(reads, reference, outputDir, args):
    '''
    Takes as input fastq reads or directory to reads,
    a fasta file with all the plasmid or directory pointing to the plasmid files. 
    '''
        
    if None not in (reads, reference, outputDir):
        
        runMedaka(reads, reference, outputDir, args.threads, args.model, args.igv)
            
        
    else:
        sys.exit('medaka needs input reads (-i), output directory (-o), and reference sequences (-r)')



def runMedaka(reads, reference, outputDir, threads, model, igv):
    '''
    Takes the input reads, reference, and output directory (after being processed through loadReads) and runs Medaka.
    Please see Dependencies.txt or https://github.com/nanoporetech/medaka for more information about installing medaka.
    Only medaka_consensus used in this pipeline.
    '''
    print('----------------------------------\n')
    print('Running medaka with the following arguments, more info in medaka_log.txt')
    print('----------------------------------\n')
    print('Input files: %s \nReference Files: %s \nOutput directory: %s \nThreads: %s \nModel: %s\n'
         % (reads, reference, outputDir, str(threads), model))
    
    if os.path.isfile(outputDir):
        sys.exit('Ouput directory exists')
    
    elif not os.path.exists(outputDir):
        subprocess.run(['mkdir %s' % outputDir], shell = True)
    
    with open('%s/medaka_log.txt' % outputDir, 'wt') as log:
        # Quietly run medaka, make a medaka_log
        #20210512 Added back -g (no fill option, want to see least helped results)
        subprocess.run(['medaka_consensus',
                        '-g',
                        '-i', str(reads),
                        '-d', str(reference),
                        '-o', str(outputDir),
                        '-t', str(threads),
                        '-m', str(model)],
                       stdout = log, stderr = log, check = True)
    
    processMedakaOutput(outputDir, reference, igv)
    
    return 

def filterMedakaBam(outputDir):
    '''
    Filter out secondary and low quality mapping from medaka's 'calls_to_draft.bam'
    '''
    medaka_alignment = outputDir+'/calls_to_draft.bam'
    filtered_bam = outputDir +'/filtered_alignment.bam'
    
    subprocess.run('samtools view -bq 10 -F 2048 %s | samtools sort -o %s' %
                   (medaka_alignment, filtered_bam), 
                   shell = True)
    
    subprocess.run('samtools index %s' % filtered_bam, shell = True)
    
    return filtered_bam

def processMedakaOutput(outputDir, reference, igv = None):
    '''
    Splits the consensus.fasta file that was generated by Medaka into individual plasmid seqs.
    
    Generates Emboss Needle global alignments
    '''
    #Filter 'calls_to_draft.bam file'
    output_alignmet = filterMedakaBam(outputDir)
    
    if os.path.isdir(outputDir):
        consensusFile = outputDir+'/consensus.fasta'
        
        needle_commands = []
        subprocess.run(['mkdir', '%s/consensus_sequences' % outputDir])
        destination = outputDir+'/consensus_sequences/'
        
        for record_consensus in SeqIO.parse(consensusFile, 'fasta'):
            for record_reference in SeqIO.parse(reference, 'fasta'):
                if record_reference.name in record_consensus.name: #Match the consensus to the reference
                    print('Consensus source  = ' + str(consensusFile))
                    print('A seq, reference = ' + record_reference.name)
                    print('B seq, consensus = ' + record_consensus)
                    SeqIO.write(record_consensus, destination + record_consensus.name+'_consensus.fasta', 'fasta')
                    SeqIO.write(record_reference, outputDir + '/' +record_reference.name+'_reference.fasta', 'fasta')
                    needle_commands.append(NeedleCommandline(asequence = outputDir 
                                                             + '/'+ record_reference.name+'_reference.fasta',
                                                         bsequence = destination + record_consensus.name+'_consensus.fasta',
                                                         gapopen = 10,
                                                         gapextend = 0.5,
                                                         outfile = destination + record_consensus.name + '.txt'))
    else:
        sys.exit('Medaka output directory not found')
            
    for command in needle_commands:
        subprocess.run(str(command) + ' -scircular1=True', shell = True)
                
    #Run take screenshots based on final_processed.bam    
    if igv != None:
        bulkPlasmidSeq.takeScreenshots(reference, outputDir, igv)