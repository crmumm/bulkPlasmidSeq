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


def pick_submodule(args):
    
    if None not in (args.input_reads, args.reference, args.output_dir):
    
        reads, reference, outputDir = loadReads(args.input_reads, args.reference, args.output_dir)
    
    else:
        sys.exit('Missing one of the following arguments, --input(-i), --reference(-r), --output_dir(-o)')
    
    if args.submod == 'Porechop':
        import porechop_wrap as submod
            
    elif args.submod == 'Medaka':
            
        import medaka_wrap as submod
    
    submod.run(reads, reference, outputDir, args)

def main():
    '''
    See dependencies in accompanying text file.
    
    BulkPlasmidSeq usage examples:
    
    For alignment, polishing, and creating consensus without binning (Medaka):
    
        python bulkPlasmidSeq.py Medaka -i 'path/to/reads' -r 'path/to/plasmids' -o outputDirectory
    
    For binning and alignment (Porechop):
    
        python bulkPlasmidSeq.py Porechop -i 'path/to/reads' -o outputDirectory -BC 'path/to/plasmids'
        
        python bulkPlasmidSeq.py Porechop -i 'path/to/reads' -o outputDirectory -BC 'path/to/plasmids' --end_size 1000 -N 3        
    '''
    #20200220 Notes, I'm trying to make the Porechop and Medaka functionalities more separate 
    #Modeling this after arq5x's poretools submodule calling structure
    args = getArgs()
            
    
def getArgs():
    '''
    Parses arguments:
        - See -h, --help
    
    '''
    ap = argparse.ArgumentParser()
    subparsers = ap.add_subparsers(title='[sub-commands]', dest='submod')
    
    #NanoFilt args
    nanofilt = ap.add_argument_group('Nanofilt arguments')
    nanofilt.add_argument('--max_length', required = False, default = 100000, help = 'Filtering reads by maximum length')
    nanofilt.add_argument('--min_length', required = False, default = 0,  help = 'Filtering reads by minimum length')
    nanofilt.add_argument('-q', '--min_quality', required = False, default = 7, help = 'Filter reads by quality score > N')
    
        
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***

    
    #Porechop Args
    porechopArgs = subparsers.add_parser('Porechop')
    
    generalArgs = porechopArgs.add_argument_group('General Arguments')
    generalArgs.add_argument('-i', '--input_reads', required = True,
                          help= 'input reads (directory to reads or .fastq)')
    generalArgs.add_argument('-r', '--reference', required = True,
                          help = 'plasmid sequences (directory or .fa)')
    generalArgs.add_argument('-o', '--output_dir', required = True,
                          help = 'output directory')
    generalArgs.add_argument('-t', '--threads', required = False, default = 2,
                          help = 'number of threads, default 2')
    
    porechop = porechopArgs.add_argument_group('Porechop Specific arguments')
    porechop.add_argument('--barcode_threshold', required = False, default = 75)
    porechop.add_argument('--end_size', required = False, default = 250)
    porechop.add_argument('-rounds', '--porechop_iterations', required = False, default = 1,
                              help = 'Number of rounds of Porechop binning')
    
    #IGV screenshot args Porchop
    igvPorechop = porechopArgs.add_argument_group('IGV arguments for screenshotting')
    igvPorechop.add_argument('--screenshot', required = False, action = 'store_true',
                     help = 'Take screenshots of IGV for each plasmid')
    igvPorechop.add_argument('--igv', required = False, default = None, help = 'Path to igv.sh')
    
    porechopArgs.set_defaults(func=pick_submodule)
            
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***
    
    #Medaka takes no additional args but I'd like to separate them
    medakaArgs = subparsers.add_parser('Medaka')
    generalArgsMedaka = medakaArgs.add_argument_group('General Arguments')
    generalArgsMedaka.add_argument('-i', '--input_reads', required=False,
                          help= 'input reads (directory to reads or .fastq)')
    generalArgsMedaka.add_argument('-r', '--reference', required=False,
                          help = 'plasmid sequences (directory or .fa)')
    generalArgsMedaka.add_argument('-o', '--output_dir', required =False,
                          help = 'output directory')
    generalArgsMedaka.add_argument('-t', '--threads', required = False, default = 2,
                          help = 'number of threads, default 2')
    
    generalArgsMedaka.add_argument('--info', required = False, help = 
                           'Medaka takes no additional args, thanks for using this submodule')
    
    #IGV screenshot args Porchop
    igvMedaka = medakaArgs.add_argument_group('IGV arguments for screenshotting')
    igvMedaka.add_argument('--screenshot', required = False, action = 'store_true',
                     help = 'Take screenshots of IGV for each plasmid')
    igvMedaka.add_argument('--igv', required = False, default = None, help = 'Path to igv.sh')
    
    medakaArgs.set_defaults(func=pick_submodule)
    
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***

    args = ap.parse_args()
    args.func(args)

    if args.screenshot is True and args.igv is None:
        sys.exit('Please specify a path to igv.sh for taking screenshots')
        
    return args

def loadReads(inputFiles, referenceFiles, outputDir):
    print(inputFiles, referenceFiles, outputDir)
    '''
    This function checks to make the input, reference, and outputDir. If input or reference arguments are
    directories, this function concatenates the .fasta or fastq files to make the input reads easier to work with
    and the plasmids into a 'Plasmid Genome.' Checks that outputDir is not a file.
    '''
    if os.path.exists(inputFiles) == False:
        sys.exit('Input reads not found')
        
    elif os.path.isdir(inputFiles):
        
        print('Directory was given for input reads, concatenating these to output_reads.fastq')
        
        subprocess.run(['cat %s/*.fastq > input_reads.fastq' % inputFiles], shell = True)
        reads = 'input_reads.fastq'
                                          
    elif os.path.isfile(inputFiles):
        #Checks if fastq format
        if str(inputFiles).lower().endswith('.fastq'):
            reads = inputFiles
        else:
            sys.exit('Error: Reads should be in .fastq format ')
    
    if os.path.exists(referenceFiles) == False:
        sys.exit('Reference not found')
    

    elif os.path.isdir(referenceFiles):

        print('Directory was given for plasmid reference, concatenating these to output_ref.fasta')

        subprocess.run(['cat %s/*.fa* > plasmid_genome_ref.fasta' % referenceFiles], shell = True)
        reference = 'plasmid_genome_ref.fasta'


    elif os.path.isfile(referenceFiles):
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


def takeScreenshots(reference, outputDir, igv):
    '''
    Uses batch screenshot functionality of IGV to take images of the alignment for each plasmid:
    
    Args:
    reference - plasmid fasta
    outoutDir - output directory for reading and writing
    igv - Path to igv.sh
    '''
    
    subprocess.run(['mkdir %s/screenshots'% outputDir], shell = True)
    plasmidList = []
    with open(reference, 'rt') as f:#Determine what screenshots to take based on reference
        for line in f:
            if line.startswith('>'):
                plasmidList.append(line.strip('> \n'))
                
    stringForIGV = ''
    for plasmid in plasmidList:
        #Sleeping between commands as a workaround of Java bug
        stringForIGV += 'goto %s\nsnapshot\nsetSleepInterval 10\n' % (plasmid)
    
    igvBatch = open('%s/igv_screenshot.bat'% outputDir, 'wt')
    
    #IGV is picky about formatting, watch spaces
    template = 'new\nsnapshotDirectory %s\ngenome %s\nload %s\nmaxPanelHeight 400' \
    '\n%sexit'% (outputDir + '/screenshots', reference, outputDir + '/final_processed.bam', stringForIGV)
    
    igvBatch.write(template)
    igvBatch.close()
    #Run igv.sh -b batchfile.bat
    subprocess.run(['%s -b %s' % (igv, outputDir + '/igv_screenshot.bat')], shell = True)
    
       
if __name__ == "__main__":
    main()