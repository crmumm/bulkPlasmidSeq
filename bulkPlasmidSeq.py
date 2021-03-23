'''
Camille Mumm - Boyle Lab rotation 2020

Uses Medaka (https://github.com/nanoporetech/medaka), Porechop (https://github.com/rrwick/Porechop),
and NanoFilt (https://github.com/wdecoster/nanofilt) to process reads from bulk plasmid sequencing. 
'''
import argparse
import os
import sys
import subprocess


def pick_submodule(args):
    
    if None not in (args.input_reads, args.reference, args.output_dir):
    
        reads, reference, outputDir = loadReads(args.input_reads, args.reference, args.output_dir, args.double)
    
    else:
        sys.exit('Missing one of the following arguments, --input(-i), --reference(-r), --output_dir(-o)')
        
    if args.filter == True:
        
        reads = filterReads(reads, outputDir, args.max_length, args.min_length, args.min_quality)
    
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
    
        python bulkPlasmidSeq.py Medaka -i 'path/to/reads' -r 'path/to/plasmids' -o outputDirectory - t 2
    
    For binning and alignment (Porechop):
    
        python bulkPlasmidSeq.py Porechop -i 'path/to/reads' -o outputDirectory -r 'path/to/plasmids.fasta'
        
        python bulkPlasmidSeq.py Porechop -i 'path/to/reads' -o outputDirectory -r 'path/to/plasmids' --end_size 500 -N 8        
    '''
    #Modeling submodule picking after arq5x's poretools structure
    #Get args, picking of submodule runs in here
    args = getArgs()
            
    
def getArgs():
    '''
    Parses arguments:
        - See -h, --help
        
    Many of the arguments (ex. input, reference, output) are redundant because they go with both submodules.
    
    '''
    ap = argparse.ArgumentParser()
    subparsers = ap.add_subparsers(title='[sub-commands]', dest='submod')

    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***
    
    #Porechop Args
    porechopArgs = subparsers.add_parser('Porechop')
    
    generalArgsPorechop = porechopArgs.add_argument_group('General Arguments')
    generalArgsPorechop.add_argument('-i', '--input_reads', required = True,
                          help= 'input reads (directory to reads or .fastq)')
    generalArgsPorechop.add_argument('-r', '--reference', required = True,
                          help = 'plasmid sequences (directory or .fa)')
    generalArgsPorechop.add_argument('-o', '--output_dir', required = True,
                          help = 'output directory')
    generalArgsPorechop.add_argument('-t', '--threads', required = False, default = 2,
                          help = 'number of threads, default 2')
    generalArgsPorechop.add_argument('--filter', required = False, action='store_true',
                       default = False, help = 'Filter reads before analysis')
    generalArgsPorechop.add_argument('--double', required = False, action = 'store_true', default = False,
                          help = 'Double the reference genome, great for visualization, less for consensus generation')
    
    porechop = porechopArgs.add_argument_group('Porechop Specific arguments')
    porechop.add_argument('--barcode_threshold', required = False, default = 75)
    porechop.add_argument('--end_size', required = False, default = 250)
    porechop.add_argument('--rounds', '--porechop_iterations', required = False, default = 1,
                              help = 'Number of rounds of Porechop binning')
    
    #IGV screenshot args Porchop
    igvPorechop = porechopArgs.add_argument_group('IGV arguments for screenshotting')
    igvPorechop.add_argument('--screenshot', required = False, action = 'store_true',
                     help = 'Take screenshots of IGV for each plasmid')
    igvPorechop.add_argument('--igv', required = False, default = None, help = 'Path to igv.sh')
    
    #NanoFilt args Porechop
    nanofiltPorechop = porechopArgs.add_argument_group('Nanofilt arguments')
    nanofiltPorechop.add_argument('--max_length', required = False, default = 100000,
                                  help = 'Filtering reads by maximum length')
    nanofiltPorechop.add_argument('--min_length', required = False, default = 0,
                                  help = 'Filtering reads by minimum length')
    nanofiltPorechop.add_argument('-q', '--min_quality', required = False, default = 7,
                                  help = 'Filter reads by quality score > N')
    
    porechopArgs.set_defaults(func=pick_submodule)
            
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***
    
    #Medaka Args, the only Medaka specific arg is -m --model, for Pore and Guppy version
    #Didn't make an argument group for that alone, add to general
    
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
    generalArgsMedaka.add_argument('--double', required = False, action = 'store_true', default = False,
                          help = 'Double the reference genome, great for visualization, less for consensus generation')
    generalArgsMedaka.add_argument('-m', '--model', required = False, default = 'r941_min_high_g344',
                          help = 'Medaka consensus model, Pore/Guppy version, use medaka tools list_models for list')
    
    generalArgsMedaka.add_argument('--filter', required = False, action='store_true',
                       default = False, help = 'Filter reads before analysis')
    
    #IGV screenshot args Porchop
    igvMedaka = medakaArgs.add_argument_group('IGV arguments for screenshotting')
    igvMedaka.add_argument('--screenshot', required = False, action = 'store_true',
                     help = 'Take screenshots of IGV for each plasmid')
    igvMedaka.add_argument('--igv', required = False, default = None, help = 'Path to igv.sh')
    
    
    #NanoFilt args
    nanofiltMedaka= medakaArgs.add_argument_group('Nanofilt arguments')
    nanofiltMedaka.add_argument('--max_length', required = False, default = 100000,
                                help = 'Filtering reads by maximum length, default is 1Mb')
    nanofiltMedaka.add_argument('--min_length', required = False, default = 0,
                                help = 'Filtering reads by minimum length, default is 0')
    nanofiltMedaka.add_argument('-q', '--min_quality', required = False, default = 7,
                                help = 'Filter reads by quality score > N, default is 7')
    
    medakaArgs.set_defaults(func=pick_submodule)
    
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***

    args = ap.parse_args()
    args.func(args)

    if args.screenshot is True and args.igv is None:
        sys.exit('Please specify a path to igv.sh for taking screenshots')
        
    return args


def loadReads(inputFiles, referenceFiles, outputDir, double = False):
    '''
    This function checks to make the input, reference, and outputDir. If input or reference arguments are
    directories, this function concatenates the .fasta or fastq files to make the input reads easier to work with
    and the plasmids into a 'Plasmid Genome.' Checks that outputDir is not a file.
    '''
    
    if os.path.isfile(outputDir):
        sys.exit('Output directory exists exists as a file.')
        
    elif not os.path.exists(outputDir):
        subprocess.run(['mkdir %s' % outputDir], shell = True)
        
    if not os.path.exists(inputFiles):
        sys.exit('Input reads not found')
        
    elif os.path.isdir(inputFiles):
        #print('Directory was given for input reads, concatenating these to output_reads.fastq \n')
        
        try:
            subprocess.run(['cat %s/*.fastq > %s/input_reads.fastq' % (inputFiles, outputDir)],
                       check = True, shell = True)
            
        except subprocess.CalledProcessError:
            sys.exit('Could find reads in directory, check for fastq files in input')
        
        reads = '%s/input_reads.fastq' % outputDir
                                          
    elif os.path.isfile(inputFiles):
        #Checks if fastq format
        if str(inputFiles).lower().endswith('.fastq'):
            reads = inputFiles
        else:
            sys.exit('Error: Reads should be in .fastq format ')
    
    if not os.path.exists(referenceFiles):
        sys.exit('Reference not found')
    

    elif os.path.isdir(referenceFiles):
        #print('Directory was given for plasmid reference, concatenating these to output_ref.fasta \n')
        
        inRefFiles = os.listdir(referenceFiles)
        keep = []
        [keep.append(referenceFiles + keepFile) for keepFile in inRefFiles \
         if keepFile.endswith('.fa') or keepFile.endswith('.fasta')]
        
        print(keep)
        
        with open('%s/plasmid_genome_ref.fasta' % outputDir, 'wb') as outfile:
            for f in keep:
                with open(f, "rb") as infile:
                    outfile.write(infile.read())
        
        #try:
            #subprocess.run(['cat %s/> %s/plasmid_genome_ref.fasta' %
            #(referenceFiles, outputDir)], shell = True, check = True)
            
        #except subprocess.CalledProcessError:
            #sys.exit('Could not find reference plasmid sequences')
            
        reference = '%s/plasmid_genome_ref.fasta' % outputDir


    elif os.path.isfile(referenceFiles):
        #Checks if reference is already assembled in to plasmid genome
        if str(referenceFiles).lower().endswith('.fa') or str(referenceFiles).lower().endswith('.fasta'):
            reference = referenceFiles
        else:
            sys.exit('Error: Reference sequence should be in .fasta or .fa format')
    
    if double:
        print('\n Double is True, duplicating the reference sequence, used for visualization  of plasmid fragmentation\n')
        #Duplicates the plasmids in the reference, used as an argument with Medaka. 
        #Great for visualization, not consensus generation/polishing
        seqs = getFasta(reference)
        double = open('%s/double_reference_genome.fasta' % outputDir, 'a+')
        for seq in seqs:
            double.write('>' + seq[0] + '\n')
            double.write(seq[1])
            double.write(seq[1] + '\n')

        double.close()
        reference = '%s/double_reference_genome.fasta' % outputDir
         
    return reads, reference, outputDir


def filterReads(reads, outputDir, maxLength, minLength, quality):
    '''
    Filters reads using NanoFilt. More options are available if you use NanoFilt itself please see
    https://github.com/wdecoster/nanofilt for more information. This program only allows for read length and q score 
    filtering. If filtering first, please remember to change your input path for further analysis. 
    '''
    print('----------------------------------\n')
    print('Filtering reads with NanoFilt: \nMin Length: %s \nMax Length: %s \nQ score: %s \n'
          % (minLength, maxLength, quality))
    print('Output: filtered_reads.fastq')
    print('----------------------------------\n')
    
    subprocess.run(['NanoFilt -q %s -l %s --maxlength %s %s> %s/filtered_reads.fastq' %
                    (quality, minLength, maxLength, reads, outputDir)], shell = True)
    
    new_reads = 'filtered_reads.fastq'
    
    return new_reads


def takeScreenshots(reference, outputDir, igv):
    '''
    Uses batch screenshot functionality of IGV to take images of the alignment for each plasmid:
    Makes a .bat file in output dir.
    Args:
    reference - plasmid fasta
    outoutDir - output directory for reading and writing
    igv - Path to igv.sh
    '''
    print('----------------------------------\n')
    print('Taking screenshots in IGV - writing images to %s/screenshots' % outputDir)
    print('----------------------------------\n')
    
    subprocess.run(['mkdir', '%s/screenshots' % outputDir])
    plasmidList = []
    with open(reference, 'rt') as f:#Determine what screenshots to take based on reference
        for line in f:
            if line.startswith('>'):
                plasmidList.append(line.strip('> \n'))
                
    stringForIGV = ''
    for plasmid in plasmidList:
        #Sleeping between commands as a workaround of IGV bug
        stringForIGV += 'goto %s\nsnapshot\nsetSleepInterval 10\n' % (plasmid)
    
    igvBatch = open('%s/igv_screenshot.bat'% outputDir, 'wt')
    
    #IGV is picky about formatting, watch spaces
    template = 'new\nsnapshotDirectory %s\ngenome %s\nload %s\nmaxPanelHeight 400' \
    '\n%sexit'% (outputDir + '/screenshots', reference, outputDir + '/final_processed.bam', stringForIGV)
    
    igvBatch.write(template)
    igvBatch.close()
    
    if not igv.lower().endswith('.sh'):
        sys.exit('\n Taking screenshots needs path to igv.sh \n')
        
    #Run igv.sh -b batchfile.bat
    #Quiet the java stdout 
    
    subprocess.run(['%s' % igv, ' -b',  ' %s' % (outputDir+ '/igv_screenshot.bat')])
                            
    
def getFasta(file):
    '''
    Returns list of tuples (name, seq) for each plasmid in .fasta 
    '''
    file = open(file, 'rt')
    name=''
    seq=''
    fastaList = []
    for line in file:
        # Capture the next header, report what we have, and update
        if line.startswith('>') and seq: #not first seq
            name = name[1:] #removes the carrot
            fastaList.append((name, seq))
            name=line.strip()
            seq=''
            # Get to the first header
        elif line.startswith('>'):  #first seq
            name=line.strip()
        # Just add sequence if it is the only thing there
        else:
            seq+=line.strip()
        # At the end, return the last entries
    if name and seq: #last seq
        name = name[1:]
        fastaList.append((name, seq))                           
    return fastaList
       
if __name__ == "__main__":
    main()