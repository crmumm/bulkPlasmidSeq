'''
Camille Mumm - Boyle Lab 2021

Uses Medaka (https://github.com/nanoporetech/medaka) and 
NanoFilt (https://github.com/wdecoster/nanofilt) to process reads from bulk plasmid sequencing. 
'''
import argparse
import os
import sys
import subprocess
from Bio import SeqIO


def pick_submodule(args):
    
    if None not in (args.input_reads, args.reference, args.output_dir):
    
        reads, reference, outputDir = loadReads(args.input_reads, args.reference, args.output_dir, args.double)
        
    else:
        sys.exit('Missing one of the following arguments, --input(-i), --reference(-r), --output_dir(-o)')
        
    if args.filter == True:
        
        reads = filterReads(reads, outputDir, args.max_length, args.min_length, args.min_quality)
    
    if args.submod == 'biobin':
        import biobinning as submod
            
    elif args.submod == 'medaka':
        import medaka_wrap as submod
    
    submod.run(reads, reference, outputDir, args)

def main():
    '''
    See dependencies in accompanying text file.
    
    BulkPlasmidSeq usage examples:
    
    For alignment, polishing, and creating consensus without binning (Medaka):
    
        python bulkPlasmidSeq.py Medaka -i 'path/to/reads' -r 'path/to/plasmids' -o outputDirectory - t 2
    
    For binning and alignment:
    
        python bulkPlasmidSeq.py biobin -i 'path/to/reads' -o outputDirectory -r 'path/to/plasmids.fasta'
             
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
    
    #biobin Args
    biobinArgs = subparsers.add_parser('biobin')
    
    generalArgsbiobin = biobinArgs.add_argument_group('General Arguments')
    generalArgsbiobin.add_argument('-i', '--input_reads', required = True,
                          help= 'input reads (directory to reads or .fastq)')
    
    generalArgsbiobin.add_argument('-r', '--reference', required = True,
                          help = 'plasmid sequences (directory or .fa)')
    
    generalArgsbiobin.add_argument('-o', '--output_dir', required = True,
                          help = 'output directory')
    
    generalArgsbiobin.add_argument('-t', '--threads', required = False, default = 2,
                          help = 'number of threads, default 2')
    
    generalArgsbiobin.add_argument('--filter', required = False, action='store_true',
                       default = False, help = 'Filter reads before analysis')
    
    generalArgsbiobin.add_argument('-m', '--model', required = False, default = 'r941_min_high_g360',
                          help = 'Medaka consensus model, Pore/Guppy version, use medaka tools list_models for list')
    
    generalArgsbiobin.add_argument('--double', required = False, action = 'store_true', default = False,
                          help = 'Double the reference genome, great for visualization, less for consensus generation')
    generalArgsbiobin.add_argument('--trim', required = False, action = 'store_true', default = False,
                          help = 'Trim adapters from reads with Porechop')
    
    biobin = biobinArgs.add_argument_group('biobin Specific arguments')
    
    biobin.add_argument('--marker_score', required = False, default = 95, help = 'Percent score for longest unique region')
    
    biobin.add_argument('-k', '--kmer_length', required = False, default = 12)
    
    biobin.add_argument('--match', required = False, default = 3)
    biobin.add_argument('--mismatch', required = False, default = -6)
    biobin.add_argument('--gap_open', required = False, default = -10)
    biobin.add_argument('--gap_extend', required = False, default = -5)
    biobin.add_argument('--context_map', required = False, default = 0.80)
    biobin.add_argument('--fine_map', required = False, default = 0.95)
    biobin.add_argument('--max_regions', required = False, default = 3, help = 'Maximum number of regions to align/score')
    
    #IGV screenshot args Porchop
    igvbiobin = biobinArgs.add_argument_group('IGV arguments for screenshotting')
    
    igvbiobin.add_argument('--igv', required = False, default = None, help = 'Path to igv.sh')
    
    #NanoFilt args biobin
    nanofiltbiobin = biobinArgs.add_argument_group('Nanofilt arguments')
    nanofiltbiobin.add_argument('--max_length', required = False, default = 1000000,
                                  help = 'Filtering reads by maximum length')
    nanofiltbiobin.add_argument('--min_length', required = False, default = 0,
                                  help = 'Filtering reads by minimum length')
    nanofiltbiobin.add_argument('-q', '--min_quality', required = False, default = 7,
                                  help = 'Filter reads by quality score > N')
    
    biobinArgs.set_defaults(func=pick_submodule)
            
    ###***###***###***###***###***###***###***###***###***###***###***###***###***###***###***
    
    #Medaka Args, the only Medaka specific arg is -m --model, for Pore and Guppy version
    #Didn't make an argument group for that alone, add to general
    
    medakaArgs = subparsers.add_parser('medaka')
    
    generalArgsmedaka = medakaArgs.add_argument_group('General Arguments')
    generalArgsmedaka.add_argument('-i', '--input_reads', required=False,
                          help= 'input reads (directory to reads or .fastq)')
    generalArgsmedaka.add_argument('-r', '--reference', required=False,
                          help = 'plasmid sequences (directory or .fa)')
    generalArgsmedaka.add_argument('-o', '--output_dir', required =False,
                          help = 'output directory')
    generalArgsmedaka.add_argument('-t', '--threads', required = False, default = 2,
                          help = 'number of threads, default 2')
    generalArgsmedaka.add_argument('--double', required = False, action = 'store_true', default = False,
                          help = 'Double the reference genome, great for visualization, less for consensus generation')
    generalArgsmedaka.add_argument('--trim', required = False, action = 'store_true', default = False,
                          help = 'Trim adapters from reads with Porechop')
    generalArgsmedaka.add_argument('-m', '--model', required = False, default = 'r941_min_high_g360',
                          help = 'Medaka consensus model, Pore/Guppy version, use medaka tools list_models for list')
    
    generalArgsmedaka.add_argument('--filter', required = False, action='store_true',
                       default = False, help = 'Filter reads before analysis')
    
    #IGV screenshot args
    igvMedaka = medakaArgs.add_argument_group('IGV arguments for screenshotting')
    
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
    try:
        func = args.func
    except AttributeError:
        ap.error("No arguments provided. Add -h for help")
    func(args)
        
    return args


def loadReads(inputFiles, referenceFiles, outputDir, double = False):
    '''
    This function checks to make the input, reference, and outputDir. If input or reference arguments are
    directories, this function concatenates the .fasta or fastq files to make the input reads easier to work with
    and the plasmids into a 'Plasmid Genome.' Checks that outputDir is not a file.
    '''
    
    if os.path.isfile(outputDir):
        sys.exit('Output directory exists')
        
    elif not os.path.exists(outputDir):
        subprocess.run(['mkdir %s' % outputDir], shell = True)
        
    if not os.path.exists(inputFiles):
        sys.exit('Input reads not found')
        
    elif os.path.isdir(inputFiles):
        try:
            subprocess.run(['cat %s/*.fastq > %s/input_reads.fastq' % (inputFiles, outputDir)],
                       check = True, shell = True)
            
        except subprocess.CalledProcessError:
            sys.exit('Could find reads in directory, check for fastq files in input')
        
        reads = '%s/input_reads.fastq' % outputDir
        
    if trim:
        reads = trim_reads(reads)
                                          
    elif os.path.isfile(inputFiles):
        #Checks if fastq format
        if str(inputFiles).lower().endswith('.fastq') or str(inputFiles).lower().endswith('.fq'):
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
        
        with open('%s/plasmid_genome_ref.fasta' % outputDir, 'wb') as outfile:
            for f in keep:
                with open(f, "rb") as infile:
                    outfile.write(infile.read())
            
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
            double.write(str(seq[1]))
            double.write(str(seq[1]) + '\n')

        double.close()
        reference = '%s/double_reference_genome.fasta' % outputDir
         
    return reads, reference, outputDir

def trim_reads(reads, outputDir):
    '''
    Trim reads with Porechop
    '''
    
    subprocess.run(['porechop -i %s -o %s/trimmed_reads.fastq' % (reads, outputDir)])
    
    trimmed_reads = reads = '%s/trimmed_reads.fastq' % outputDir
    
    return trimmed_reads


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
        
    head, tail = os.path.split(outputDir)
    output_abs = os.path.dirname(os.path.abspath(outputDir))
    output_absolute_dir = os.path.join(output_abs, tail)
    
    igvBatch = open('%s/igv_screenshot.bat'% output_absolute_dir, 'wt')
    print('************************************************')
    head, tail = os.path.split(reference)
    reference_abs = os.path.dirname(os.path.abspath(reference))
    reference_dir = os.path.join(reference_abs, tail)
    
    #IGV is picky about formatting, watch spaces
    template = 'new\nsnapshotDirectory %s\ngenome %s\nload %s\nmaxPanelHeight 400' \
    '\n%sexit'% (output_absolute_dir + '/screenshots', reference_dir,
                 output_absolute_dir + '/filtered_alignment.bam', stringForIGV)
    
    igvBatch.write(template)
    igvBatch.close()
    
    if not igv.lower().endswith('.sh'):
        sys.exit('\n Taking screenshots needs path to igv.sh \n')
        
    #Run igv.sh -b batchfile.bat
    #Quiet the java stdout
    subprocess.run(['%s -b %s' % (igv, output_absolute_dir+ '/igv_screenshot.bat')], shell = True)
                            
    
def getFasta(file):
    '''
    Returns list of tuples (name, seq) for each plasmid in .fasta 
    '''
    
    fastaList = []
    
    for plasmid_ref in SeqIO.parse(file, "fasta"):
        fastaList.append((plasmid_ref.id, plasmid_ref.seq))
                         
    return fastaList
       
if __name__ == "__main__":
    main()