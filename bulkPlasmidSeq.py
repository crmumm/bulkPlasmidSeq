'''
Camille Mumm - Boyle Lab 2021

Uses Medaka (https://github.com/nanoporetech/medaka) and 
NanoFilt (https://github.com/wdecoster/nanofilt) to process reads from bulk plasmid sequencing. 
'''
import argparse
import os
import sys
import subprocess
import yaml
from Bio import SeqIO,  Restriction
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def pick_submodule(args):
    
    if None not in (args.input_reads, args.reference, args.output_dir):
        #removed args.double - 20210728
        reads, reference, outputDir = loadReads(args.input_reads, args.reference, args.output_dir,
                                                args.restriction_enzyme, args.restriction_enzyme_table,
                                                args.trim)
        
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
    
    biobin.add_argument('-re', '--restriction_enzyme', required = False, action = 'store_true', default = False,
                                   help = 'Prompts user to input restriction enzyme cut sites to rotate reference')
    biobin.add_argument('--restriction_enzyme_table', required = False,
                        help = 'Provide a yaml with cut sites for all plasmids')
    
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
    generalArgsmedaka.add_argument('--trim', required = False, action = 'store_true', default = False,
                          help = 'Trim adapters from reads with Porechop')
    generalArgsmedaka.add_argument('-m', '--model', required = False, default = 'r941_min_high_g360',
                          help = 'Medaka consensus model, Pore/Guppy version, use medaka tools list_models for list')
    
    generalArgsmedaka.add_argument('--filter', required = False, action='store_true',
                       default = False, help = 'Filter reads before analysis')
    
    generalArgsmedaka.add_argument('-re', '--restriction_enzyme', required = False, action = 'store_true',
                                   help = 'Prompts user to input restriction enzyme cut sites to rotate reference')
    generalArgsmedaka.add_argument('--restriction_enzyme_table', required = False,
                        help = 'Provide a yaml with cut sites for all plasmids')
    
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


def loadReads(inputFiles, referenceFiles, outputDir, restriction_enzyme, restriction_enzyme_table, trim):
    '''
    This function checks to make the input, reference, and outputDir. If input or reference arguments are
    directories, this function concatenates the .fasta or fastq files to make the input reads easier to work with
    and the plasmids into a 'Plasmid Genome.' Checks that outputDir is not a file.
    
    #Removed double - 20210728
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
                                          
    elif os.path.isfile(inputFiles):
        #Checks if fastq (or fq) format 
        if str(inputFiles).lower().endswith('.fastq') or str(inputFiles).lower.endswith('.fq'):
            reads = inputFiles
        else:
            sys.exit('Error: Reads should be in .fastq (or .fq) format ')
    
    if trim:
        reads = trim_reads(reads, outputDir)
    
    if not os.path.exists(referenceFiles):
        sys.exit('Reference not found')
    

    elif os.path.isdir(referenceFiles):
        #Concatenates all files ending with .fa or .fasta into one plasmid genome
        inRefFiles = os.listdir(referenceFiles)
        with open('%s/plasmid_genome_ref.fasta' % outputDir, 'w') as outfile:
            for keepFile in inRefFiles:
                full_path_ref = os.path.join(referenceFiles, keepFile)
                if os.path.isfile(full_path_ref):
                    if full_path_ref.endswith('.fa') or full_path_ref.endswith('.fasta'):
                        with open(full_path_ref, "r") as infile:
                            for seq_record in SeqIO.parse(infile, "fasta"):
                                SeqIO.write(seq_record, outfile, "fasta")
                                
            
        outfile.close()
        infile.close()
            
        reference = '%s/plasmid_genome_ref.fasta' % outputDir


    elif os.path.isfile(referenceFiles):
        #Checks if reference is already assembled in to plasmid genome
        if str(referenceFiles).lower().endswith('.fa') or str(referenceFiles).lower().endswith('.fasta'):
            with open('%s/plasmid_genome_ref.fasta' % outputDir, 'w') as outfile:
                with open(referenceFiles, "r") as infile:
                    for seq_record in SeqIO.parse(infile, "fasta"):
                        SeqIO.write(seq_record, outfile, "fasta")
                        
            outfile.close()
            infile.close()
            
            reference = '%s/plasmid_genome_ref.fasta' % outputDir

        else:
            sys.exit('Error: Reference sequence should be in .fasta or .fa format')
            
    if restriction_enzyme:
        reference = rotate_refs(reference, outputDir, None)
        
    if restriction_enzyme_table:
        reference = rotate_refs(reference, outputDir, restriction_enzyme_table)
         
    return reads, reference, outputDir

def rotate_refs(reference, outputDir, restriction_enzyme_table):
    """
    Prompts users for restriciton enzyme cut sites and writes a new reference file
    """
    rotated_plasmids = []
    
    if restriction_enzyme_table != None:
        
        re_table = open(restriction_enzyme_table, 'r')
        plasmid_cut_site = yaml.safe_load(re_table)
        for plasmid_ref in SeqIO.parse(reference, "fasta"):
            
            try:
                enzyme = plasmid_cut_site[plasmid_ref.name]['enzyme']
                cut_site = plasmid_cut_site[plasmid_ref.name]['cut-site']
                
            except KeyError:
                #Catches if the plasmid was not provided. KeyError for enzyme
                cut_site = 0
                enzyme = None
             
            if enzyme == None and cut_site == None:
                cut_site = 0
                
            if enzyme != None:
                batch = Restriction.RestrictionBatch([enzyme])
                for provided_enzyme in batch:
                    search_site = batch.search(plasmid_ref.seq, linear = False)
                    if len(search_site[provided_enzyme]) == 0:
                        sys.exit('No cut sites found in: ' + str(plasmid_ref.name) + '. For provided enzyme: ' + str(provided_enzyme))
                    if len(search_site[provided_enzyme]) > 1:
                        sys.exit('Provided enzyme has more than 1 cut site')
                    else:
                        cut_site = int(search_site[provided_enzyme][0])
                        
            rotated_seq = plasmid_ref.seq[cut_site:] + plasmid_ref.seq[:cut_site]
            new_plasmid = SeqRecord(rotated_seq, 
                                    id = plasmid_ref.name,
                                    description = '')

            rotated_plasmids.append(new_plasmid)
            
        re_table.close()
        
    else:
        #Import all commercial enzymes
        commercial_enzymes = Restriction.CommOnly
        #for plasmid_ref in SeqIO.parse(reference, "fasta"):
        for plasmid_ref in SeqIO.parse(reference, "fasta"):
            #for single_enzyme in enzymes:
            result = commercial_enzymes.search(plasmid_ref.seq, linear = False)
            single_cutters = Restriction.RestrictionBatch()
            #Find all the single cutters
            for single_enzyme in result:
                if len(result[single_enzyme]) == 1:
                    single_cutters.add(single_enzyme)
                    print(single_enzyme)
                                
            found_enzyme = False
            attempt = 0
            cut_site = None
            while found_enzyme == False and attempt <= 2:
                selected_enzyme = input('Please enter one of the available unique cutting enzymes above for: ' \
                          + plasmid_ref.name + ': ')
                
                for available in single_cutters:
                    if selected_enzyme == str(available):
                        found_enzyme = True
                        selected_enzyme = available
                        cut_site = int(result[available][0])
                    
                attempt += 1
                
            if cut_site == None:
                sys.exit('Available enzyme not selected, exiting')
            
            print(plasmid_ref.seq)
            rotated_seq = selected_enzyme.catalyse(plasmid_ref.seq, linear = False)
            #rotated_seq.name = pl
            #rotated_seq = plasmid_ref.seq[cut_site-1:] + plasmid_ref.seq[:cut_site]
            new_plasmid = SeqRecord(rotated_seq)
            rotated_plasmids.append(new_plasmid)

    output_reference = outputDir + '/rotated_reference.fasta'
    #with open('%s/rotated_reference.fasta' % outputDir, 'wb') as rotated:
    SeqIO.write(rotated_plasmids, output_reference, 'fasta')
        
    return output_reference
        
    
def trim_reads(reads, outputDir):
    '''
    Trim reads with Porechop
    '''
    
    subprocess.run(['porechop -i %s -o %s/trimmed_reads.fastq' % (reads, outputDir)], shell = True)
    
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