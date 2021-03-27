from Bio import Align,SeqIO
from collections import defaultdict, Counter
import sys
#import matplotlib.pyplot as plt
import medaka_wrap

def run(fastq_reads, reference, output_directory, args):
    
    k = int(args.kmer_length)
    
    run_medaka_binned(fastq_reads, reference, k, output_directory, args, args.match, args.mismatch, args.gap_open, args.gap_extend)
    
    return

def define_aligner(match, mismatch, open_gap, extend):
    #Nanopore read Scoring
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 3
    aligner.mismatch_score = -6
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -5
    
    return aligner

def getFasta(file):
    name=''
    seq=''
    plasmids = defaultdict(list)
    with open(file, 'rt') as f:    
        for line in f:
            if line.startswith('>') and seq: 
                name = name[1:]
                plasmids[name] = seq
                name=line.strip()
                seq=''
            elif line.startswith('>'):
                name=line.strip()
            else:
                seq+=line.strip()
        if name and seq: #last seq
            name = name[1:]
            plasmids[name] = seq
                
    return plasmids

def reverse_complement(seq):
    """
    Swapping adenine with thymine and guanine with cytosine.
    Reversing newly generated string
    """
    mapping = str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)[::-1]

def build_kmers(file, k):
    seqs = getFasta(file)
    plasmids = defaultdict(list)
    for x in seqs:
        plasmids[x] = [(seqs[x][i : i + k],i) for i in range(0, len(seqs[x]) - k + 1)]
    return plasmids

def unique_kmers(file, k):
    plasmids = build_kmers(file, k)
    unique_kmers = defaultdict(list)
    for x in plasmids:
        op = [bc[0] for other in plasmids for bc in plasmids[other] if other != x]
        for marker in plasmids[x]:
            if marker[0] not in op and reverse_complement(marker[0]) not in op:
                unique_kmers[x].append(marker)
    return unique_kmers

def define_markers(file, k):
    kmers = unique_kmers(file, k)
    seqs = getFasta(file)
    intervals = defaultdict(list)
    for plasmid in kmers:
        print('New plasmid')
        print(kmers[plasmid])
        save_start = kmers[plasmid][0][1]
        print('Start: ' + str(save_start))
        save_end = kmers[plasmid][-1][1]
        print('End: '+ str(save_end))
        moving = save_start
        for i in range(len(kmers[plasmid])):
            print('Len = ' + str(len(kmers[plasmid])))
            print('Kmer index ' + str(kmers[plasmid][i][1]))
            if moving == kmers[plasmid][i][1]:
                print(moving)
                moving +=1
                if moving == save_end:
                    print("At end interval")
                    intervals[plasmid].append((save_start, moving))
            else:
                print('In else')
                print(moving)
                intervals[plasmid].append((save_start, moving))
                print(intervals)
                save_start = kmers[plasmid][i+1][1]
                moving = kmers[plasmid][i+1][1]

    #Get some maker definitions
    best_markers = defaultdict(list)
    for x in intervals:
        #Keep the longest unique interval 
        #Format = {'plasmid_name': [(start, end), len, seq]}
        best_markers[x]= [(0, 1), 1]
        for y in intervals[x]:
            length_marker = y[1]-y[0]
            #Calculate unique length based on k
            #Len + 1 = k + uniq - 1
            uniq_length = length_marker + 2 - k
            #print(best_markers[x][1])
            if uniq_length > best_markers[x][1]:
                marker_start = round(y[0] + k*0.5)
                marker_end = round(y[1]+k*0.5)
                marker_seq = seqs[x][marker_start:marker_end]
                best_markers[x] = [(marker_start, marker_end), uniq_length, marker_seq]
    print(best_markers)                         
    return best_markers

def align_reads(fastq_reads, reference, k, match = 3, mismatch = -6, gap_open = -10, gap_extend = -5):
    # Reads binned by BC
    reads_dict = defaultdict(list)
    aligner = define_aligner(match, mismatch, gap_open, gap_extend)
    markers = define_markers(reference, k)
    #record_dict = SeqIO.index(read_libray, 'fastq')
    for record in SeqIO.parse(fastq_reads, 'fastq'):
        for plasmid in markers:
            fwd_marker = markers[plasmid][2]
            best_score = aligner.score(fwd_marker, fwd_marker)
            top_score = max(aligner.score(record.seq, fwd_marker), aligner.score(record.seq, reverse_complement(fwd_marker)))
            ##print(top_score)
            if top_score >= best_score*0.95:
                reads_dict[plasmid].append(record)
                          
    return reads_dict

def write_bins(fastq_reads, reference, k, output_directory, match = 3, mismatch = -6, gap_open = -10, gap_extend = -5):
    
    reads_dict = align_reads(fastq_reads, reference, k, match, mismatch, gap_open, gap_extend)
    output_directory = output_directory.strip('/')
    output_bins = []
    for plasmid in reads_dict:
        output_bins.append(output_directory+'/%s' % plasmid)
        SeqIO.write(reads_dict[plasmid], output_directory+'/%s.fastq' % plasmid, "fastq")
        
    return output_bins

def run_medaka_binned(fastq_reads, reference, k, output_directory, args,
                      match = 3, mismatch = -6, gap_open = -10, gap_extend = -5):
    
    output_bins = write_bins(fastq_reads, reference, k, output_directory,
                             match = 3, mismatch = -6, gap_open = -10, gap_extend = -5)
    
    with open(reference) as ref_list:
        #for record, plasmid_bin in zip(SeqIO.parse(ref_list, 'fasta'), output_bins):
        for record in SeqIO.parse(ref_list, 'fasta'):
            for plasmid_bin in output_bins:
                plasmid_name = plasmid_bin.split('/')[-1]
                if plasmid_name == record.name:
                    #for plasmid_bin in output_bins:
                    single_reference = SeqIO.write(record, '%s.fasta' % plasmid_bin, 'fasta')
                    medaka_wrap.runMedaka(plasmid_bin + '.fastq', '%s.fasta'% plasmid_bin,
                                      plasmid_bin+'_consensus',
                                      args.threads,
                                      args.model, args.screenshot, args.igv)
                                       