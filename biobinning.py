from Bio import Align,SeqIO
from collections import defaultdict, Counter
import sys
#import matplotlib.pyplot as plt
import medaka_wrap
import bulkPlasmidSeq

def run(fastq_reads, reference, output_directory, args):
    
    k = int(args.kmer_length)
    
    run_medaka_binned(fastq_reads, reference, k, output_directory, args,
                      args.match, args.mismatch, args.gap_open, args.gap_extend, args.context_map, args.fine_map)
    
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
    #print('Found %d plasmids' % len(unique_kmers))
    for plasmid in kmers:
        save_start = kmers[plasmid][0][1]
        save_end = kmers[plasmid][-1][1]
        moving = save_start
        for i in range(len(kmers[plasmid])):
            if moving == kmers[plasmid][i][1]:
                moving +=1
                if moving == save_end:
                    intervals[plasmid].append((save_start, moving))
            else:
                intervals[plasmid].append((save_start, moving))
                moving = kmers[plasmid][i][1]

    #Get some maker definitions
    best_markers = defaultdict(list)
    for x in intervals:
        #Keep the longest unique interval 
        #Format = {'plasmid_name': [(start, end), len, unique_seq, context_seq]}
        best_markers[x]= [(0, 1), 1, 'None', 'None']
        for y in intervals[x]:
            length_marker = y[1]-y[0]
            #Calculate unique length based on k
            #Len + 1 = k + uniq - 1
            uniq_length = length_marker + 2 - k
            if uniq_length > best_markers[x][1]:
                marker_start = round(y[0]+k-2)
                marker_end = round(y[1])
                marker_seq = seqs[x][marker_start:marker_end]
                if marker_start-10 < 0 or marker_end+10 > len(seqs[x]):
                    print('End of plasmid')
                else:
                    context_seq = seqs[x][marker_start-10:marker_end+10]
                
                best_markers[x] = [(marker_start, marker_end), uniq_length, marker_seq, context_seq]
                    
    return best_markers

def align_reads(fastq_reads, fasta_ref, k,  match, mismatch, gap_open, gap_extend, context_score, fine_alignment_score):
    reads_dict = defaultdict(list)
    aligner = define_aligner(match, mismatch, gap_open, gap_extend)
    #Get the context and unique sequence from define makers
    markers = define_markers(fasta_ref, k)
    for plasmid in markers:
        context_fwd_marker = markers[plasmid][3]
        context_rev_marker = reverse_complement(context_fwd_marker)
        #Match with context
        best_score = aligner.score(context_fwd_marker, context_fwd_marker)
        
        #Iterate through each plasmid, then each read
        for record in SeqIO.parse(fastq_reads, 'fastq'):
            #Determines if forward or backward is better without doing alignment
            scores = [aligner.score(record.seq, context_fwd_marker),
                     aligner.score(record.seq, context_rev_marker)]
            top_score = (max(scores), scores.index(max(scores)))
            
            if top_score[0] >= context_score*best_score:
                #Define top possible score for unique vs unique
                best_fine_score = aligner.score(markers[plasmid][2], markers[plasmid][2])
                
                #If the forward marker had the best score
                if top_score[1] == 0:
                    top_fwd_alignment = sorted(aligner.align(record.seq, context_fwd_marker))[0]
                    #Get the sequence that the aligned to the larger region
                    subseq = record.seq[top_fwd_alignment.aligned[0][0][0]:top_fwd_alignment.aligned[0][-1][-1]]
                    #Possible add length subseq check
                    if aligner.score(subseq, markers[plasmid][2]) >= fine_alignment_score*best_fine_score:
                        reads_dict[plasmid].append(record)
                        
                #If the reverse marker has the best score
                if top_score[1] == 1:
                    #Get the sequence that the aligned to the larger region
                    top_rev_alignment = sorted(aligner.align(record.seq, context_fwd_marker))[0]
                    subseq = record.seq[top_rev_alignment.aligned[0][0][0]:top_rev_alignment.aligned[0][-1][-1]]
                    if aligner.score(subseq, reverse_complement(markers[plasmid][2])) >= fine_alignment_score*best_fine_score:
                        reads_dict[plasmid].append(record)
                
    return reads_dict

def write_bins(fastq_reads, reference, k, output_directory, match = 3,
               mismatch = -6, gap_open = -10, gap_extend = -5, context_map = 0.80, fine_map = 0.95):
    
    reads_dict = align_reads(fastq_reads, reference, k,  match, mismatch, gap_open, gap_extend, context_map, fine_map)
    output_directory = output_directory.strip('/')
    output_bins = []
    for plasmid in reads_dict:
        output_bins.append(output_directory+'/%s' % plasmid)
        SeqIO.write(reads_dict[plasmid], output_directory+'/%s.fastq' % plasmid, "fastq")
        
    return output_bins

def run_medaka_binned(fastq_reads, reference, k, output_directory, args,
                      match = 3, mismatch = -6, gap_open = -10, gap_extend = -5, context_map = 0.80, fine_map = 0.95):
    
    output_bins = write_bins(fastq_reads, reference, k, output_directory, match, mismatch, gap_open, gap_extend,
                            context_map, fine_map)
    
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
                                       