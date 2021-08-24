from Bio import Align,SeqIO
from collections import defaultdict, Counter
import sys
from operator import itemgetter
import medaka_wrap
import bulkPlasmidSeq

def run(fastq_reads, reference, output_directory, args):
    
    k = int(args.kmer_length)
    run_medaka_binned(fastq_reads, reference, k, output_directory, args,
                      args.match, args.mismatch, args.gap_open, args.gap_extend,
                      args.context_map, args.fine_map, int(args.max_regions))
    
    return

def define_aligner(match = 3, mismatch = -6, open_gap = -10, extend = -5):
    '''
    Nanopore read scoring
     - match = 3
     - mismatch = -6
     - open_gap = -10
     - extend_gap = -5
    
    Returns Biopython PairwiseAligner
     - mode = local
    '''
    
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = int(match) 
    aligner.mismatch_score = int(mismatch)
    aligner.open_gap_score = int(open_gap)
    aligner.extend_gap_score = int(extend)
    
    return aligner

def getFasta(file):
    
    plasmids = defaultdict(list)
    for plasmid_ref in SeqIO.parse(file, "fasta"):
        plasmids[plasmid_ref.id] = plasmid_ref.seq
                
    return plasmids

def reverse_complement(seq):
    """
    Swapping adenine with thymine and guanine with cytosine.
    Reversing newly generated string
    """
    mapping = str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)[::-1]

def build_kmers(file, k):
    '''
    Called by unique_kmers. 
    
    Returns dictionary of all kmers in plasmids.
    
    '''
    seqs = getFasta(file)
    plasmids = defaultdict(list)
        
    for x in seqs:
        plasmids[x] = [(seqs[x][i : i + k],i) for i in range(0, len(seqs[x]) - k + 1)]
    return plasmids

def unique_kmers(file, k):
    '''
    Called by define markers, calls build_kmers. 
    
    Returns a default dict of unique kmers
    '''
    
    plasmids = build_kmers(file, k)
    unique_kmers = defaultdict(list)
    for x in plasmids:
#         for other in plasmids:
#             for a in plasmids[other]:
#                 print(a[0])
            
        op = [bc[0] for other in plasmids for bc in plasmids[other] if other != x]
        for marker in plasmids[x]:
            if marker[0] not in op and reverse_complement(str(marker[0])) not in op:
                unique_kmers[x].append(marker)
    
    return unique_kmers

def define_markers(file, k, max_regions):
    '''
    Called by align_reads. Calls unique_kmers to return all the unique regions for plasmids.
    Processes the unique regions to get the start, end, length, sequence, and context sequence for the 
    regions to align to. These region are then sorted by length and the top [max_regions] are selected.
    
    If the regions are long enough (> ~50bp) would be better to use medaka to distinguish plasmids.
    Due to high identity scores required for binning, longer regions return fewer reads.
         - Lower score context and fine map score 
         - Use medaka
    
    Number reads * max_regions 
    
    Returns default dict (list) of best [max_regions, default 3] regions.
    Format {'plasmid_name': [(start, end), len, unique_seq, context_seq]}
    
    '''
    kmers = unique_kmers(file, k)
    seqs = getFasta(file)
    intervals = defaultdict(list)
    #print('Found %d plasmids' % len(unique_kmers))
    for plasmid in kmers:
        save_start = kmers[plasmid][0][1]
        save_end = kmers[plasmid][-1][1]
        moving = save_start
        for i in range(len(kmers[plasmid])):
            #print(save_start, moving, save_end, kmers[plasmid][i][1])
            if moving == kmers[plasmid][i][1]:
                moving +=1
                if moving == save_end:
                    intervals[plasmid].append((save_start, moving))
            else:
                intervals[plasmid].append((save_start, moving))
                moving = kmers[plasmid][i+1][1]
                save_start = moving
                            
    #Get some marker definitions
    best_markers = defaultdict(list)
    
    for x in intervals:
        #Keep the longest unique interval 
        #Format = {'plasmid_name': [(start, end), len, unique_seq, context_seq]}
        #best_markers[x]= [(0, 1), 1, 'None', 'None']
        for y in intervals[x]:
            length_marker = y[1]-y[0]
#             print(length_marker)
#             print(y)
            #print(k)
            #Calculate unique length based on k
            uniq_length = length_marker - k + 4
            if uniq_length >= 3:
                #print(uniq_length)
                marker_start = round(y[0]+k-4)
                marker_end = round(y[1])
                marker_seq = seqs[x][marker_start:marker_end]

                if marker_start-10 < 0:
                    print('Warning! Limited context available, markers near end of provided reference')
                    context_seq = seqs[x][0:marker_end+10]

                elif marker_end+10 > len(seqs[x]):
                    print('Warning! Limited context available, markers near end of provided reference')
                    context_seq = seqs[x][marker_start-10:]

                else:
                    context_seq = seqs[x][marker_start-10:marker_end+10]
                #print('Adding ' + context_seq + 'here')
                best_markers[x].append([(marker_start, marker_end), uniq_length, marker_seq, context_seq])
    
    for item in best_markers:
        unsorted_markers_list = best_markers[item]
        length_sorted_markers = sorted(unsorted_markers_list, key = itemgetter(2))
        
        if len(length_sorted_markers) >= max_regions:
            best_markers[item] = length_sorted_markers[:max_regions]
        
        else:
            best_markers[item] = length_sorted_markers
    
    return best_markers

def align_reads(fastq_reads, fasta_ref, k,  match, mismatch, gap_open, gap_extend,
                context_score, fine_alignment_score, max_regions):
    '''
    Called by write bins. Calls the define_aligner and define_markers. Gets the regions to score across and 
    uses those for Biopython PairwiseAligner() alignment to all of the provided reads.
    
    Returns dictionary with plasmids (references) as keys and reads as values
    
    If there is only one plasmid, (or if the length of unique sequence is too long), return None to run medaka.
    '''
    reads_dict = defaultdict(list)
    aligner = define_aligner(match, mismatch, gap_open, gap_extend)
    
    #Get the context and unique sequence from define makers
    markers_dict = define_markers(fasta_ref, k, max_regions)
    
    #Leave the biobin steps here if there is only 1 plasmid
    if len(markers_dict) == 1:
        print('Single plasmid found, passing all reads and reference to medaka')
        return None
    
    for plasmid in markers_dict:
        for marker_region in markers_dict[plasmid]:
            context_fwd_marker = marker_region[3]
            context_rev_marker = reverse_complement(str(context_fwd_marker))
            #Match with context
        
            best_score = aligner.score(context_fwd_marker, context_fwd_marker)
            #Iterate through each plasmid, then each read
            for record in SeqIO.parse(fastq_reads, 'fastq'):
                #Determines if forward or backward is better without doing alignment
                scores = [aligner.score(record.seq, context_fwd_marker),
                         aligner.score(record.seq, context_rev_marker)]
                
                top_score = (max(scores), scores.index(max(scores)))

                if top_score[0] >= float(context_score)*best_score:
                    #Define top possible score for unique vs unique
                    best_fine_score = aligner.score(marker_region[2], marker_region[2])

                    #If the forward marker had the best score
                    if top_score[1] == 0:
                        top_fwd_alignment = sorted(aligner.align(record.seq, context_fwd_marker))[0]
                        #Get the sequence that the aligned to the larger region
                        subseq = record.seq[top_fwd_alignment.aligned[0][0][0]:top_fwd_alignment.aligned[0][-1][-1]]
                        #Possible add length subseq check
                        if aligner.score(subseq, marker_region[2]) >= float(fine_alignment_score)*best_fine_score:
                            reads_dict[plasmid].append(record)

                    #If the reverse marker has the best score
                    if top_score[1] == 1:
                        #Get the sequence that the aligned to the larger region
                        top_rev_alignment = sorted(aligner.align(record.seq, context_fwd_marker))[0]
                        subseq = record.seq[top_rev_alignment.aligned[0][0][0]:top_rev_alignment.aligned[0][-1][-1]]
                        if aligner.score(subseq, reverse_complement(str(marker_region[2]))) >= float(fine_alignment_score)*best_fine_score:
                            reads_dict[plasmid].append(record)
                
    return reads_dict

def write_bins(fastq_reads, reference, k, output_directory, match = 3,
               mismatch = -6, gap_open = -10, gap_extend = -5,
               context_map = 0.80, fine_map = 0.95, max_regions = 3):
    
    reads_dict = align_reads(fastq_reads, reference, k,  match, mismatch, gap_open, gap_extend, context_map, fine_map, max_regions)
    
    if reads_dict == None:
        return None
    
#     output_directory = output_directory.strip('/')
    output_bins = []
    for plasmid in reads_dict:
        output_bins.append(output_directory+'/%s' % plasmid)
        with open(output_directory+'/%s.fastq' % plasmid, 'w') as out_bin:
            SeqIO.write(reads_dict[plasmid], out_bin, "fastq")
            
        out_bin.close()

    return output_bins

def run_medaka_binned(fastq_reads, reference, k, output_directory, args,
                      match = 3, mismatch = -6, gap_open = -10, gap_extend = -5,
                      context_map = 0.80, fine_map = 0.95, max_regions = 3):
    
    output_bins = write_bins(fastq_reads, reference, k, output_directory, match, mismatch, gap_open, gap_extend,
                            context_map, fine_map, max_regions)
    
    if output_bins == None:
        #If the alignment step was broken for single plasmid or (long unique seq), just pass the reads and referense to medaka
        medaka_wrap.runMedaka(fastq_reads,
                              reference,
                              output_directory,
                              args.threads,
                              args.model,
                              args.igv)
        return
    
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
                                      args.model, args.igv)
    return