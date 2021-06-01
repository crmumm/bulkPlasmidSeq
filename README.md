Dependencies:
    
    All require Python 3, please see software sources for more detailed requirements. I recommend installing
    all of these inside the Medaka conda environment. Samtools and minimap2 already installed with Medaka environment.
    
    Medaka
    
        Prefered installation - Using conda environment, takes care of medaka dependencies:

            conda create -n medaka -c conda-forge -c bioconda medaka
            
            conda activate medaka

            --------------------------------------------------------------

            virtualenv medaka --python=python3 --prompt "(medaka) "
            . medaka/bin/activate
            pip install medaka
            
        Source: https://github.com/nanoporetech/medaka
        
    NanoFilt
    
        One of several available long read filtering softwares, NanoFilt has the advantage of specifying
        the maximimum read length which is good for filtering out extra long reads. 
            
            pip install nanofilt
            pip install nanofilt --upgrade
           
            --------------------------------------------------------------
            
            conda install -c bioconda nanofilt
            
        Source: https://github.com/wdecoster/nanofilt
   
    Samtools
    
    Minimap2
    
    Integrative Genomics Browser - IGV 
    
        Specify --igv path/to/igv.sh to take screenshots of alignments.
        
    Biopython
        
        Used for IO and marker based binning. 
        
        Source: https://biopython.org/wiki/Download
    
    Emboss
        
        Needle pairwise alignments for consensus alignments.
        
        Source: http://emboss.sourceforge.net/download/   
           
______________________________________________________________________
BulkPlasmidSeq usage examples:
    
    For binning sequences based on unique sequences in reference - biobin:
  
        python bulkPlasmidSeq.py biobin -i path/to/reads -r path/to/plasmids -o output_directory
    
    For generating concensus sequences - medaka:
        
        python bulkPlasmidSeq.py medaka -i my_reads.fastq -r my_plasmid_genome.fa -o output_directory -t 4
        
    Filter reads with nanofilt. 
    
        python bulkPlasmidSeq.py medaka --filter \
            -i unfiltered_reads.fastq \
            -r my_plasmid_genome.fa \
            -q 7 \
            -o output_directroy
    
    