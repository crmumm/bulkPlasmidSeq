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
        
    Porechop
    
        Porechop is officially unsupported by Nanopore, but is one of the better tools if using custom barcodes. 
        For custom barcoding, using sequence features other than the ligated ONT barcodes, please use crmumm's 
        fork of Porechop. This integrates barcoding.py into Porechop. 
    
            git clone https://github.com/rrwick/Porechop.git
            cd Porechop
            python3 setup.py install
            
        Source: https://github.com/rrwick/Porechop
            
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
    
    
    ______________________________________________________________________
    BulkPlasmidSeq usage examples:
    
    For binning sequences based on reference - Porechop:
  
        python bulkPlasmidSeq.py Porechop -i path/to/reads -r path/to/plasmids -o output_directory
        
        python bulkPlasmidSeq.py Porechop -i path/to/reads -r path/to/plasmids -o output_directory --screenshot \
            --igv path/to/igv.sh
            
    
    For generating concensus sequences - Medaka:
        
        python bulkPlasmidSeq.py Medaka -i my_reads.fastq -r my_plasmid_genome.fa -o output_directory -t 4
        
    For unfiltered reads - Guppy filters to Q7 ~85% basecalling accuracy:
    
        python bulkPlasmidSeq.py -i my_reads.fastq -r my_plasmid_genome.fa -q 8 
