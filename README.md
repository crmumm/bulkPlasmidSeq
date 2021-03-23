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
    
______________________________________________________________________
BulkPlasmidSeq usage examples:
    
    For binning sequences based on reference - biobin:
  
        python bulkPlasmidSeq.py biobin -i path/to/reads -r path/to/plasmids -o output_directory
    
    For generating concensus sequences - Medaka:
        
        python bulkPlasmidSeq.py Medaka -i my_reads.fastq -r my_plasmid_genome.fa -o output_directory -t 4
        
    For unfiltered reads - Guppy filters to Q7 ~85% basecalling accuracy:
    
        python bulkPlasmidSeq.py --filter -i unfiltered_reads.fastq -r my_plasmid_genome.fa -q 8 

______________________________________________________________________

Troubleshooting:

     Common Problems:
	
	Taking screenshots with IGV crashes or gives a Null Pointer error.
	
		Make sure a .genome file has been created for this 'Plasmid Genome' and that
		this genome is available in the genome dropdown list. This pipeline does not
		create it the .genome file. To create this file open igv, go to Genomes
		dropdown menu and use Create .genomes. Give the genome a name and point it to
		the plasmid_genome.fasta created by the pileline or your reference fasta.  
