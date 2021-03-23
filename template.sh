#!bin/bash/
READS=
REFERENCE=
OUTPUTDIR=
THREADS=

END_SIZE=
BARCODE_THRESHOLD=

PATH_TO_IGV_SH=

#Don't forget to activate medaka environment for using Medaka

#For Medaka
python bulkPlasmidSeq Medaka -i ${READS} -r ${REFERENCE} -o ${OUTPUTDIR} -t ${THREADS}

