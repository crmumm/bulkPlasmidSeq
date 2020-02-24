#!bin/bash/
READS=
REFERENCE=
OUTPUTDIR=
THREADS=

END_SIZE=
BARCODE_THRESHOLD=

PATH_TO_IGV_SH=

#For Porechop with no filtering or igv screenshots
python bulkPlasmidSeq Porechop -i ${READS} -BC ${REFERENCE} -o ${OUTPUTDIR} -t ${THREADS}

#For Porechop with filtering and screenshots
python bulkPlasmidSeq Porechop -i ${READS} -BC ${REFERENCE} -o ${OUTPUTDIR} -t ${THREADS} --filter \
	--min_length ${MIN} --max_length ${MAX} -q 7 \
	--igv ${PATH_TO_IGV} --screenshot

#Don't forget to activate medaka environment for using Medaka

#For Medaka
python bulkPlasmidSeq Medaka -i ${READS} -r ${REFERENCE} -o ${OUTPITDIR} -t ${THREADS}

