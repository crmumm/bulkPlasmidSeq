#!bin/bash/
READS=
REFERENCE=
OUTPUTDIR=
THREADS=

END_SIZE=
BARCODE_THRESHOLD=

PATH_TO_IGV=

#For Porechop with no filtering or igv screenshots
python bulkPlasmidSeq -i ${READS} -BC ${REFERENCE} -o ${OUTPUTDIR} -t ${THREADS} -P

#For Porechop with filtering and screenshots
python bulkPlasmidSeq -i ${READS} -BC ${REFERENCE} -o ${OUTPUTDIR} -t ${THREADS} -PF \
	--min_length ${MIN} --max_length ${MAX} -q 8 \
	--igv ${PATH_TO_IGV} --screenshot

#___Don't forget to activate medaka environment for using Medaka__

#For Medaka
python bulkPlasmidSeq -i ${READS} -r ${REFERENCE} -o ${OUTPITDIR} -t ${THREADS} -M


