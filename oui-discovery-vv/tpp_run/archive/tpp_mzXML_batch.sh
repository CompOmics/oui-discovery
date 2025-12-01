#!bin/bash

#general variabls
WD=/home/vvshazia/tpp_pride_reanalysis
FD=/mnt/d/vasv1101/Documents_discC/msfragger_pride_reanalysis/mzXML
DBD=/mnt/d/vasv1101/Documents_discC/msfragger_pride_reanalysis/concatenated-fastas
#DBF=concat.fasta
DBF=concat_maped.fasta
OUD=$WD/tpp_pride_reanalysis_arch230525

#database variables

#DB=Human_2023_01_canonical_ionbot_db
#DBP=canon
DB=Human_2023_01_trembl_ionbot_db
DBP=trembl
#DB=20241102-openprot-database_ionbot_db
#DBP=openprot

#dataset variabls

DS=PXD014258

#run Comet

for file in $FD/${DS}_mzXML/*.mzXML
do
	/usr/local/tpp/bin/comet -P$WD/comet.params.new -D$DBD/$DB/$DBF ${file}
done

#move files to separate directory

mkdir $OUD/$DS/
mkdir $OUD/$DS/$DS-$DBP
mv $FD/${DS}_mzXML/*.pep.xml* $OUD/$DS/${DS}-$DBP/
cd $OUD/$DS/${DS}-$DBP

#batches variables

BATCH1=ESC-HF-Sample-BT474
BATCH2=ESC-HF-Sample-MCF
BATCH3=ESC-HF-SampleHela

#run PeptideProphet

/usr/local/tpp/bin/xinteract -D$DBD/$DB/$DBF -N${BATCH1}_${DBP}_PeptideProphet.pep.xml $OUD/$DS/${DS}-$DBP/*${BATCH1}*.pep.xml
/usr/local/tpp/bin/xinteract -D$DBD/$DB/$DBF -N${BATCH2}_${DBP}_PeptideProphet.pep.xml $OUD/$DS/${DS}-$DBP/*${BATCH2}*.pep.xml
/usr/local/tpp/bin/xinteract -D$DBD/$DB/$DBF -N${BATCH3}_${DBP}_PeptideProphet.pep.xml $OUD/$DS/${DS}-$DBP/*${BATCH3}*.pep.xml


#run ProteinProphet

/usr/local/tpp/bin/ProteinProphet ${BATCH1}_${DBP}_PeptideProphet.pep.xml ${BATCH1}_${DBP}_ProteinProphet.pep.xml
/usr/local/tpp/bin/ProteinProphet ${BATCH2}_${DBP}_PeptideProphet.pep.xml ${BATCH2}_${DBP}_ProteinProphet.pep.xml
/usr/local/tpp/bin/ProteinProphet ${BATCH3}_${DBP}_PeptideProphet.pep.xml ${BATCH3}_${DBP}_ProteinProphet.pep.xml


cd $OUD/
