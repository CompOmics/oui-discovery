#!bin/bash

#general variabls
WD=/home/vvshazia/tpp_pride_reanalysis
FD=/mnt/d/vasv1101/Documents_discC/msfragger_pride_reanalysis/mzXML
DBD=/mnt/d/vasv1101/Documents_discC/pride_reanalysis/entrapment/entrapment-concat-fasta/entrapment-concat-fasta
DBF=concat.fasta
OUD=$WD/tpp_pride_reanalysis_entrp_arch040825

#database variables

#DB=Human_2023_01_canonical_ionbot_db
#DBP=canon
#DB=Human_2023_01_trembl_ionbot_db
#DBP=trembl
DB=peptide_openprot_ionbot_db
DBP=openprot

#dataset variabls

DS=PXD014258

#run Comet

for file in $FD/${DS}_mzXML/*.mzXML
do
	/usr/local/tpp/bin/comet -P$WD/comet.params_entrapment.new -D$DBD/$DB/$DBF ${file}
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

#run PeptideProphet with nonspecific enzyme and no mapping to database

/usr/local/tpp/bin/xinteract -eN -nR -D$DBD/$DB/$DBF -N${BATCH1}_${DBP}_PeptideProphet.pep.xml $OUD/$DS/${DS}-$DBP/*${BATCH1}*.pep.xml
/usr/local/tpp/bin/xinteract -eN -nR -D$DBD/$DB/$DBF -N${BATCH2}_${DBP}_PeptideProphet.pep.xml $OUD/$DS/${DS}-$DBP/*${BATCH2}*.pep.xml
/usr/local/tpp/bin/xinteract -eN -nR -D$DBD/$DB/$DBF -N${BATCH3}_${DBP}_PeptideProphet.pep.xml $OUD/$DS/${DS}-$DBP/*${BATCH3}*.pep.xml

cd $OUD/

