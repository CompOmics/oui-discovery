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

DS=PXD002057
#DS=PXD005833
#DS=PXD014258

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

#run PeptideProphet --- should I put -d?

for file in $OUD/$DS/${DS}-$DBP/*.pep.xml
do
	filename="${file##*/}"
	BASE="${filename%.pep.xml}"
	/usr/local/tpp/bin/xinteract -D$DBD/$DB/$DBF -N${BASE}_${DBP}_PeptideProphet.pep.xml $OUD/$DS/${DS}-$DBP/${BASE}.pep.xml

done

#run ProteinProphet

for file in $OUD/$DS/${DS}-$DBP/*_${DBP}_PeptideProphet.pep.xml
do
	filename="${file##*/}"
	BASE="${filename%_${DBP}_PeptideProphet.pep.xml}"
	/usr/local/tpp/bin/ProteinProphet ${BASE}_${DBP}_PeptideProphet.pep.xml ${BASE}_${DBP}_ProteinProphet.pep.xml
done

cd $OUD/
