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

#DS=PXD002057
DS=PXD005833

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

#run PeptideProphet with nonspecific enzyme and no mapping to database

for file in $OUD/$DS/${DS}-$DBP/*.pep.xml
do
	filename="${file##*/}"
	BASE="${filename%.pep.xml}"
	/usr/local/tpp/bin/xinteract -eN -nR -D$DBD/$DB/$DBF -N${BASE}_${DBP}_PeptideProphet.pep.xml $OUD/$DS/${DS}-$DBP/${BASE}.pep.xml

done

cd $WD
