#!bin/bash

#general variabls
PHILOSOPHER=/mnt/d/vasv1101/Documents_discC/FragPipe-23.0/tools/Philosopher/philosopher-v5.1.1.exe
WD=/mnt/d/vasv1101/Documents_discC/msfragger_pride_reanalysis/msfragger_pride_reanalysis_entrapment/msfragger_pepXML
DBD=/mnt/d/vasv1101/Documents_discC/pride_reanalysis/entrapment/entrapment-concat-fasta/entrapment-concat-fasta
DBF=concat.fasta

#database variables

DB=Human_2023_01_canonical_ionbot_db
DBP=canon
#DB=Human_2023_01_trembl_ionbot_db
#DBP=trembl
#DB=peptide_openprot_ionbot_db
#DBP=openprot

#dataset variabls

#DS=PXD002057
DS=PXD005833

#go to workspace

cd $WD/$DS/${DS}-$DBP

#copy database to workspace

cp $DBD/$DB/$DBF .

#init workspace and annotate database

$PHILOSOPHER workspace --init
$PHILOSOPHER database --annotate $DBF --prefix decoy

#run PeptideProphet

for file in $WD/$DS/${DS}-$DBP/*.pepXML
do
    echo $file
    win_file=$(wslpath -w "$file")
    echo $win_file
    $PHILOSOPHER peptideprophet \
        --nonparam --expectscore --decoyprobs --masswidth 1000.0 \
        --clevel -0 --enzyme nonspecific --decoy decoy --database "$DBF" \
        --output $DBP \
        "$win_file"
done

#remove database from workspace

rm $DBF

#go back to root dir

cd $WD
