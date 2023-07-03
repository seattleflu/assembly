
WORK_DIR=work/assembly
CONSENSUS_DIR="${WORK_DIR}/consensus_genomes"
for REFERENCE in `ls $CONSENSUS_DIR`
  do echo $REFERENCE
  for FILE in `ls $CONSENSUS_DIR/$REFERENCE/*.masked_consensus.fasta`
    do SAMPLE=`basename $FILE .masked_consensus.fasta`
    echo $SAMPLE
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$CONSENSUS_DIR/$REFERENCE/$SAMPLE.masked_consensus.fasta" > "$CONSENSUS_DIR/$REFERENCE/$SAMPLE.masked_consensus.single_line.fasta"
    done
  for FILE in `ls $CONSENSUS_DIR/$REFERENCE/*.masked_consensus.single_line.fasta`
    do (cat "$FILE"; echo) >> $CONSENSUS_DIR/$REFERENCE.all-samples.masked_consensus.fasta
    done
  done
