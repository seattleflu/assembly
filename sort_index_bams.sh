WORK_DIR=work/assembly
MAPPED_DIR="${WORK_DIR}/process/mapped"
for REFERENCE in `ls $MAPPED_DIR`
  do echo $REFERENCE
  for FILE in `ls $MAPPED_DIR/$REFERENCE/*.bam`
  #for SAMPLE in {list of sample IDs}
    do SAMPLE=`basename $FILE .bam`
    echo $SAMPLE
    samtools sort $MAPPED_DIR/$REFERENCE/$SAMPLE.bam -o $MAPPED_DIR/$REFERENCE/$SAMPLE.sorted.bam
    samtools index $MAPPED_DIR/$REFERENCE/$SAMPLE.sorted.bam
    done
  done

