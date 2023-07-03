# Unlike the rest of the pipeline, this should not be run in the seattle-flu Conda environment.

PICARD_PATH="/data/mambaforge/pkgs/picard-3.0.0-hdfd78af_0/share/picard-3.0.0-0/picard.jar"
WORK_DIR="work/assembly"
set -x

for REFERENCE in `ls "$WORK_DIR"/process/mapped`
do \
  mkdir -p "$WORK_DIR"/summary/picard/$REFERENCE
  # for SAMPLE in {list of sample IDs}
  for BAM in `ls "$WORK_DIR"/process/mapped/$REFERENCE/*.sorted.bam`
  do \
    echo $BAM
    SAMPLE=`basename "$BAM" .sorted.bam`
    echo $SAMPLE
    if test -f "$WORK_DIR"/summary/picard/$REFERENCE/$SAMPLE.picard.txt
      then echo "FILE EXISTS"
      else java -jar "$PICARD_PATH" CollectWgsMetrics \
        REFERENCE_SEQUENCE=./references/criseq/$REFERENCE.fasta \
        INPUT="$WORK_DIR"/process/mapped/$REFERENCE/$SAMPLE.sorted.bam \
        OUTPUT="$WORK_DIR"/summary/picard/$REFERENCE/$SAMPLE.picard.txt
      fi
  done
done

