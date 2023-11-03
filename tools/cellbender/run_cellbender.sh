input=$1
out=$2

/usr/bin/singularity exec --cleanenv --nv /vol/projects/BIIM/resources/tools/singularity/cellbender.sif cellbender remove-background \
  --input $input \
  --output $out \
  --expected-cells 8000 \
  --cuda \
  --cells-posterior-reg-calc 10 \
  --posterior-batch-size 2 \
  --epochs 150
