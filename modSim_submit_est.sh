#!/bin/bash

mkdir "$TMPDIR"/modSim/

cd "$HOME"/modSim

for i in `seq 1 150`;
do
	sed s/iter/$i/g modSim_jobs_est.sh > CUR_submit.sh
    qsub CUR_submit.sh
done

rm -r "$TMPDIR"/modSim/
rm -r CUR_submit.sh



