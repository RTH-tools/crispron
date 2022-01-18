#!/bin/bash
pip install -r requirements.txt
if [ ! -d ".cache" ]; then
	mkdir .cache
	mkdir -p data/model
	wget -O .cache/crisproff.tar.gz https://rth.dk/resources/crispr/crisproff/downloads/crisproff-1.1.2.tar.gz
	tar -xf .cache/crisproff.tar.gz -C .cache/
	# added CRISPRspec_CRISPRoff_pipeline.py to git. No need to download
	# mv .cache/crisproff-1.1.2/CRISPRspec_CRISPRoff_pipeline.py bin/CRISPRspec_CRISPRoff_pipeline.py
	mv .cache/crisproff-1.1.2/energy_dics.pkl data/model/energy_dics.pkl
fi 

which diff || exit 1

./bin/CRISPRon.sh test/seq.fa test/outdir.test

if (diff -r test/outdir.test test/outdir.original | grep -e Only -e diff); then
	echo TEST failed
	exit 1
fi

echo TEST ok
