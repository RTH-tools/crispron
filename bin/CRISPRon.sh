#!/bin/bash
BINDIR=`dirname $0`
DATADIR=$BINDIR/../data
CRISPRON_FASTA=$1
OUTDIR=$2

if [[ -z $OUTDIR ]]; then
	echo "Needs two arguments to run e.g. $0 test/seq.fa test/outdir" 1>&2
	exit 1
fi

mkdir -p $OUTDIR || exit 1

if [[ ! -s $CRISPRON_FASTA ]]; then
	echo "Needs a fasta file as first input paramete to run e.g. $0 indir/test.fa outdir" 1>&2
	exit 1
fi

which python3 || exit 1
RNAfold --version || exit 1

$BINDIR/get_30mers_from_fa.py  -f ${CRISPRON_FASTA} -m $OUTDIR/30mers.fa -g $OUTDIR/23mers.fa || exit 1

if [[ ! -s $OUTDIR/30mers.fa ]]; then
	echo failed to get target sequences 1>&2
	exit 1
fi

echo "#Running CRISPROff pipeline"
$BINDIR/CRISPRspec_CRISPRoff_pipeline.py \
	--guides $OUTDIR/23mers.fa \
	--specificity_report $OUTDIR/CRISPRspec.tsv \
	--guide_params_out $OUTDIR/CRISPRparams.tsv \
	--duplex_energy_params $DATADIR/model/energy_dics.pkl \
	--no_azimuth || exit 1

echo "#Running CRISPRon evaluation"
$BINDIR/DeepCRISPRon_eval.py $OUTDIR $OUTDIR/30mers.fa $OUTDIR/CRISPRparams.tsv $DATADIR/deep_models/best/*/  2>&1 | grep -v tracing

echo "#All done" 1>&2
