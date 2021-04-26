#!/bin/bash

which diff || exit 1

./bin/CRISPRon.sh test/seq.fa test/outdir.test

if (diff -r test/outdir.test test/outdir.original | grep -e Only -e diff); then
	echo TEST failed
	exit 1
fi

echo TEST ok
