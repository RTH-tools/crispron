# CRISPRon
CRISPRon: enhanced data-driven on-target CRISPR-Cas9 gRNA efficiency prediction

## Webserver
If you are just looking for ways to design or analyse your CRISPR sgRNAs, you may find our webserver at https://rth.dk/resources/crispr/crispron/ an easier starting point.

## Installation

### Prerequisites
The software has been tested with the following versions, but later versions of
CRISPRoff, python, biopthon, and tensorflow should also work.

* biopython  : 1.78 
* python     : 3.9.2
* tensorflow : 2.4.1
* viennarna  : 2.2.5
* CRISPRoff  : 1.1.2

### CRISPRoff
You need to install CRISPRoff (version 1.1.2 or later) which can be downloaded from 

https://rth.dk/resources/crispr/crisproff/download

After downloading and un-packing, you should move
CRISPRspec_CRISPRoff_pipeline.py to the bin folder of CRISPRon and
energy_dics.pkl to data/model/energy_dics.pkl

### Other requirements

The easiest way to get the needed prerequisites to run CRISPRon is through
conda. If you have conda installed already you can skip this step, otherwise go
to https://docs.conda.io/en/latest/miniconda.html to learn how to install conda
on your system. Once conda is correctly installed. You need to install the
crispron requirements with

	conda create -y -c bioconda -c conda-forge -n crispron python=3 tensorflow=2
	conda install -y -c bioconda -c conda-forge -n crispron biopython viennarna=2.2.5

Later versions are also expected to work. However, the program depends on
RNAfold and versions other than 2.2.5 of the ViennaRNA package will give
slightly different results.

## Testing the software
Assuming you have installed the prerequisites in a conda environment called
crispron, you can run the built-in software test

	conda activate crispron
	./bin/test.sh

Which should end with

	TEST ok

Note that the test requires diff from diffutils which must be installed as a
normal package on your system.

## Running the software
Assuming you have installed the prerequisites in a conda environment called
crispron, you can run the software on the test data

	conda activate crispron
	./bin/CRISPRon.sh test/seq.fa test/outdir


**Output files**

The script CRISPRon.sh outputs the following files:
- 23mers.fa: the 23 nt target + PAM sequences extracted from the input FASTA file
- 30mers.fa: the 30 nt prefix + target + PAM + suffix sequences extracted from the input FASTA file
- CRISPRparams.tsv: a tab-separated table containing the free energy changes computed by the CRISPRoff software. These are the RNA-DNA hybridisation energy (both unweighted and weighted), the DNA-DNA opening energy, the RNA-RNA spacer self-folding energy, and the CRISPRoff score. For details, read about CRISPRoff at https://rth.dk/resources/crispr/
- crispron.csv: a comma-separated table containing the 30 nt prefix + target + PAM + suffix sequences  and the CRISPRon predicted indel frequencies


To run the software on your own data, first construct a FASTA file with all the
sequences you want to have tested. See test/seq.fa for FASTA format. Just
remember that the program needs at least 30 nucleotides to fit the full target

	prefix (4nt) -- target (20nt) -- PAM (3nt, NGG) -- suffix (3nt)

And then run the program with your own fasta file and an appropriate output
directory.

### Example run

Running on the sequence the 

	>test
	ACTGAACTTGAAAAGCAAAAAGAAACTGGCCATACTTTCGAAGAAATGCTACTGACTG

You will get the following output in crispron.csv

| ID        | 30mer                                    | CRISPRon  |
| --------- | ---------------------------------------- | ---------:|
| test_p_7  | `TGAACTTGAAAAGCAAAAAGAAAC`**`TGG`**`CCA` | 11.48     |
| test_m_30 | `GTCAGTAGCATTTCTTCGAAAGTA`**`TGG`**`CCA` | 37.06     |

Which means that CRISPRon finds two possible 30mers with targets + PAM sequence
inside. In the first, test_p_7, the target starts on position 7 on the plus
strand. In the second, test_m_30, the target is on the negative strand and the
PAM + target is found on position 30 in the test sequence as the reverse
complement to 

	CCATACTTTCGAAGAAATGCTAC

which is

	GTAGCATTTCTTCGAAAGTATGG

## Data and training

The data used for training CRISPRon v. 1.0 may be downloaded from the official CRISPRon page at (https://rth.dk/resources/crispr/crispron/download)

The DeepCRISPRon_train.py is run like this

	export OPT=adam
	export LEARN=0.0001
	export EPOCHS=3000
	export N_VAL=6
	export N_MOD=5
	export SEQ_C=30mer_gRNA
	export VAL_C=Quant_norm_efficiency
	export VAL_G=CRISPRoff
	export BATCH_SIZE=500
	export PREF="validation_set"
	export DT= DeepCRISPRon_train.py
	export SEED=0
	export TYPE=CG
	export M=2

	python3 $DT $OPT $LEARN $EPOCHS $SEQ_C $VAL_C $VAL_G $N_VAL $M $BATCH_SIZE $SEED $TYPE  $PREF*

With the options above, you would need a partition of the data in 6 datasets,
where the 6th will be held out. The data sets should be named
validation_set1.csv ..  ..6.csv, and the sixth will be held out. With M=2, the
second data set would be used for validation (early stopping).

How to partition the data is described in the paper

Xiang, X., Corsi, G.I., Anthon, C. et al. Enhancing CRISPR-Cas9 gRNA efficiency
prediction by data integration and deep learning. Nat Commun 12, 3238 (2021).
https://doi.org/10.1038/s41467-021-23576-0


## Copyright

Copyright 2021 by the contributors (see AUTHORS file)

This is a free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This software is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with
this software, see LICENSE. If not, see http://www.gnu.org/licenses/.

## Citations

If you use CRISPRon in your publication please cite

**CRISPRon: enhanced data-driven on-target CRISPR-Cas9 gRNA efficiency prediction.** Xiang X, Corsi GI, Anthon C, Qu K, , Pan X, Liang X, Han P, Dong Z, Liu L, Zhong J, Ma T, Wang J, Zhang X, Jiang H, Xu F, Liu X, Xu X, Wang J, Yang H, Bolund L, Church GM, Lin L, Gorodkin J, Luo Y. Submitted.

The data used for training of CRISPRon comes in part from the paper from Kim
2019, please cite that as well

**SpCas9 activity prediction by DeepSpCas9, a deep learning-based model with high generalization performance.** Kim, H.K. et al. Sci Adv 5, eaax9249 (2019).

If you use CRISPRspec / CRISPRoff in your publication please cite

**CRISPR-Cas9 off-targeting assessment with nucleic acid duplex energy parameters.** Alkan F, Wenzel A, Anthon C, Havgaard JH, Gorodkin J Genome Biol.
2018 Oct 26;19(1):177


## Contact

In case of problems or bug reports, please contact <software+crispron@rth.dk>

