# CRISPRon
Enhanced on-target CRISPR-Cas9 gRNA efficiency prediction by data generation, integration and deep learning

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

https://rth.dk/resources/crispr/crisproff/downloads

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

To run the software on your own data, first construct a FASTA file with all the
sequences you want to have tested. See test/seq.fa for FASTA format. Just
remember that the program needs at least 30 nucleotides to fit the full target

	prefix (4nt) -- target (20nt) -- PAM (3nt, NGG) -- suffix (4nt)

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

**CRISPR-Cas9 off-targeting assessment with nucleic acid duplex energy
parameters.** Alkan F, Wenzel A, Anthon C, Havgaard JH, Gorodkin J Genome Biol.
2018 Oct 26;19(1):177


## Contact

In case of problems or bug reports, please contact <software+crispron@rth.dk>

