# CRISPRon
CRISPRon: enhanced data-driven on-target CRISPR-Cas9 gRNA efficiency prediction

## Requirements

- Docker version 20.10.9
- docker-compose version 1.29.2
- [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)


### Building
To build with GPU run:

```
docker-compose up -d --build
```

## Testing the software
Enter the docker container and run 

```
docker exec -it $(docker-compose ps -q) bash
./bin/test.sh
```

Which should end with

	TEST ok


## Running the software
Assuming you have installed the prerequisites in a conda environment called
crispron, you can run the software on the test data

	conda activate crispron
	./bin/CRISPRon.sh test/seq.fa test/outdir

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

