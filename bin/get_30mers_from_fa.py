#!/usr/bin/env python3
########################################################################
#
#  Copyright (c) 2021 by the contributors (see AUTHORS file)
#
#  This file is part of the CRISPRon
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  It is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this script, see file COPYING.
#  If not, see <http://www.gnu.org/licenses/>.
##########################################################################
import sys
import os
from subprocess import *
from Bio import SeqIO
import argparse as ap
import io
from pathlib import Path

REV_NT_MAP = {'-':'','a':'t','A':'T','c':'g','C':'G','g':'c','G':'C', 't':'a','T':'A','u':'a','U':'A','n':'n','N':'N'}

PRE_GUIDE=4
GUIDE=20
PRE_PAM=1
PAM='GG'
nPAM=len(PAM)
POST_PAM=3
TOTAL=PRE_GUIDE+GUIDE+PRE_PAM+nPAM+POST_PAM
PAMoffset=PRE_GUIDE+GUIDE+PRE_PAM


def rev_comp_seq(seq):
    seqx = [REV_NT_MAP[s] for s in seq]
    seqx.reverse()
    return ''.join(seqx)

#check if file is fasta
def is_fasta(fasta_s):
    fasta = SeqIO.parse(io.StringIO(fasta_s), "fasta")
    return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

# Read 30mer sequences from fasta file
def construct_30mers_fasta(fasta_s):
    ontargets={}
    for record in SeqIO.parse(io.StringIO(fasta_s), "fasta"):
        seq = str(record.seq)
        if len(str(record.seq))<TOTAL:
            print("seq is %d and too short, at least 30nt needed" % len(seq), file=sys.stderr)
            sys.exit(1)
        else:
            # get valid TARGETS + context
            x=0
            grna_seq=""
            for i in range(len(seq)):
                if seq[i].upper() not in ["A","C","G","T","U"]:
                    grna_seq =""
                elif len(grna_seq)<TOTAL-1:
                    grna_seq = grna_seq + seq[i]
                else:
                    grna_seq = grna_seq + seq[i]
                    if grna_seq[PAMoffset:PAMoffset+nPAM].upper() == PAM:
                        ontargets[record.id+"_p_"+str(i-(TOTAL-PRE_GUIDE-2))] = (i-TOTAL+1, grna_seq)
                    grna_seq = grna_seq[1:]

            grna_seq=""
            revcseq = rev_comp_seq(seq)
            for i in range(len(revcseq)):
                if revcseq[i].upper() not in ["A","C","G","T","U"]:
                    grna_seq =""
                elif len(grna_seq)<TOTAL-1:
                    grna_seq = grna_seq + revcseq[i]
                else:
                    grna_seq = grna_seq + revcseq[i]
                    if grna_seq[PAMoffset:PAMoffset+len(PAM)].upper() == PAM:
                        ontargets[record.id+"_m_"+str(len(revcseq) + POST_PAM - i)] = (i-TOTAL+1, grna_seq)
                    grna_seq = grna_seq[1:]
    return ontargets

if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument('-f', help='Path to fasta input', type=str, required=True, default='')
    parser.add_argument('-m', help='Path to output 30mers in fasta format', type=str, required=False, default='30mers.fa')
    parser.add_argument('-g', help='Path to output guides + PAM', type=str, required=False, default='23mers.fa')
    args = parser.parse_args()

    if args.f and os.path.exists(args.f):
        seqfasta = Path(args.f).read_text()
    else:
        print('Needs a fasta file (-f) as input')
        sys.exit(1)

    if is_fasta(seqfasta):
        Seqs = construct_30mers_fasta(seqfasta)
    else:
        print('Provided input does not appear to be fasta')
        sys.exit(1)

    with open(args.m,"wt") as out_mk, open(args.g,"wt") as out_gk:
        for k,v in Seqs.items():
            out_mk.write(">"+k+"\n"+v[1]+"\n")
            out_gk.write(">"+k+"\n"+v[1][PRE_GUIDE:TOTAL-POST_PAM]+"\n")
