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
from collections import OrderedDict 
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
from tensorflow.keras import models
from scipy import stats
from Bio import SeqIO


np.set_printoptions(threshold=np.inf)

#length of input seq
eLENGTH=30
#depth of onehot encoding
eDEPTH=4

def read_data_files(fn, crisproff):
    d = OrderedDict()
    fasta = SeqIO.parse(fn, "fasta")
    for r in fasta:
        id, seq = r.id, r.seq
        if id in d:
            print('sequence identifiers most be unique', file=sys.stderr)
            sys.exit(1)
        elif not len(seq) == eLENGTH:
            print('"%s" is not a sequence of length %d' % (s, eLENGTH), file=sys.stderr)
            sys.exit(1)
        else:
            d[id] = list()
            d[id].append(seq)
    with open(crisproff, 'rt') as f:
        head = f.readline().rstrip().split('\t')
        for i,l in enumerate(f):
            v = l.rstrip().split('\t')
            (id, g) = (v[0], v[6])
            try:
                g = float(g)
            except:
                print('not a float in line %d', (s, i+1), file=sys.stderr)
                sys.exit(1)
            d[id].insert(0, g)
    return d

def onehot(x):
    z = list()
    for y in list(x):
        if y in "Aa":  z.append(0)
        elif y in "Cc": z.append(1)
        elif y in "Gg": z.append(2)
        elif y in "TtUu": z.append(3)
        else:
            print("Non-ATGCU character " + data[l], file=sys.stderr)
            raise Exception
    return z

def set_data(DX, s):
    for j,x in enumerate(onehot(s)):
        DX[j][x] = 1


def preprocess_seq(data):
    
    dkeys = data.keys()
    DATA_X = np.zeros((len(dkeys),eLENGTH,eDEPTH), dtype=np.float32) # onehot
    DATA_G = np.zeros((len(dkeys)), dtype=np.float32) #deltaGb

    seqs = list()


    for l, id in enumerate(dkeys):
        d = data[id]
        set_data(DATA_X[l], d[1])
        DATA_G[l] = d[0]
        seqs.append(d[1])

    return (seqs, DATA_X, DATA_G,)


OUTD = sys.argv[1]
fasta = sys.argv[2]
crisproff = sys.argv[3]
MODELD=sys.argv[4:] #paths to model(s)


d = read_data_files(fasta, crisproff)
(s, x, g) = preprocess_seq(d)
tinput = list()
tinput.append(x)
tinput.append(g)

#evaluate per model
v = list()
for i, md in enumerate(MODELD):
    model = models.load_model(md)
    model.compile(loss='mse', metrics=['mae', 'mse'])
    predicted = model.predict(tinput, use_multiprocessing=True, workers=32)
    predicted = [p[0] for p in predicted]
    v.append([list(d.keys()), s, predicted])

#average model outputs
ids = v[0][0]
sequences = v[0][1]
predicted = [0 for i in range(len(sequences))]
nmodels = len(MODELD)
for m in range(nmodels):
    predicted = [predicted[i] + v[m][2][i]/nmodels for i in range(len(sequences))]
with open(OUTD + '/crispron.csv', 'wt') as f:
    print(','.join(['ID', '30mer', 'CRISPRon']), file=f)
    for si, ss, sv in zip(ids, sequences, predicted):
        print(','.join([si, str(ss), '%.2f' % sv]), file=f)
