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
import shutil
import traceback
from collections import OrderedDict 
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
np.set_printoptions(threshold=np.inf)

OPT = sys.argv[1] #optimizer (eg. adam)
LEARN = float(sys.argv[2]) #learning rage eg.0.001
EPOCHS = int(sys.argv[3]) # N. of epochs, eg 100
SEQ_C = sys.argv[4] # get seq from this column
VAL_C = sys.argv[5] # get value to train from this column 
VAL_G = sys.argv[6] # get value to train from this column 
N_VAL = int(sys.argv[7]) # number of validation set
TEST_N = int(sys.argv[8]) # skip this test set
BATCH_SIZE=int(sys.argv[9]) # batch test 
SEED = int(sys.argv[10]) # common random seed (0 equals randomly set)
TYPE = sys.argv[11] # common random seed (0 equals randomly set)
TRAIN_FP = sys.argv[12:] # path to train sets

print(
'OPT=%s' % OPT,
'LEARN=%f' % LEARN,
'EPOCHS=%i' % EPOCHS,
'SEQ_C=%s' % SEQ_C,
'VAL_C=%s' % VAL_C,
'VAL_G=%s' % VAL_G,
'N_VAL=%i' % N_VAL,
'TEST_N=%i' % TEST_N,
'BATCH_SIZE=%i' % BATCH_SIZE,
'SEED=%i' % SEED,
'TYPE=%s' % TYPE,
)
print(', '.join(TRAIN_FP))

#length of input seq
eLENGTH=30
#depth of onehot encoding
eDEPTH=4


def read_seq_files(fns, seq_c0, val_c0, val_g0, n_val, test_n, test=False):
    d = OrderedDict()
    for j, fn in enumerate(fns):
        if j+1 == n_val:
            #independent validation set
            continue
        if test and j+1 != test_n:
            #not internal test set, skip when not loading test data
            continue
        if (not test) and j+1 == test_n:
            #internal test set, skip when not loading training data
            continue 
        with open(fn, 'rt') as f:
            head = f.readline().rstrip().split('\t')
            if seq_c0.isnumeric():
                seq_c = int(seq_c0) -1
            else:
                if seq_c0 in head:
                    seq_c = head.index(seq_c0)
                else:
                    raise Exception

            if val_c0.isnumeric():
                val_c = int(val_c0) -1
            else:
                if val_c0 in head:
                    val_c = head.index(val_c0)
                else:
                    raise Exception

            if val_g0.isnumeric():
                val_g = int(val_g0) -1
            else:
                if val_g0 in head:
                    val_g = head.index(val_g0)
                else:
                    raise Exception
            f.seek(0)

            for i,l in enumerate(f):
                v = l.rstrip().split('\t')
                s = v[seq_c]
                if not len(s) == eLENGTH:
                    print('"%s" is not a string of length %d in line %d' % (s, eLENGTH, i+1))
                    continue

                if s in d:
                    print('"%s" is not unique in line %d', (s, i+1))
                    continue

                try:
                    e = float(v[val_c])
                except:
                    print('no float value for "%s" in line %d', (s, i+1))
                    e = 0.0
                    continue
                try:
                    g = float(v[val_g])
                except:
                    print('no float value for "%s" in line %d', (s, i+1))
                    g = 0.0
                    continue

                d[s] = [g, e]
    return d

def onehot(x):
    z = list()
    for y in list(x):
        if y in "Aa":  z.append(0)
        elif y in "Cc": z.append(1)
        elif y in "Gg": z.append(2)
        elif y in "TtUu": z.append(3)
        else:
            print("Non-ATGCU character " + data[l])
            raise Exception
    return z

def set_data(DX, s):
    for j,x in enumerate(onehot(s)):
        DX[j][x] = 1


def preprocess_seq(data):
    
    dkeys = data.keys()
    DATA_X = np.zeros((len(dkeys),eLENGTH,eDEPTH), dtype=np.float32) # onehot
    DATA_G = np.zeros((len(dkeys)), dtype=np.float32) #deltaGb

    DATA_Y = np.zeros((len(dkeys)), dtype=np.float32) # efficiency

    for l, s in enumerate(dkeys):
        d = data[s]
        set_data(DATA_X[l], s)
        DATA_G[l] = d[0]
        DATA_Y[l] = d[1]

    return (list(dkeys), DATA_X, DATA_G, DATA_Y)


import tensorflow as tf
from tensorflow.keras import datasets, layers, models, callbacks, Model, optimizers, Input, utils
from tensorflow.keras.layers import Conv1D, Dropout, AveragePooling1D, Flatten, Dense, concatenate, SpatialDropout1D
from scipy import stats
from random import randint
import sys

# check of inputs and setting a few base parameters
if SEED == 0:
    SEED = randint(0, sys.maxsize)
print('seed:', SEED)
tf.random.set_seed(SEED)
if OPT =='adam':
    optimizer = optimizers.Adam(LEARN)
elif OPT == 'rmsprop':
    optimizer = optimizers.RMSprop(LEARN)
else:
    raise Exception

#inputs
#one hot
inputs = list()
if TYPE.find('C') > -1:
    input_c = Input(shape=(eLENGTH, eDEPTH,), name="input_onehot")
    inputs.append(input_c)

#delta Gb
if TYPE.find('G') > -1:
    input_g = Input(shape=(1,), name="input_dGB")
    inputs.append(input_g)

#
for_dense = list()

#first convolution layer
if TYPE.find('C') > -1:
    conv1_out = Conv1D(100, 3, activation='relu', input_shape=(eLENGTH,4,), name="conv_3")(input_c)
    conv1_dropout_out = Dropout(0.3, name="drop_3")(conv1_out)
    conv1_pool_out = AveragePooling1D(2, padding='SAME', name="pool_3")(conv1_dropout_out)
    conv1_flatten_out = Flatten(name="flatten_3")(conv1_pool_out)
    for_dense.append(conv1_flatten_out)

#second convolution layer
    conv2_out = Conv1D(70, 5, activation='relu', input_shape=(eLENGTH,4,), name="conv_5")(input_c)
    conv2_dropout_out = Dropout(0.3, name="drop_5")(conv2_out)
    conv2_pool_out = AveragePooling1D(2, padding='SAME', name="pool_5")(conv2_dropout_out)
    conv2_flatten_out = Flatten(name="flatten_5")(conv2_pool_out)
    for_dense.append(conv2_flatten_out)

#third convolution layer
    conv3_out = Conv1D(40, 7, activation='relu', input_shape=(eLENGTH,4,), name="conv_7")(input_c)
    conv3_dropout_out = Dropout(0.3, name="drop_7")(conv3_out)
    conv3_pool_out = AveragePooling1D(2, padding='SAME', name="pool_7")(conv3_dropout_out)
    conv3_flatten_out = Flatten(name="flatten_7")(conv3_pool_out)
    for_dense.append(conv3_flatten_out)

#concatenation of conv layers and deltaGb layer
if len(for_dense) == 1:
    concat_out = for_dense[0]
else:
    concat_out = concatenate(for_dense)

for_dense1 = list()

#first dense (fully connected) layer
dense0_out = Dense(80, activation='relu', name="dense_0")(concat_out)
dense0_dropout_out = Dropout(0.3, name="drop_d0")(dense0_out)
for_dense1.append(dense0_dropout_out)


#Gb input used raw
if TYPE.find('G') > -1:
    for_dense1.append(input_g)


if len(for_dense1) == 1:
    concat1_out = for_dense1[0]
else:
    concat1_out = concatenate(for_dense1)


#first dense (fully connected) layer
dense1_out = Dense(80, activation='relu', name="dense_1")(concat1_out)
dense1_dropout_out = Dropout(0.3, name="drop_d1")(dense1_out)

#second dense (fully connected) layer
dense2_out = Dense(60, activation='relu', name="dense_2")(dense1_dropout_out)
dense2_dropout_out = Dropout(0.3, name="drop_d2")(dense2_out)

#output layer
output = Dense(1, name="output")(dense2_dropout_out)

#model construction
model= Model(inputs=inputs, outputs=[output])
model.summary()
model.compile(loss='mse', optimizer=optimizer, metrics=['mae', 'mse'])
utils.plot_model(model, to_file=str(TEST_N) + '.model.png', show_shapes=True, dpi=600)


es = callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=150)
mc = callbacks.ModelCheckpoint(str(TEST_N) + '.model.best', verbose=1, save_best_only=True)

tinput = list()
vinput = list()

print('reading training data')
d = read_seq_files(TRAIN_FP, SEQ_C, VAL_C, VAL_G, N_VAL, TEST_N, test=False)
(s, x, g, y) = preprocess_seq(d)

print('reading validation data')
db = read_seq_files(TRAIN_FP, SEQ_C, VAL_C, VAL_G, N_VAL, TEST_N, test=True)
(sv, xv, gv, yv) = preprocess_seq(db)


if TYPE.find('C') > -1:
    tinput.append(x)
    vinput.append(xv)

if TYPE.find('G') > -1:
    tinput.append(g)
    vinput.append(gv)

tinput = tuple(tinput)
vinput = tuple(vinput)
history = model.fit(tinput, y, validation_data=(vinput, yv), batch_size=BATCH_SIZE, \
        epochs=EPOCHS, use_multiprocessing=True, workers=16, verbose=2, callbacks=[es, mc])
