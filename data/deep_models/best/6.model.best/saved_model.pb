Ϫ
??
8
Const
output"dtype"
valuetensor"
dtypetype

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype?
?
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ?
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape?"serve*2.2.02unknown8??
z
conv_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*
shared_nameconv_3/kernel
s
!conv_3/kernel/Read/ReadVariableOpReadVariableOpconv_3/kernel*"
_output_shapes
:d*
dtype0
n
conv_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*
shared_nameconv_3/bias
g
conv_3/bias/Read/ReadVariableOpReadVariableOpconv_3/bias*
_output_shapes
:d*
dtype0
z
conv_5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:F*
shared_nameconv_5/kernel
s
!conv_5/kernel/Read/ReadVariableOpReadVariableOpconv_5/kernel*"
_output_shapes
:F*
dtype0
n
conv_5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:F*
shared_nameconv_5/bias
g
conv_5/bias/Read/ReadVariableOpReadVariableOpconv_5/bias*
_output_shapes
:F*
dtype0
z
conv_7/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*
shared_nameconv_7/kernel
s
!conv_7/kernel/Read/ReadVariableOpReadVariableOpconv_7/kernel*"
_output_shapes
:(*
dtype0
n
conv_7/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*
shared_nameconv_7/bias
g
conv_7/bias/Read/ReadVariableOpReadVariableOpconv_7/bias*
_output_shapes
:(*
dtype0
y
dense_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	?P*
shared_namedense_0/kernel
r
"dense_0/kernel/Read/ReadVariableOpReadVariableOpdense_0/kernel*
_output_shapes
:	?P*
dtype0
p
dense_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:P*
shared_namedense_0/bias
i
 dense_0/bias/Read/ReadVariableOpReadVariableOpdense_0/bias*
_output_shapes
:P*
dtype0
x
dense_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:QP*
shared_namedense_1/kernel
q
"dense_1/kernel/Read/ReadVariableOpReadVariableOpdense_1/kernel*
_output_shapes

:QP*
dtype0
p
dense_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:P*
shared_namedense_1/bias
i
 dense_1/bias/Read/ReadVariableOpReadVariableOpdense_1/bias*
_output_shapes
:P*
dtype0
x
dense_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:P<*
shared_namedense_2/kernel
q
"dense_2/kernel/Read/ReadVariableOpReadVariableOpdense_2/kernel*
_output_shapes

:P<*
dtype0
p
dense_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:<*
shared_namedense_2/bias
i
 dense_2/bias/Read/ReadVariableOpReadVariableOpdense_2/bias*
_output_shapes
:<*
dtype0
v
output/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:<*
shared_nameoutput/kernel
o
!output/kernel/Read/ReadVariableOpReadVariableOpoutput/kernel*
_output_shapes

:<*
dtype0
n
output/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameoutput/bias
g
output/bias/Read/ReadVariableOpReadVariableOpoutput/bias*
_output_shapes
:*
dtype0
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
x
Adam/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/learning_rate
q
&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*
_output_shapes
: *
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0
b
total_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_2
[
total_2/Read/ReadVariableOpReadVariableOptotal_2*
_output_shapes
: *
dtype0
b
count_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_2
[
count_2/Read/ReadVariableOpReadVariableOpcount_2*
_output_shapes
: *
dtype0
?
Adam/conv_3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*%
shared_nameAdam/conv_3/kernel/m
?
(Adam/conv_3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv_3/kernel/m*"
_output_shapes
:d*
dtype0
|
Adam/conv_3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*#
shared_nameAdam/conv_3/bias/m
u
&Adam/conv_3/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv_3/bias/m*
_output_shapes
:d*
dtype0
?
Adam/conv_5/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:F*%
shared_nameAdam/conv_5/kernel/m
?
(Adam/conv_5/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv_5/kernel/m*"
_output_shapes
:F*
dtype0
|
Adam/conv_5/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:F*#
shared_nameAdam/conv_5/bias/m
u
&Adam/conv_5/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv_5/bias/m*
_output_shapes
:F*
dtype0
?
Adam/conv_7/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*%
shared_nameAdam/conv_7/kernel/m
?
(Adam/conv_7/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv_7/kernel/m*"
_output_shapes
:(*
dtype0
|
Adam/conv_7/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*#
shared_nameAdam/conv_7/bias/m
u
&Adam/conv_7/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv_7/bias/m*
_output_shapes
:(*
dtype0
?
Adam/dense_0/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	?P*&
shared_nameAdam/dense_0/kernel/m
?
)Adam/dense_0/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_0/kernel/m*
_output_shapes
:	?P*
dtype0
~
Adam/dense_0/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:P*$
shared_nameAdam/dense_0/bias/m
w
'Adam/dense_0/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_0/bias/m*
_output_shapes
:P*
dtype0
?
Adam/dense_1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:QP*&
shared_nameAdam/dense_1/kernel/m

)Adam/dense_1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_1/kernel/m*
_output_shapes

:QP*
dtype0
~
Adam/dense_1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:P*$
shared_nameAdam/dense_1/bias/m
w
'Adam/dense_1/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_1/bias/m*
_output_shapes
:P*
dtype0
?
Adam/dense_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:P<*&
shared_nameAdam/dense_2/kernel/m

)Adam/dense_2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_2/kernel/m*
_output_shapes

:P<*
dtype0
~
Adam/dense_2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:<*$
shared_nameAdam/dense_2/bias/m
w
'Adam/dense_2/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_2/bias/m*
_output_shapes
:<*
dtype0
?
Adam/output/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:<*%
shared_nameAdam/output/kernel/m
}
(Adam/output/kernel/m/Read/ReadVariableOpReadVariableOpAdam/output/kernel/m*
_output_shapes

:<*
dtype0
|
Adam/output/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/output/bias/m
u
&Adam/output/bias/m/Read/ReadVariableOpReadVariableOpAdam/output/bias/m*
_output_shapes
:*
dtype0
?
Adam/conv_3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*%
shared_nameAdam/conv_3/kernel/v
?
(Adam/conv_3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv_3/kernel/v*"
_output_shapes
:d*
dtype0
|
Adam/conv_3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*#
shared_nameAdam/conv_3/bias/v
u
&Adam/conv_3/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv_3/bias/v*
_output_shapes
:d*
dtype0
?
Adam/conv_5/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:F*%
shared_nameAdam/conv_5/kernel/v
?
(Adam/conv_5/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv_5/kernel/v*"
_output_shapes
:F*
dtype0
|
Adam/conv_5/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:F*#
shared_nameAdam/conv_5/bias/v
u
&Adam/conv_5/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv_5/bias/v*
_output_shapes
:F*
dtype0
?
Adam/conv_7/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*%
shared_nameAdam/conv_7/kernel/v
?
(Adam/conv_7/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv_7/kernel/v*"
_output_shapes
:(*
dtype0
|
Adam/conv_7/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*#
shared_nameAdam/conv_7/bias/v
u
&Adam/conv_7/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv_7/bias/v*
_output_shapes
:(*
dtype0
?
Adam/dense_0/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	?P*&
shared_nameAdam/dense_0/kernel/v
?
)Adam/dense_0/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_0/kernel/v*
_output_shapes
:	?P*
dtype0
~
Adam/dense_0/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:P*$
shared_nameAdam/dense_0/bias/v
w
'Adam/dense_0/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_0/bias/v*
_output_shapes
:P*
dtype0
?
Adam/dense_1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:QP*&
shared_nameAdam/dense_1/kernel/v

)Adam/dense_1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_1/kernel/v*
_output_shapes

:QP*
dtype0
~
Adam/dense_1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:P*$
shared_nameAdam/dense_1/bias/v
w
'Adam/dense_1/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_1/bias/v*
_output_shapes
:P*
dtype0
?
Adam/dense_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:P<*&
shared_nameAdam/dense_2/kernel/v

)Adam/dense_2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_2/kernel/v*
_output_shapes

:P<*
dtype0
~
Adam/dense_2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:<*$
shared_nameAdam/dense_2/bias/v
w
'Adam/dense_2/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_2/bias/v*
_output_shapes
:<*
dtype0
?
Adam/output/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:<*%
shared_nameAdam/output/kernel/v
}
(Adam/output/kernel/v/Read/ReadVariableOpReadVariableOpAdam/output/kernel/v*
_output_shapes

:<*
dtype0
|
Adam/output/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/output/bias/v
u
&Adam/output/bias/v/Read/ReadVariableOpReadVariableOpAdam/output/bias/v*
_output_shapes
:*
dtype0

NoOpNoOp
?l
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*?k
value?kB?k B?k
?
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer-4
layer-5
layer-6
layer-7
	layer-8

layer-9
layer-10
layer-11
layer-12
layer-13
layer_with_weights-3
layer-14
layer-15
layer-16
layer-17
layer_with_weights-4
layer-18
layer-19
layer_with_weights-5
layer-20
layer-21
layer_with_weights-6
layer-22
	optimizer
regularization_losses
trainable_variables
	variables
	keras_api

signatures
 
h

kernel
bias
 regularization_losses
!trainable_variables
"	variables
#	keras_api
h

$kernel
%bias
&regularization_losses
'trainable_variables
(	variables
)	keras_api
h

*kernel
+bias
,regularization_losses
-trainable_variables
.	variables
/	keras_api
R
0regularization_losses
1trainable_variables
2	variables
3	keras_api
R
4regularization_losses
5trainable_variables
6	variables
7	keras_api
R
8regularization_losses
9trainable_variables
:	variables
;	keras_api
R
<regularization_losses
=trainable_variables
>	variables
?	keras_api
R
@regularization_losses
Atrainable_variables
B	variables
C	keras_api
R
Dregularization_losses
Etrainable_variables
F	variables
G	keras_api
R
Hregularization_losses
Itrainable_variables
J	variables
K	keras_api
R
Lregularization_losses
Mtrainable_variables
N	variables
O	keras_api
R
Pregularization_losses
Qtrainable_variables
R	variables
S	keras_api
R
Tregularization_losses
Utrainable_variables
V	variables
W	keras_api
h

Xkernel
Ybias
Zregularization_losses
[trainable_variables
\	variables
]	keras_api
R
^regularization_losses
_trainable_variables
`	variables
a	keras_api
 
R
bregularization_losses
ctrainable_variables
d	variables
e	keras_api
h

fkernel
gbias
hregularization_losses
itrainable_variables
j	variables
k	keras_api
R
lregularization_losses
mtrainable_variables
n	variables
o	keras_api
h

pkernel
qbias
rregularization_losses
strainable_variables
t	variables
u	keras_api
R
vregularization_losses
wtrainable_variables
x	variables
y	keras_api
h

zkernel
{bias
|regularization_losses
}trainable_variables
~	variables
	keras_api
?
	?iter
?beta_1
?beta_2

?decay
?learning_ratem?m?$m?%m?*m?+m?Xm?Ym?fm?gm?pm?qm?zm?{m?v?v?$v?%v?*v?+v?Xv?Yv?fv?gv?pv?qv?zv?{v?
 
f
0
1
$2
%3
*4
+5
X6
Y7
f8
g9
p10
q11
z12
{13
f
0
1
$2
%3
*4
+5
X6
Y7
f8
g9
p10
q11
z12
{13
?
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?layer_metrics
?layers
regularization_losses
trainable_variables
	variables
 
YW
VARIABLE_VALUEconv_3/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEconv_3/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE
 

0
1

0
1
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
 regularization_losses
!trainable_variables
"	variables
YW
VARIABLE_VALUEconv_5/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEconv_5/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE
 

$0
%1

$0
%1
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
&regularization_losses
'trainable_variables
(	variables
YW
VARIABLE_VALUEconv_7/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEconv_7/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE
 

*0
+1

*0
+1
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
,regularization_losses
-trainable_variables
.	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
0regularization_losses
1trainable_variables
2	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
4regularization_losses
5trainable_variables
6	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
8regularization_losses
9trainable_variables
:	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
<regularization_losses
=trainable_variables
>	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
@regularization_losses
Atrainable_variables
B	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Dregularization_losses
Etrainable_variables
F	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Hregularization_losses
Itrainable_variables
J	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Lregularization_losses
Mtrainable_variables
N	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Pregularization_losses
Qtrainable_variables
R	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Tregularization_losses
Utrainable_variables
V	variables
ZX
VARIABLE_VALUEdense_0/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_0/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE
 

X0
Y1

X0
Y1
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Zregularization_losses
[trainable_variables
\	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
^regularization_losses
_trainable_variables
`	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
bregularization_losses
ctrainable_variables
d	variables
ZX
VARIABLE_VALUEdense_1/kernel6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_1/bias4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUE
 

f0
g1

f0
g1
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
hregularization_losses
itrainable_variables
j	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
lregularization_losses
mtrainable_variables
n	variables
ZX
VARIABLE_VALUEdense_2/kernel6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_2/bias4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUE
 

p0
q1

p0
q1
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
rregularization_losses
strainable_variables
t	variables
 
 
 
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
vregularization_losses
wtrainable_variables
x	variables
YW
VARIABLE_VALUEoutput/kernel6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEoutput/bias4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUE
 

z0
{1

z0
{1
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
|regularization_losses
}trainable_variables
~	variables
HF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
ZX
VARIABLE_VALUEAdam/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE

?0
?1
?2
 
 
 
?
0
1
2
3
4
5
6
7
	8

9
10
11
12
13
14
15
16
17
18
19
20
21
22
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
8

?total

?count
?	variables
?	keras_api
I

?total

?count
?
_fn_kwargs
?	variables
?	keras_api
I

?total

?count
?
_fn_kwargs
?	variables
?	keras_api
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

?0
?1

?	variables
QO
VARIABLE_VALUEtotal_14keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_14keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE
 

?0
?1

?	variables
QO
VARIABLE_VALUEtotal_24keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_24keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUE
 

?0
?1

?	variables
|z
VARIABLE_VALUEAdam/conv_3/kernel/mRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv_3/bias/mPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/conv_5/kernel/mRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv_5/bias/mPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/conv_7/kernel/mRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv_7/bias/mPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_0/kernel/mRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_0/bias/mPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_1/kernel/mRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_1/bias/mPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_2/kernel/mRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_2/bias/mPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/output/kernel/mRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/output/bias/mPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/conv_3/kernel/vRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv_3/bias/vPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/conv_5/kernel/vRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv_5/bias/vPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/conv_7/kernel/vRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv_7/bias/vPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_0/kernel/vRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_0/bias/vPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_1/kernel/vRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_1/bias/vPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_2/kernel/vRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_2/bias/vPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/output/kernel/vRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/output/bias/vPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
|
serving_default_input_dGBPlaceholder*'
_output_shapes
:?????????*
dtype0*
shape:?????????
?
serving_default_input_onehotPlaceholder*+
_output_shapes
:?????????*
dtype0* 
shape:?????????
?
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_dGBserving_default_input_onehotconv_7/kernelconv_7/biasconv_5/kernelconv_5/biasconv_3/kernelconv_3/biasdense_0/kerneldense_0/biasdense_1/kerneldense_1/biasdense_2/kerneldense_2/biasoutput/kerneloutput/bias*
Tin
2*
Tout
2*'
_output_shapes
:?????????*0
_read_only_resource_inputs
	
**
config_proto

CPU

GPU 2J 8*-
f(R&
$__inference_signature_wrapper_377744
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
?
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename!conv_3/kernel/Read/ReadVariableOpconv_3/bias/Read/ReadVariableOp!conv_5/kernel/Read/ReadVariableOpconv_5/bias/Read/ReadVariableOp!conv_7/kernel/Read/ReadVariableOpconv_7/bias/Read/ReadVariableOp"dense_0/kernel/Read/ReadVariableOp dense_0/bias/Read/ReadVariableOp"dense_1/kernel/Read/ReadVariableOp dense_1/bias/Read/ReadVariableOp"dense_2/kernel/Read/ReadVariableOp dense_2/bias/Read/ReadVariableOp!output/kernel/Read/ReadVariableOpoutput/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal_2/Read/ReadVariableOpcount_2/Read/ReadVariableOp(Adam/conv_3/kernel/m/Read/ReadVariableOp&Adam/conv_3/bias/m/Read/ReadVariableOp(Adam/conv_5/kernel/m/Read/ReadVariableOp&Adam/conv_5/bias/m/Read/ReadVariableOp(Adam/conv_7/kernel/m/Read/ReadVariableOp&Adam/conv_7/bias/m/Read/ReadVariableOp)Adam/dense_0/kernel/m/Read/ReadVariableOp'Adam/dense_0/bias/m/Read/ReadVariableOp)Adam/dense_1/kernel/m/Read/ReadVariableOp'Adam/dense_1/bias/m/Read/ReadVariableOp)Adam/dense_2/kernel/m/Read/ReadVariableOp'Adam/dense_2/bias/m/Read/ReadVariableOp(Adam/output/kernel/m/Read/ReadVariableOp&Adam/output/bias/m/Read/ReadVariableOp(Adam/conv_3/kernel/v/Read/ReadVariableOp&Adam/conv_3/bias/v/Read/ReadVariableOp(Adam/conv_5/kernel/v/Read/ReadVariableOp&Adam/conv_5/bias/v/Read/ReadVariableOp(Adam/conv_7/kernel/v/Read/ReadVariableOp&Adam/conv_7/bias/v/Read/ReadVariableOp)Adam/dense_0/kernel/v/Read/ReadVariableOp'Adam/dense_0/bias/v/Read/ReadVariableOp)Adam/dense_1/kernel/v/Read/ReadVariableOp'Adam/dense_1/bias/v/Read/ReadVariableOp)Adam/dense_2/kernel/v/Read/ReadVariableOp'Adam/dense_2/bias/v/Read/ReadVariableOp(Adam/output/kernel/v/Read/ReadVariableOp&Adam/output/bias/v/Read/ReadVariableOpConst*B
Tin;
927	*
Tout
2*
_output_shapes
: * 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*(
f#R!
__inference__traced_save_378535
?	
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameconv_3/kernelconv_3/biasconv_5/kernelconv_5/biasconv_7/kernelconv_7/biasdense_0/kerneldense_0/biasdense_1/kerneldense_1/biasdense_2/kerneldense_2/biasoutput/kerneloutput/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_ratetotalcounttotal_1count_1total_2count_2Adam/conv_3/kernel/mAdam/conv_3/bias/mAdam/conv_5/kernel/mAdam/conv_5/bias/mAdam/conv_7/kernel/mAdam/conv_7/bias/mAdam/dense_0/kernel/mAdam/dense_0/bias/mAdam/dense_1/kernel/mAdam/dense_1/bias/mAdam/dense_2/kernel/mAdam/dense_2/bias/mAdam/output/kernel/mAdam/output/bias/mAdam/conv_3/kernel/vAdam/conv_3/bias/vAdam/conv_5/kernel/vAdam/conv_5/bias/vAdam/conv_7/kernel/vAdam/conv_7/bias/vAdam/dense_0/kernel/vAdam/dense_0/bias/vAdam/dense_1/kernel/vAdam/dense_1/bias/vAdam/dense_2/kernel/vAdam/dense_2/bias/vAdam/output/kernel/vAdam/output/bias/v*A
Tin:
826*
Tout
2*
_output_shapes
: * 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*+
f&R$
"__inference__traced_restore_378706??
?
Z
.__inference_concatenate_1_layer_call_fn_378235
inputs_0
inputs_1
identity?
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*'
_output_shapes
:?????????Q* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*R
fMRK
I__inference_concatenate_1_layer_call_and_return_conditional_losses_3773192
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:?????????Q2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????P:?????????:Q M
'
_output_shapes
:?????????P
"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????
"
_user_specified_name
inputs/1
?
a
E__inference_flatten_3_layer_call_and_return_conditional_losses_377201

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"????x  2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????
2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????d:S O
+
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
a
C__inference_drop_d1_layer_call_and_return_conditional_losses_378272

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:?????????P2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:?????????P2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
?
&__inference_model_layer_call_fn_377612
input_onehot
	input_dgb
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinput_onehot	input_dgbunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:?????????*0
_read_only_resource_inputs
	
**
config_proto

CPU

GPU 2J 8*J
fERC
A__inference_model_layer_call_and_return_conditional_losses_3775812
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:Y U
+
_output_shapes
:?????????
&
_user_specified_nameinput_onehot:RN
'
_output_shapes
:?????????
#
_user_specified_name	input_dGB:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?L
?
A__inference_model_layer_call_and_return_conditional_losses_377523
input_onehot
	input_dgb
conv_7_377473
conv_7_377475
conv_5_377478
conv_5_377480
conv_3_377483
conv_3_377485
dense_0_377498
dense_0_377500
dense_1_377505
dense_1_377507
dense_2_377511
dense_2_377513
output_377517
output_377519
identity??conv_3/StatefulPartitionedCall?conv_5/StatefulPartitionedCall?conv_7/StatefulPartitionedCall?dense_0/StatefulPartitionedCall?dense_1/StatefulPartitionedCall?dense_2/StatefulPartitionedCall?output/StatefulPartitionedCall?
conv_7/StatefulPartitionedCallStatefulPartitionedCallinput_onehotconv_7_377473conv_7_377475*
Tin
2*
Tout
2*+
_output_shapes
:?????????(*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_7_layer_call_and_return_conditional_losses_3770272 
conv_7/StatefulPartitionedCall?
conv_5/StatefulPartitionedCallStatefulPartitionedCallinput_onehotconv_5_377478conv_5_377480*
Tin
2*
Tout
2*+
_output_shapes
:?????????F*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_5_layer_call_and_return_conditional_losses_3770002 
conv_5/StatefulPartitionedCall?
conv_3/StatefulPartitionedCallStatefulPartitionedCallinput_onehotconv_3_377483conv_3_377485*
Tin
2*
Tout
2*+
_output_shapes
:?????????d*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_3_layer_call_and_return_conditional_losses_3769732 
conv_3/StatefulPartitionedCall?
drop_7/PartitionedCallPartitionedCall'conv_7/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_7_layer_call_and_return_conditional_losses_3771192
drop_7/PartitionedCall?
drop_5/PartitionedCallPartitionedCall'conv_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_5_layer_call_and_return_conditional_losses_3771492
drop_5/PartitionedCall?
drop_3/PartitionedCallPartitionedCall'conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_3_layer_call_and_return_conditional_losses_3771792
drop_3/PartitionedCall?
pool_7/PartitionedCallPartitionedCalldrop_7/PartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_7_layer_call_and_return_conditional_losses_3770762
pool_7/PartitionedCall?
pool_5/PartitionedCallPartitionedCalldrop_5/PartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_5_layer_call_and_return_conditional_losses_3770612
pool_5/PartitionedCall?
pool_3/PartitionedCallPartitionedCalldrop_3/PartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_3_layer_call_and_return_conditional_losses_3770462
pool_3/PartitionedCall?
flatten_3/PartitionedCallPartitionedCallpool_3/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_3_layer_call_and_return_conditional_losses_3772012
flatten_3/PartitionedCall?
flatten_5/PartitionedCallPartitionedCallpool_5/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_5_layer_call_and_return_conditional_losses_3772152
flatten_5/PartitionedCall?
flatten_7/PartitionedCallPartitionedCallpool_7/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_7_layer_call_and_return_conditional_losses_3772292
flatten_7/PartitionedCall?
concatenate/PartitionedCallPartitionedCall"flatten_3/PartitionedCall:output:0"flatten_5/PartitionedCall:output:0"flatten_7/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*P
fKRI
G__inference_concatenate_layer_call_and_return_conditional_losses_3772452
concatenate/PartitionedCall?
dense_0/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_0_377498dense_0_377500*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_0_layer_call_and_return_conditional_losses_3772662!
dense_0/StatefulPartitionedCall?
drop_d0/PartitionedCallPartitionedCall(dense_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d0_layer_call_and_return_conditional_losses_3772992
drop_d0/PartitionedCall?
concatenate_1/PartitionedCallPartitionedCall drop_d0/PartitionedCall:output:0	input_dgb*
Tin
2*
Tout
2*'
_output_shapes
:?????????Q* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*R
fMRK
I__inference_concatenate_1_layer_call_and_return_conditional_losses_3773192
concatenate_1/PartitionedCall?
dense_1/StatefulPartitionedCallStatefulPartitionedCall&concatenate_1/PartitionedCall:output:0dense_1_377505dense_1_377507*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_1_layer_call_and_return_conditional_losses_3773392!
dense_1/StatefulPartitionedCall?
drop_d1/PartitionedCallPartitionedCall(dense_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d1_layer_call_and_return_conditional_losses_3773722
drop_d1/PartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCall drop_d1/PartitionedCall:output:0dense_2_377511dense_2_377513*
Tin
2*
Tout
2*'
_output_shapes
:?????????<*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_2_layer_call_and_return_conditional_losses_3773962!
dense_2/StatefulPartitionedCall?
drop_d2/PartitionedCallPartitionedCall(dense_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*'
_output_shapes
:?????????<* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d2_layer_call_and_return_conditional_losses_3774292
drop_d2/PartitionedCall?
output/StatefulPartitionedCallStatefulPartitionedCall drop_d2/PartitionedCall:output:0output_377517output_377519*
Tin
2*
Tout
2*'
_output_shapes
:?????????*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_output_layer_call_and_return_conditional_losses_3774522 
output/StatefulPartitionedCall?
IdentityIdentity'output/StatefulPartitionedCall:output:0^conv_3/StatefulPartitionedCall^conv_5/StatefulPartitionedCall^conv_7/StatefulPartitionedCall ^dense_0/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall ^dense_2/StatefulPartitionedCall^output/StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????::::::::::::::2@
conv_3/StatefulPartitionedCallconv_3/StatefulPartitionedCall2@
conv_5/StatefulPartitionedCallconv_5/StatefulPartitionedCall2@
conv_7/StatefulPartitionedCallconv_7/StatefulPartitionedCall2B
dense_0/StatefulPartitionedCalldense_0/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2@
output/StatefulPartitionedCalloutput/StatefulPartitionedCall:Y U
+
_output_shapes
:?????????
&
_user_specified_nameinput_onehot:RN
'
_output_shapes
:?????????
#
_user_specified_name	input_dGB:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
b
C__inference_drop_d1_layer_call_and_return_conditional_losses_378267

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:?????????P2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:?????????P*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:?????????P2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:?????????P2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:?????????P2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
a
B__inference_drop_5_layer_call_and_return_conditional_losses_378085

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:?????????F2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:?????????F*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:?????????F2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:?????????F2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:?????????F2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:?????????F2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????F:S O
+
_output_shapes
:?????????F
 
_user_specified_nameinputs
?
?
C__inference_dense_1_layer_call_and_return_conditional_losses_377339

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:QP*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:P*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????Q:::O K
'
_output_shapes
:?????????Q
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
a
B__inference_drop_5_layer_call_and_return_conditional_losses_377144

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:?????????F2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:?????????F*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:?????????F2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:?????????F2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:?????????F2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:?????????F2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????F:S O
+
_output_shapes
:?????????F
 
_user_specified_nameinputs
?
`
B__inference_drop_5_layer_call_and_return_conditional_losses_378090

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:?????????F2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:?????????F2

Identity_1"!

identity_1Identity_1:output:0**
_input_shapes
:?????????F:S O
+
_output_shapes
:?????????F
 
_user_specified_nameinputs
?
a
C__inference_drop_d0_layer_call_and_return_conditional_losses_378212

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:?????????P2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:?????????P2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
b
C__inference_drop_d0_layer_call_and_return_conditional_losses_377294

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:?????????P2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:?????????P*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:?????????P2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:?????????P2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:?????????P2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
?
$__inference_signature_wrapper_377744
	input_dgb
input_onehot
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinput_onehot	input_dgbunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:?????????*0
_read_only_resource_inputs
	
**
config_proto

CPU

GPU 2J 8**
f%R#
!__inference__wrapped_model_3769562
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:R N
'
_output_shapes
:?????????
#
_user_specified_name	input_dGB:YU
+
_output_shapes
:?????????
&
_user_specified_nameinput_onehot:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
`
'__inference_drop_7_layer_call_fn_378122

inputs
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_7_layer_call_and_return_conditional_losses_3771142
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:?????????(2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????(22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:?????????(
 
_user_specified_nameinputs
?
F
*__inference_flatten_7_layer_call_fn_378160

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_7_layer_call_and_return_conditional_losses_3772292
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????(:S O
+
_output_shapes
:?????????(
 
_user_specified_nameinputs
?
b
C__inference_drop_d0_layer_call_and_return_conditional_losses_378207

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:?????????P2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:?????????P*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:?????????P2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:?????????P2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:?????????P2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
??
?
"__inference__traced_restore_378706
file_prefix"
assignvariableop_conv_3_kernel"
assignvariableop_1_conv_3_bias$
 assignvariableop_2_conv_5_kernel"
assignvariableop_3_conv_5_bias$
 assignvariableop_4_conv_7_kernel"
assignvariableop_5_conv_7_bias%
!assignvariableop_6_dense_0_kernel#
assignvariableop_7_dense_0_bias%
!assignvariableop_8_dense_1_kernel#
assignvariableop_9_dense_1_bias&
"assignvariableop_10_dense_2_kernel$
 assignvariableop_11_dense_2_bias%
!assignvariableop_12_output_kernel#
assignvariableop_13_output_bias!
assignvariableop_14_adam_iter#
assignvariableop_15_adam_beta_1#
assignvariableop_16_adam_beta_2"
assignvariableop_17_adam_decay*
&assignvariableop_18_adam_learning_rate
assignvariableop_19_total
assignvariableop_20_count
assignvariableop_21_total_1
assignvariableop_22_count_1
assignvariableop_23_total_2
assignvariableop_24_count_2,
(assignvariableop_25_adam_conv_3_kernel_m*
&assignvariableop_26_adam_conv_3_bias_m,
(assignvariableop_27_adam_conv_5_kernel_m*
&assignvariableop_28_adam_conv_5_bias_m,
(assignvariableop_29_adam_conv_7_kernel_m*
&assignvariableop_30_adam_conv_7_bias_m-
)assignvariableop_31_adam_dense_0_kernel_m+
'assignvariableop_32_adam_dense_0_bias_m-
)assignvariableop_33_adam_dense_1_kernel_m+
'assignvariableop_34_adam_dense_1_bias_m-
)assignvariableop_35_adam_dense_2_kernel_m+
'assignvariableop_36_adam_dense_2_bias_m,
(assignvariableop_37_adam_output_kernel_m*
&assignvariableop_38_adam_output_bias_m,
(assignvariableop_39_adam_conv_3_kernel_v*
&assignvariableop_40_adam_conv_3_bias_v,
(assignvariableop_41_adam_conv_5_kernel_v*
&assignvariableop_42_adam_conv_5_bias_v,
(assignvariableop_43_adam_conv_7_kernel_v*
&assignvariableop_44_adam_conv_7_bias_v-
)assignvariableop_45_adam_dense_0_kernel_v+
'assignvariableop_46_adam_dense_0_bias_v-
)assignvariableop_47_adam_dense_1_kernel_v+
'assignvariableop_48_adam_dense_1_bias_v-
)assignvariableop_49_adam_dense_2_kernel_v+
'assignvariableop_50_adam_dense_2_bias_v,
(assignvariableop_51_adam_output_kernel_v*
&assignvariableop_52_adam_output_bias_v
identity_54??AssignVariableOp?AssignVariableOp_1?AssignVariableOp_10?AssignVariableOp_11?AssignVariableOp_12?AssignVariableOp_13?AssignVariableOp_14?AssignVariableOp_15?AssignVariableOp_16?AssignVariableOp_17?AssignVariableOp_18?AssignVariableOp_19?AssignVariableOp_2?AssignVariableOp_20?AssignVariableOp_21?AssignVariableOp_22?AssignVariableOp_23?AssignVariableOp_24?AssignVariableOp_25?AssignVariableOp_26?AssignVariableOp_27?AssignVariableOp_28?AssignVariableOp_29?AssignVariableOp_3?AssignVariableOp_30?AssignVariableOp_31?AssignVariableOp_32?AssignVariableOp_33?AssignVariableOp_34?AssignVariableOp_35?AssignVariableOp_36?AssignVariableOp_37?AssignVariableOp_38?AssignVariableOp_39?AssignVariableOp_4?AssignVariableOp_40?AssignVariableOp_41?AssignVariableOp_42?AssignVariableOp_43?AssignVariableOp_44?AssignVariableOp_45?AssignVariableOp_46?AssignVariableOp_47?AssignVariableOp_48?AssignVariableOp_49?AssignVariableOp_5?AssignVariableOp_50?AssignVariableOp_51?AssignVariableOp_52?AssignVariableOp_6?AssignVariableOp_7?AssignVariableOp_8?AssignVariableOp_9?	RestoreV2?RestoreV2_1?
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:5*
dtype0*?
value?B?5B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE2
RestoreV2/tensor_names?
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:5*
dtype0*}
valuetBr5B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slices?
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*?
_output_shapes?
?:::::::::::::::::::::::::::::::::::::::::::::::::::::*C
dtypes9
725	2
	RestoreV2X
IdentityIdentityRestoreV2:tensors:0*
T0*
_output_shapes
:2

Identity?
AssignVariableOpAssignVariableOpassignvariableop_conv_3_kernelIdentity:output:0*
_output_shapes
 *
dtype02
AssignVariableOp\

Identity_1IdentityRestoreV2:tensors:1*
T0*
_output_shapes
:2

Identity_1?
AssignVariableOp_1AssignVariableOpassignvariableop_1_conv_3_biasIdentity_1:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_1\

Identity_2IdentityRestoreV2:tensors:2*
T0*
_output_shapes
:2

Identity_2?
AssignVariableOp_2AssignVariableOp assignvariableop_2_conv_5_kernelIdentity_2:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_2\

Identity_3IdentityRestoreV2:tensors:3*
T0*
_output_shapes
:2

Identity_3?
AssignVariableOp_3AssignVariableOpassignvariableop_3_conv_5_biasIdentity_3:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_3\

Identity_4IdentityRestoreV2:tensors:4*
T0*
_output_shapes
:2

Identity_4?
AssignVariableOp_4AssignVariableOp assignvariableop_4_conv_7_kernelIdentity_4:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_4\

Identity_5IdentityRestoreV2:tensors:5*
T0*
_output_shapes
:2

Identity_5?
AssignVariableOp_5AssignVariableOpassignvariableop_5_conv_7_biasIdentity_5:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_5\

Identity_6IdentityRestoreV2:tensors:6*
T0*
_output_shapes
:2

Identity_6?
AssignVariableOp_6AssignVariableOp!assignvariableop_6_dense_0_kernelIdentity_6:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_6\

Identity_7IdentityRestoreV2:tensors:7*
T0*
_output_shapes
:2

Identity_7?
AssignVariableOp_7AssignVariableOpassignvariableop_7_dense_0_biasIdentity_7:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_7\

Identity_8IdentityRestoreV2:tensors:8*
T0*
_output_shapes
:2

Identity_8?
AssignVariableOp_8AssignVariableOp!assignvariableop_8_dense_1_kernelIdentity_8:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_8\

Identity_9IdentityRestoreV2:tensors:9*
T0*
_output_shapes
:2

Identity_9?
AssignVariableOp_9AssignVariableOpassignvariableop_9_dense_1_biasIdentity_9:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_9_
Identity_10IdentityRestoreV2:tensors:10*
T0*
_output_shapes
:2
Identity_10?
AssignVariableOp_10AssignVariableOp"assignvariableop_10_dense_2_kernelIdentity_10:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_10_
Identity_11IdentityRestoreV2:tensors:11*
T0*
_output_shapes
:2
Identity_11?
AssignVariableOp_11AssignVariableOp assignvariableop_11_dense_2_biasIdentity_11:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_11_
Identity_12IdentityRestoreV2:tensors:12*
T0*
_output_shapes
:2
Identity_12?
AssignVariableOp_12AssignVariableOp!assignvariableop_12_output_kernelIdentity_12:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_12_
Identity_13IdentityRestoreV2:tensors:13*
T0*
_output_shapes
:2
Identity_13?
AssignVariableOp_13AssignVariableOpassignvariableop_13_output_biasIdentity_13:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_13_
Identity_14IdentityRestoreV2:tensors:14*
T0	*
_output_shapes
:2
Identity_14?
AssignVariableOp_14AssignVariableOpassignvariableop_14_adam_iterIdentity_14:output:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_14_
Identity_15IdentityRestoreV2:tensors:15*
T0*
_output_shapes
:2
Identity_15?
AssignVariableOp_15AssignVariableOpassignvariableop_15_adam_beta_1Identity_15:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_15_
Identity_16IdentityRestoreV2:tensors:16*
T0*
_output_shapes
:2
Identity_16?
AssignVariableOp_16AssignVariableOpassignvariableop_16_adam_beta_2Identity_16:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_16_
Identity_17IdentityRestoreV2:tensors:17*
T0*
_output_shapes
:2
Identity_17?
AssignVariableOp_17AssignVariableOpassignvariableop_17_adam_decayIdentity_17:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_17_
Identity_18IdentityRestoreV2:tensors:18*
T0*
_output_shapes
:2
Identity_18?
AssignVariableOp_18AssignVariableOp&assignvariableop_18_adam_learning_rateIdentity_18:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_18_
Identity_19IdentityRestoreV2:tensors:19*
T0*
_output_shapes
:2
Identity_19?
AssignVariableOp_19AssignVariableOpassignvariableop_19_totalIdentity_19:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_19_
Identity_20IdentityRestoreV2:tensors:20*
T0*
_output_shapes
:2
Identity_20?
AssignVariableOp_20AssignVariableOpassignvariableop_20_countIdentity_20:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_20_
Identity_21IdentityRestoreV2:tensors:21*
T0*
_output_shapes
:2
Identity_21?
AssignVariableOp_21AssignVariableOpassignvariableop_21_total_1Identity_21:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_21_
Identity_22IdentityRestoreV2:tensors:22*
T0*
_output_shapes
:2
Identity_22?
AssignVariableOp_22AssignVariableOpassignvariableop_22_count_1Identity_22:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_22_
Identity_23IdentityRestoreV2:tensors:23*
T0*
_output_shapes
:2
Identity_23?
AssignVariableOp_23AssignVariableOpassignvariableop_23_total_2Identity_23:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_23_
Identity_24IdentityRestoreV2:tensors:24*
T0*
_output_shapes
:2
Identity_24?
AssignVariableOp_24AssignVariableOpassignvariableop_24_count_2Identity_24:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_24_
Identity_25IdentityRestoreV2:tensors:25*
T0*
_output_shapes
:2
Identity_25?
AssignVariableOp_25AssignVariableOp(assignvariableop_25_adam_conv_3_kernel_mIdentity_25:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_25_
Identity_26IdentityRestoreV2:tensors:26*
T0*
_output_shapes
:2
Identity_26?
AssignVariableOp_26AssignVariableOp&assignvariableop_26_adam_conv_3_bias_mIdentity_26:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_26_
Identity_27IdentityRestoreV2:tensors:27*
T0*
_output_shapes
:2
Identity_27?
AssignVariableOp_27AssignVariableOp(assignvariableop_27_adam_conv_5_kernel_mIdentity_27:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_27_
Identity_28IdentityRestoreV2:tensors:28*
T0*
_output_shapes
:2
Identity_28?
AssignVariableOp_28AssignVariableOp&assignvariableop_28_adam_conv_5_bias_mIdentity_28:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_28_
Identity_29IdentityRestoreV2:tensors:29*
T0*
_output_shapes
:2
Identity_29?
AssignVariableOp_29AssignVariableOp(assignvariableop_29_adam_conv_7_kernel_mIdentity_29:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_29_
Identity_30IdentityRestoreV2:tensors:30*
T0*
_output_shapes
:2
Identity_30?
AssignVariableOp_30AssignVariableOp&assignvariableop_30_adam_conv_7_bias_mIdentity_30:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_30_
Identity_31IdentityRestoreV2:tensors:31*
T0*
_output_shapes
:2
Identity_31?
AssignVariableOp_31AssignVariableOp)assignvariableop_31_adam_dense_0_kernel_mIdentity_31:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_31_
Identity_32IdentityRestoreV2:tensors:32*
T0*
_output_shapes
:2
Identity_32?
AssignVariableOp_32AssignVariableOp'assignvariableop_32_adam_dense_0_bias_mIdentity_32:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_32_
Identity_33IdentityRestoreV2:tensors:33*
T0*
_output_shapes
:2
Identity_33?
AssignVariableOp_33AssignVariableOp)assignvariableop_33_adam_dense_1_kernel_mIdentity_33:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_33_
Identity_34IdentityRestoreV2:tensors:34*
T0*
_output_shapes
:2
Identity_34?
AssignVariableOp_34AssignVariableOp'assignvariableop_34_adam_dense_1_bias_mIdentity_34:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_34_
Identity_35IdentityRestoreV2:tensors:35*
T0*
_output_shapes
:2
Identity_35?
AssignVariableOp_35AssignVariableOp)assignvariableop_35_adam_dense_2_kernel_mIdentity_35:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_35_
Identity_36IdentityRestoreV2:tensors:36*
T0*
_output_shapes
:2
Identity_36?
AssignVariableOp_36AssignVariableOp'assignvariableop_36_adam_dense_2_bias_mIdentity_36:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_36_
Identity_37IdentityRestoreV2:tensors:37*
T0*
_output_shapes
:2
Identity_37?
AssignVariableOp_37AssignVariableOp(assignvariableop_37_adam_output_kernel_mIdentity_37:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_37_
Identity_38IdentityRestoreV2:tensors:38*
T0*
_output_shapes
:2
Identity_38?
AssignVariableOp_38AssignVariableOp&assignvariableop_38_adam_output_bias_mIdentity_38:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_38_
Identity_39IdentityRestoreV2:tensors:39*
T0*
_output_shapes
:2
Identity_39?
AssignVariableOp_39AssignVariableOp(assignvariableop_39_adam_conv_3_kernel_vIdentity_39:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_39_
Identity_40IdentityRestoreV2:tensors:40*
T0*
_output_shapes
:2
Identity_40?
AssignVariableOp_40AssignVariableOp&assignvariableop_40_adam_conv_3_bias_vIdentity_40:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_40_
Identity_41IdentityRestoreV2:tensors:41*
T0*
_output_shapes
:2
Identity_41?
AssignVariableOp_41AssignVariableOp(assignvariableop_41_adam_conv_5_kernel_vIdentity_41:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_41_
Identity_42IdentityRestoreV2:tensors:42*
T0*
_output_shapes
:2
Identity_42?
AssignVariableOp_42AssignVariableOp&assignvariableop_42_adam_conv_5_bias_vIdentity_42:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_42_
Identity_43IdentityRestoreV2:tensors:43*
T0*
_output_shapes
:2
Identity_43?
AssignVariableOp_43AssignVariableOp(assignvariableop_43_adam_conv_7_kernel_vIdentity_43:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_43_
Identity_44IdentityRestoreV2:tensors:44*
T0*
_output_shapes
:2
Identity_44?
AssignVariableOp_44AssignVariableOp&assignvariableop_44_adam_conv_7_bias_vIdentity_44:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_44_
Identity_45IdentityRestoreV2:tensors:45*
T0*
_output_shapes
:2
Identity_45?
AssignVariableOp_45AssignVariableOp)assignvariableop_45_adam_dense_0_kernel_vIdentity_45:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_45_
Identity_46IdentityRestoreV2:tensors:46*
T0*
_output_shapes
:2
Identity_46?
AssignVariableOp_46AssignVariableOp'assignvariableop_46_adam_dense_0_bias_vIdentity_46:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_46_
Identity_47IdentityRestoreV2:tensors:47*
T0*
_output_shapes
:2
Identity_47?
AssignVariableOp_47AssignVariableOp)assignvariableop_47_adam_dense_1_kernel_vIdentity_47:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_47_
Identity_48IdentityRestoreV2:tensors:48*
T0*
_output_shapes
:2
Identity_48?
AssignVariableOp_48AssignVariableOp'assignvariableop_48_adam_dense_1_bias_vIdentity_48:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_48_
Identity_49IdentityRestoreV2:tensors:49*
T0*
_output_shapes
:2
Identity_49?
AssignVariableOp_49AssignVariableOp)assignvariableop_49_adam_dense_2_kernel_vIdentity_49:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_49_
Identity_50IdentityRestoreV2:tensors:50*
T0*
_output_shapes
:2
Identity_50?
AssignVariableOp_50AssignVariableOp'assignvariableop_50_adam_dense_2_bias_vIdentity_50:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_50_
Identity_51IdentityRestoreV2:tensors:51*
T0*
_output_shapes
:2
Identity_51?
AssignVariableOp_51AssignVariableOp(assignvariableop_51_adam_output_kernel_vIdentity_51:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_51_
Identity_52IdentityRestoreV2:tensors:52*
T0*
_output_shapes
:2
Identity_52?
AssignVariableOp_52AssignVariableOp&assignvariableop_52_adam_output_bias_vIdentity_52:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_52?
RestoreV2_1/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2_1/tensor_names?
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*
valueB
B 2
RestoreV2_1/shape_and_slices?
RestoreV2_1	RestoreV2file_prefix!RestoreV2_1/tensor_names:output:0%RestoreV2_1/shape_and_slices:output:0
^RestoreV2"/device:CPU:0*
_output_shapes
:*
dtypes
22
RestoreV2_19
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp?	
Identity_53Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_53?	
Identity_54IdentityIdentity_53:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9
^RestoreV2^RestoreV2_1*
T0*
_output_shapes
: 2
Identity_54"#
identity_54Identity_54:output:0*?
_input_shapes?
?: :::::::::::::::::::::::::::::::::::::::::::::::::::::2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_5AssignVariableOp_52*
AssignVariableOp_50AssignVariableOp_502*
AssignVariableOp_51AssignVariableOp_512*
AssignVariableOp_52AssignVariableOp_522(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_92
	RestoreV2	RestoreV22
RestoreV2_1RestoreV2_1:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: : 

_output_shapes
: :!

_output_shapes
: :"

_output_shapes
: :#

_output_shapes
: :$

_output_shapes
: :%

_output_shapes
: :&

_output_shapes
: :'

_output_shapes
: :(

_output_shapes
: :)

_output_shapes
: :*

_output_shapes
: :+

_output_shapes
: :,

_output_shapes
: :-

_output_shapes
: :.

_output_shapes
: :/

_output_shapes
: :0

_output_shapes
: :1

_output_shapes
: :2

_output_shapes
: :3

_output_shapes
: :4

_output_shapes
: :5

_output_shapes
: 
?
`
B__inference_drop_3_layer_call_and_return_conditional_losses_378063

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:?????????d2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:?????????d2

Identity_1"!

identity_1Identity_1:output:0**
_input_shapes
:?????????d:S O
+
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
a
E__inference_flatten_3_layer_call_and_return_conditional_losses_378133

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"????x  2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????
2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????d:S O
+
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
D
(__inference_drop_d0_layer_call_fn_378222

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d0_layer_call_and_return_conditional_losses_3772992
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
?
&__inference_model_layer_call_fn_377700
input_onehot
	input_dgb
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinput_onehot	input_dgbunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:?????????*0
_read_only_resource_inputs
	
**
config_proto

CPU

GPU 2J 8*J
fERC
A__inference_model_layer_call_and_return_conditional_losses_3776692
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:Y U
+
_output_shapes
:?????????
&
_user_specified_nameinput_onehot:RN
'
_output_shapes
:?????????
#
_user_specified_name	input_dGB:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
u
I__inference_concatenate_1_layer_call_and_return_conditional_losses_378229
inputs_0
inputs_1
identity\
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concat/axis?
concatConcatV2inputs_0inputs_1concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????Q2
concatc
IdentityIdentityconcat:output:0*
T0*'
_output_shapes
:?????????Q2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????P:?????????:Q M
'
_output_shapes
:?????????P
"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????
"
_user_specified_name
inputs/1
?
a
(__inference_drop_d0_layer_call_fn_378217

inputs
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d0_layer_call_and_return_conditional_losses_3772942
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????P22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
D
(__inference_drop_d1_layer_call_fn_378282

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d1_layer_call_and_return_conditional_losses_3773722
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
?
C__inference_dense_2_layer_call_and_return_conditional_losses_378293

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:P<*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:<*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:?????????<2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:?????????<2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????P:::O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
C
'__inference_drop_3_layer_call_fn_378073

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_3_layer_call_and_return_conditional_losses_3771792
PartitionedCallp
IdentityIdentityPartitionedCall:output:0*
T0*+
_output_shapes
:?????????d2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????d:S O
+
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
^
B__inference_pool_3_layer_call_and_return_conditional_losses_377046

inputs
identityb
ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
ExpandDims/dim?

ExpandDims
ExpandDimsinputsExpandDims/dim:output:0*
T0*A
_output_shapes/
-:+???????????????????????????2

ExpandDims?
AvgPoolAvgPoolExpandDims:output:0*
T0*A
_output_shapes/
-:+???????????????????????????*
ksize
*
paddingSAME*
strides
2	
AvgPool?
SqueezeSqueezeAvgPool:output:0*
T0*=
_output_shapes+
):'???????????????????????????*
squeeze_dims
2	
Squeezez
IdentityIdentitySqueeze:output:0*
T0*=
_output_shapes+
):'???????????????????????????2

Identity"
identityIdentity:output:0*<
_input_shapes+
):'???????????????????????????:e a
=
_output_shapes+
):'???????????????????????????
 
_user_specified_nameinputs
?
b
C__inference_drop_d1_layer_call_and_return_conditional_losses_377367

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:?????????P2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:?????????P*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:?????????P2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:?????????P2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:?????????P2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
^
B__inference_pool_7_layer_call_and_return_conditional_losses_377076

inputs
identityb
ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
ExpandDims/dim?

ExpandDims
ExpandDimsinputsExpandDims/dim:output:0*
T0*A
_output_shapes/
-:+???????????????????????????2

ExpandDims?
AvgPoolAvgPoolExpandDims:output:0*
T0*A
_output_shapes/
-:+???????????????????????????*
ksize
*
paddingSAME*
strides
2	
AvgPool?
SqueezeSqueezeAvgPool:output:0*
T0*=
_output_shapes+
):'???????????????????????????*
squeeze_dims
2	
Squeezez
IdentityIdentitySqueeze:output:0*
T0*=
_output_shapes+
):'???????????????????????????2

Identity"
identityIdentity:output:0*<
_input_shapes+
):'???????????????????????????:e a
=
_output_shapes+
):'???????????????????????????
 
_user_specified_nameinputs
?
F
*__inference_flatten_3_layer_call_fn_378138

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_3_layer_call_and_return_conditional_losses_3772012
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????d:S O
+
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
C
'__inference_drop_5_layer_call_fn_378100

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_5_layer_call_and_return_conditional_losses_3771492
PartitionedCallp
IdentityIdentityPartitionedCall:output:0*
T0*+
_output_shapes
:?????????F2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????F:S O
+
_output_shapes
:?????????F
 
_user_specified_nameinputs
?
?
C__inference_dense_0_layer_call_and_return_conditional_losses_377266

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?P*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:P*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*/
_input_shapes
:??????????:::P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?

G__inference_concatenate_layer_call_and_return_conditional_losses_377245

inputs
inputs_1
inputs_2
identity\
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concat/axis?
concatConcatV2inputsinputs_1inputs_2concat/axis:output:0*
N*
T0*(
_output_shapes
:??????????2
concatd
IdentityIdentityconcat:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*O
_input_shapes>
<:??????????
:??????????:??????????:P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs:PL
(
_output_shapes
:??????????
 
_user_specified_nameinputs:PL
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
^
B__inference_pool_5_layer_call_and_return_conditional_losses_377061

inputs
identityb
ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
ExpandDims/dim?

ExpandDims
ExpandDimsinputsExpandDims/dim:output:0*
T0*A
_output_shapes/
-:+???????????????????????????2

ExpandDims?
AvgPoolAvgPoolExpandDims:output:0*
T0*A
_output_shapes/
-:+???????????????????????????*
ksize
*
paddingSAME*
strides
2	
AvgPool?
SqueezeSqueezeAvgPool:output:0*
T0*=
_output_shapes+
):'???????????????????????????*
squeeze_dims
2	
Squeezez
IdentityIdentitySqueeze:output:0*
T0*=
_output_shapes+
):'???????????????????????????2

Identity"
identityIdentity:output:0*<
_input_shapes+
):'???????????????????????????:e a
=
_output_shapes+
):'???????????????????????????
 
_user_specified_nameinputs
?
?
B__inference_conv_7_layer_call_and_return_conditional_losses_377027

inputs/
+conv1d_expanddims_1_readvariableop_resource#
biasadd_readvariableop_resource
identity?p
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
conv1d/ExpandDims/dim?
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*8
_output_shapes&
$:"??????????????????2
conv1d/ExpandDims?
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:(*
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dim?
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:(2
conv1d/ExpandDims_1?
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*8
_output_shapes&
$:"??????????????????(*
paddingVALID*
strides
2
conv1d?
conv1d/SqueezeSqueezeconv1d:output:0*
T0*4
_output_shapes"
 :??????????????????(*
squeeze_dims
2
conv1d/Squeeze?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*4
_output_shapes"
 :??????????????????(2	
BiasAdde
ReluReluBiasAdd:output:0*
T0*4
_output_shapes"
 :??????????????????(2
Relus
IdentityIdentityRelu:activations:0*
T0*4
_output_shapes"
 :??????????????????(2

Identity"
identityIdentity:output:0*;
_input_shapes*
(:??????????????????:::\ X
4
_output_shapes"
 :??????????????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?w
?
!__inference__wrapped_model_376956
input_onehot
	input_dgb<
8model_conv_7_conv1d_expanddims_1_readvariableop_resource0
,model_conv_7_biasadd_readvariableop_resource<
8model_conv_5_conv1d_expanddims_1_readvariableop_resource0
,model_conv_5_biasadd_readvariableop_resource<
8model_conv_3_conv1d_expanddims_1_readvariableop_resource0
,model_conv_3_biasadd_readvariableop_resource0
,model_dense_0_matmul_readvariableop_resource1
-model_dense_0_biasadd_readvariableop_resource0
,model_dense_1_matmul_readvariableop_resource1
-model_dense_1_biasadd_readvariableop_resource0
,model_dense_2_matmul_readvariableop_resource1
-model_dense_2_biasadd_readvariableop_resource/
+model_output_matmul_readvariableop_resource0
,model_output_biasadd_readvariableop_resource
identity??
"model/conv_7/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2$
"model/conv_7/conv1d/ExpandDims/dim?
model/conv_7/conv1d/ExpandDims
ExpandDimsinput_onehot+model/conv_7/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????2 
model/conv_7/conv1d/ExpandDims?
/model/conv_7/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp8model_conv_7_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:(*
dtype021
/model/conv_7/conv1d/ExpandDims_1/ReadVariableOp?
$model/conv_7/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2&
$model/conv_7/conv1d/ExpandDims_1/dim?
 model/conv_7/conv1d/ExpandDims_1
ExpandDims7model/conv_7/conv1d/ExpandDims_1/ReadVariableOp:value:0-model/conv_7/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:(2"
 model/conv_7/conv1d/ExpandDims_1?
model/conv_7/conv1dConv2D'model/conv_7/conv1d/ExpandDims:output:0)model/conv_7/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:?????????(*
paddingVALID*
strides
2
model/conv_7/conv1d?
model/conv_7/conv1d/SqueezeSqueezemodel/conv_7/conv1d:output:0*
T0*+
_output_shapes
:?????????(*
squeeze_dims
2
model/conv_7/conv1d/Squeeze?
#model/conv_7/BiasAdd/ReadVariableOpReadVariableOp,model_conv_7_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02%
#model/conv_7/BiasAdd/ReadVariableOp?
model/conv_7/BiasAddBiasAdd$model/conv_7/conv1d/Squeeze:output:0+model/conv_7/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:?????????(2
model/conv_7/BiasAdd?
model/conv_7/ReluRelumodel/conv_7/BiasAdd:output:0*
T0*+
_output_shapes
:?????????(2
model/conv_7/Relu?
"model/conv_5/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2$
"model/conv_5/conv1d/ExpandDims/dim?
model/conv_5/conv1d/ExpandDims
ExpandDimsinput_onehot+model/conv_5/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????2 
model/conv_5/conv1d/ExpandDims?
/model/conv_5/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp8model_conv_5_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:F*
dtype021
/model/conv_5/conv1d/ExpandDims_1/ReadVariableOp?
$model/conv_5/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2&
$model/conv_5/conv1d/ExpandDims_1/dim?
 model/conv_5/conv1d/ExpandDims_1
ExpandDims7model/conv_5/conv1d/ExpandDims_1/ReadVariableOp:value:0-model/conv_5/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:F2"
 model/conv_5/conv1d/ExpandDims_1?
model/conv_5/conv1dConv2D'model/conv_5/conv1d/ExpandDims:output:0)model/conv_5/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:?????????F*
paddingVALID*
strides
2
model/conv_5/conv1d?
model/conv_5/conv1d/SqueezeSqueezemodel/conv_5/conv1d:output:0*
T0*+
_output_shapes
:?????????F*
squeeze_dims
2
model/conv_5/conv1d/Squeeze?
#model/conv_5/BiasAdd/ReadVariableOpReadVariableOp,model_conv_5_biasadd_readvariableop_resource*
_output_shapes
:F*
dtype02%
#model/conv_5/BiasAdd/ReadVariableOp?
model/conv_5/BiasAddBiasAdd$model/conv_5/conv1d/Squeeze:output:0+model/conv_5/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:?????????F2
model/conv_5/BiasAdd?
model/conv_5/ReluRelumodel/conv_5/BiasAdd:output:0*
T0*+
_output_shapes
:?????????F2
model/conv_5/Relu?
"model/conv_3/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2$
"model/conv_3/conv1d/ExpandDims/dim?
model/conv_3/conv1d/ExpandDims
ExpandDimsinput_onehot+model/conv_3/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????2 
model/conv_3/conv1d/ExpandDims?
/model/conv_3/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp8model_conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:d*
dtype021
/model/conv_3/conv1d/ExpandDims_1/ReadVariableOp?
$model/conv_3/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2&
$model/conv_3/conv1d/ExpandDims_1/dim?
 model/conv_3/conv1d/ExpandDims_1
ExpandDims7model/conv_3/conv1d/ExpandDims_1/ReadVariableOp:value:0-model/conv_3/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:d2"
 model/conv_3/conv1d/ExpandDims_1?
model/conv_3/conv1dConv2D'model/conv_3/conv1d/ExpandDims:output:0)model/conv_3/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:?????????d*
paddingVALID*
strides
2
model/conv_3/conv1d?
model/conv_3/conv1d/SqueezeSqueezemodel/conv_3/conv1d:output:0*
T0*+
_output_shapes
:?????????d*
squeeze_dims
2
model/conv_3/conv1d/Squeeze?
#model/conv_3/BiasAdd/ReadVariableOpReadVariableOp,model_conv_3_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02%
#model/conv_3/BiasAdd/ReadVariableOp?
model/conv_3/BiasAddBiasAdd$model/conv_3/conv1d/Squeeze:output:0+model/conv_3/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:?????????d2
model/conv_3/BiasAdd?
model/conv_3/ReluRelumodel/conv_3/BiasAdd:output:0*
T0*+
_output_shapes
:?????????d2
model/conv_3/Relu?
model/drop_7/IdentityIdentitymodel/conv_7/Relu:activations:0*
T0*+
_output_shapes
:?????????(2
model/drop_7/Identity?
model/drop_5/IdentityIdentitymodel/conv_5/Relu:activations:0*
T0*+
_output_shapes
:?????????F2
model/drop_5/Identity?
model/drop_3/IdentityIdentitymodel/conv_3/Relu:activations:0*
T0*+
_output_shapes
:?????????d2
model/drop_3/Identity|
model/pool_7/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
model/pool_7/ExpandDims/dim?
model/pool_7/ExpandDims
ExpandDimsmodel/drop_7/Identity:output:0$model/pool_7/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????(2
model/pool_7/ExpandDims?
model/pool_7/AvgPoolAvgPool model/pool_7/ExpandDims:output:0*
T0*/
_output_shapes
:?????????(*
ksize
*
paddingSAME*
strides
2
model/pool_7/AvgPool?
model/pool_7/SqueezeSqueezemodel/pool_7/AvgPool:output:0*
T0*+
_output_shapes
:?????????(*
squeeze_dims
2
model/pool_7/Squeeze|
model/pool_5/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
model/pool_5/ExpandDims/dim?
model/pool_5/ExpandDims
ExpandDimsmodel/drop_5/Identity:output:0$model/pool_5/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????F2
model/pool_5/ExpandDims?
model/pool_5/AvgPoolAvgPool model/pool_5/ExpandDims:output:0*
T0*/
_output_shapes
:?????????F*
ksize
*
paddingSAME*
strides
2
model/pool_5/AvgPool?
model/pool_5/SqueezeSqueezemodel/pool_5/AvgPool:output:0*
T0*+
_output_shapes
:?????????F*
squeeze_dims
2
model/pool_5/Squeeze|
model/pool_3/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
model/pool_3/ExpandDims/dim?
model/pool_3/ExpandDims
ExpandDimsmodel/drop_3/Identity:output:0$model/pool_3/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????d2
model/pool_3/ExpandDims?
model/pool_3/AvgPoolAvgPool model/pool_3/ExpandDims:output:0*
T0*/
_output_shapes
:?????????d*
ksize
*
paddingSAME*
strides
2
model/pool_3/AvgPool?
model/pool_3/SqueezeSqueezemodel/pool_3/AvgPool:output:0*
T0*+
_output_shapes
:?????????d*
squeeze_dims
2
model/pool_3/Squeeze
model/flatten_3/ConstConst*
_output_shapes
:*
dtype0*
valueB"????x  2
model/flatten_3/Const?
model/flatten_3/ReshapeReshapemodel/pool_3/Squeeze:output:0model/flatten_3/Const:output:0*
T0*(
_output_shapes
:??????????
2
model/flatten_3/Reshape
model/flatten_5/ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
model/flatten_5/Const?
model/flatten_5/ReshapeReshapemodel/pool_5/Squeeze:output:0model/flatten_5/Const:output:0*
T0*(
_output_shapes
:??????????2
model/flatten_5/Reshape
model/flatten_7/ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
model/flatten_7/Const?
model/flatten_7/ReshapeReshapemodel/pool_7/Squeeze:output:0model/flatten_7/Const:output:0*
T0*(
_output_shapes
:??????????2
model/flatten_7/Reshape?
model/concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
model/concatenate/concat/axis?
model/concatenate/concatConcatV2 model/flatten_3/Reshape:output:0 model/flatten_5/Reshape:output:0 model/flatten_7/Reshape:output:0&model/concatenate/concat/axis:output:0*
N*
T0*(
_output_shapes
:??????????2
model/concatenate/concat?
#model/dense_0/MatMul/ReadVariableOpReadVariableOp,model_dense_0_matmul_readvariableop_resource*
_output_shapes
:	?P*
dtype02%
#model/dense_0/MatMul/ReadVariableOp?
model/dense_0/MatMulMatMul!model/concatenate/concat:output:0+model/dense_0/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
model/dense_0/MatMul?
$model/dense_0/BiasAdd/ReadVariableOpReadVariableOp-model_dense_0_biasadd_readvariableop_resource*
_output_shapes
:P*
dtype02&
$model/dense_0/BiasAdd/ReadVariableOp?
model/dense_0/BiasAddBiasAddmodel/dense_0/MatMul:product:0,model/dense_0/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
model/dense_0/BiasAdd?
model/dense_0/ReluRelumodel/dense_0/BiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
model/dense_0/Relu?
model/drop_d0/IdentityIdentity model/dense_0/Relu:activations:0*
T0*'
_output_shapes
:?????????P2
model/drop_d0/Identity?
model/concatenate_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2!
model/concatenate_1/concat/axis?
model/concatenate_1/concatConcatV2model/drop_d0/Identity:output:0	input_dgb(model/concatenate_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????Q2
model/concatenate_1/concat?
#model/dense_1/MatMul/ReadVariableOpReadVariableOp,model_dense_1_matmul_readvariableop_resource*
_output_shapes

:QP*
dtype02%
#model/dense_1/MatMul/ReadVariableOp?
model/dense_1/MatMulMatMul#model/concatenate_1/concat:output:0+model/dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
model/dense_1/MatMul?
$model/dense_1/BiasAdd/ReadVariableOpReadVariableOp-model_dense_1_biasadd_readvariableop_resource*
_output_shapes
:P*
dtype02&
$model/dense_1/BiasAdd/ReadVariableOp?
model/dense_1/BiasAddBiasAddmodel/dense_1/MatMul:product:0,model/dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
model/dense_1/BiasAdd?
model/dense_1/ReluRelumodel/dense_1/BiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
model/dense_1/Relu?
model/drop_d1/IdentityIdentity model/dense_1/Relu:activations:0*
T0*'
_output_shapes
:?????????P2
model/drop_d1/Identity?
#model/dense_2/MatMul/ReadVariableOpReadVariableOp,model_dense_2_matmul_readvariableop_resource*
_output_shapes

:P<*
dtype02%
#model/dense_2/MatMul/ReadVariableOp?
model/dense_2/MatMulMatMulmodel/drop_d1/Identity:output:0+model/dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2
model/dense_2/MatMul?
$model/dense_2/BiasAdd/ReadVariableOpReadVariableOp-model_dense_2_biasadd_readvariableop_resource*
_output_shapes
:<*
dtype02&
$model/dense_2/BiasAdd/ReadVariableOp?
model/dense_2/BiasAddBiasAddmodel/dense_2/MatMul:product:0,model/dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2
model/dense_2/BiasAdd?
model/dense_2/ReluRelumodel/dense_2/BiasAdd:output:0*
T0*'
_output_shapes
:?????????<2
model/dense_2/Relu?
model/drop_d2/IdentityIdentity model/dense_2/Relu:activations:0*
T0*'
_output_shapes
:?????????<2
model/drop_d2/Identity?
"model/output/MatMul/ReadVariableOpReadVariableOp+model_output_matmul_readvariableop_resource*
_output_shapes

:<*
dtype02$
"model/output/MatMul/ReadVariableOp?
model/output/MatMulMatMulmodel/drop_d2/Identity:output:0*model/output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
model/output/MatMul?
#model/output/BiasAdd/ReadVariableOpReadVariableOp,model_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02%
#model/output/BiasAdd/ReadVariableOp?
model/output/BiasAddBiasAddmodel/output/MatMul:product:0+model/output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
model/output/BiasAddq
IdentityIdentitymodel/output/BiasAdd:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????:::::::::::::::Y U
+
_output_shapes
:?????????
&
_user_specified_nameinput_onehot:RN
'
_output_shapes
:?????????
#
_user_specified_name	input_dGB:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
|
'__inference_conv_3_layer_call_fn_376983

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*4
_output_shapes"
 :??????????????????d*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_3_layer_call_and_return_conditional_losses_3769732
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*4
_output_shapes"
 :??????????????????d2

Identity"
identityIdentity:output:0*;
_input_shapes*
(:??????????????????::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
?
&__inference_model_layer_call_fn_378012
inputs_0
inputs_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:?????????*0
_read_only_resource_inputs
	
**
config_proto

CPU

GPU 2J 8*J
fERC
A__inference_model_layer_call_and_return_conditional_losses_3775812
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:U Q
+
_output_shapes
:?????????
"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????
"
_user_specified_name
inputs/1:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
a
B__inference_drop_3_layer_call_and_return_conditional_losses_377174

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:?????????d2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:?????????d*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:?????????d2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:?????????d2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:?????????d2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:?????????d2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????d:S O
+
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
a
B__inference_drop_3_layer_call_and_return_conditional_losses_378058

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:?????????d2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:?????????d*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:?????????d2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:?????????d2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:?????????d2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:?????????d2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????d:S O
+
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
`
'__inference_drop_3_layer_call_fn_378068

inputs
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_3_layer_call_and_return_conditional_losses_3771742
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:?????????d2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????d22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:?????????d
 
_user_specified_nameinputs
?k
?
A__inference_model_layer_call_and_return_conditional_losses_377978
inputs_0
inputs_16
2conv_7_conv1d_expanddims_1_readvariableop_resource*
&conv_7_biasadd_readvariableop_resource6
2conv_5_conv1d_expanddims_1_readvariableop_resource*
&conv_5_biasadd_readvariableop_resource6
2conv_3_conv1d_expanddims_1_readvariableop_resource*
&conv_3_biasadd_readvariableop_resource*
&dense_0_matmul_readvariableop_resource+
'dense_0_biasadd_readvariableop_resource*
&dense_1_matmul_readvariableop_resource+
'dense_1_biasadd_readvariableop_resource*
&dense_2_matmul_readvariableop_resource+
'dense_2_biasadd_readvariableop_resource)
%output_matmul_readvariableop_resource*
&output_biasadd_readvariableop_resource
identity?~
conv_7/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
conv_7/conv1d/ExpandDims/dim?
conv_7/conv1d/ExpandDims
ExpandDimsinputs_0%conv_7/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????2
conv_7/conv1d/ExpandDims?
)conv_7/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_7_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:(*
dtype02+
)conv_7/conv1d/ExpandDims_1/ReadVariableOp?
conv_7/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
conv_7/conv1d/ExpandDims_1/dim?
conv_7/conv1d/ExpandDims_1
ExpandDims1conv_7/conv1d/ExpandDims_1/ReadVariableOp:value:0'conv_7/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:(2
conv_7/conv1d/ExpandDims_1?
conv_7/conv1dConv2D!conv_7/conv1d/ExpandDims:output:0#conv_7/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:?????????(*
paddingVALID*
strides
2
conv_7/conv1d?
conv_7/conv1d/SqueezeSqueezeconv_7/conv1d:output:0*
T0*+
_output_shapes
:?????????(*
squeeze_dims
2
conv_7/conv1d/Squeeze?
conv_7/BiasAdd/ReadVariableOpReadVariableOp&conv_7_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
conv_7/BiasAdd/ReadVariableOp?
conv_7/BiasAddBiasAddconv_7/conv1d/Squeeze:output:0%conv_7/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:?????????(2
conv_7/BiasAddq
conv_7/ReluReluconv_7/BiasAdd:output:0*
T0*+
_output_shapes
:?????????(2
conv_7/Relu~
conv_5/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
conv_5/conv1d/ExpandDims/dim?
conv_5/conv1d/ExpandDims
ExpandDimsinputs_0%conv_5/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????2
conv_5/conv1d/ExpandDims?
)conv_5/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_5_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:F*
dtype02+
)conv_5/conv1d/ExpandDims_1/ReadVariableOp?
conv_5/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
conv_5/conv1d/ExpandDims_1/dim?
conv_5/conv1d/ExpandDims_1
ExpandDims1conv_5/conv1d/ExpandDims_1/ReadVariableOp:value:0'conv_5/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:F2
conv_5/conv1d/ExpandDims_1?
conv_5/conv1dConv2D!conv_5/conv1d/ExpandDims:output:0#conv_5/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:?????????F*
paddingVALID*
strides
2
conv_5/conv1d?
conv_5/conv1d/SqueezeSqueezeconv_5/conv1d:output:0*
T0*+
_output_shapes
:?????????F*
squeeze_dims
2
conv_5/conv1d/Squeeze?
conv_5/BiasAdd/ReadVariableOpReadVariableOp&conv_5_biasadd_readvariableop_resource*
_output_shapes
:F*
dtype02
conv_5/BiasAdd/ReadVariableOp?
conv_5/BiasAddBiasAddconv_5/conv1d/Squeeze:output:0%conv_5/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:?????????F2
conv_5/BiasAddq
conv_5/ReluReluconv_5/BiasAdd:output:0*
T0*+
_output_shapes
:?????????F2
conv_5/Relu~
conv_3/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
conv_3/conv1d/ExpandDims/dim?
conv_3/conv1d/ExpandDims
ExpandDimsinputs_0%conv_3/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????2
conv_3/conv1d/ExpandDims?
)conv_3/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:d*
dtype02+
)conv_3/conv1d/ExpandDims_1/ReadVariableOp?
conv_3/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
conv_3/conv1d/ExpandDims_1/dim?
conv_3/conv1d/ExpandDims_1
ExpandDims1conv_3/conv1d/ExpandDims_1/ReadVariableOp:value:0'conv_3/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:d2
conv_3/conv1d/ExpandDims_1?
conv_3/conv1dConv2D!conv_3/conv1d/ExpandDims:output:0#conv_3/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:?????????d*
paddingVALID*
strides
2
conv_3/conv1d?
conv_3/conv1d/SqueezeSqueezeconv_3/conv1d:output:0*
T0*+
_output_shapes
:?????????d*
squeeze_dims
2
conv_3/conv1d/Squeeze?
conv_3/BiasAdd/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02
conv_3/BiasAdd/ReadVariableOp?
conv_3/BiasAddBiasAddconv_3/conv1d/Squeeze:output:0%conv_3/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:?????????d2
conv_3/BiasAddq
conv_3/ReluReluconv_3/BiasAdd:output:0*
T0*+
_output_shapes
:?????????d2
conv_3/Relu
drop_7/IdentityIdentityconv_7/Relu:activations:0*
T0*+
_output_shapes
:?????????(2
drop_7/Identity
drop_5/IdentityIdentityconv_5/Relu:activations:0*
T0*+
_output_shapes
:?????????F2
drop_5/Identity
drop_3/IdentityIdentityconv_3/Relu:activations:0*
T0*+
_output_shapes
:?????????d2
drop_3/Identityp
pool_7/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
pool_7/ExpandDims/dim?
pool_7/ExpandDims
ExpandDimsdrop_7/Identity:output:0pool_7/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????(2
pool_7/ExpandDims?
pool_7/AvgPoolAvgPoolpool_7/ExpandDims:output:0*
T0*/
_output_shapes
:?????????(*
ksize
*
paddingSAME*
strides
2
pool_7/AvgPool?
pool_7/SqueezeSqueezepool_7/AvgPool:output:0*
T0*+
_output_shapes
:?????????(*
squeeze_dims
2
pool_7/Squeezep
pool_5/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
pool_5/ExpandDims/dim?
pool_5/ExpandDims
ExpandDimsdrop_5/Identity:output:0pool_5/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????F2
pool_5/ExpandDims?
pool_5/AvgPoolAvgPoolpool_5/ExpandDims:output:0*
T0*/
_output_shapes
:?????????F*
ksize
*
paddingSAME*
strides
2
pool_5/AvgPool?
pool_5/SqueezeSqueezepool_5/AvgPool:output:0*
T0*+
_output_shapes
:?????????F*
squeeze_dims
2
pool_5/Squeezep
pool_3/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
pool_3/ExpandDims/dim?
pool_3/ExpandDims
ExpandDimsdrop_3/Identity:output:0pool_3/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????d2
pool_3/ExpandDims?
pool_3/AvgPoolAvgPoolpool_3/ExpandDims:output:0*
T0*/
_output_shapes
:?????????d*
ksize
*
paddingSAME*
strides
2
pool_3/AvgPool?
pool_3/SqueezeSqueezepool_3/AvgPool:output:0*
T0*+
_output_shapes
:?????????d*
squeeze_dims
2
pool_3/Squeezes
flatten_3/ConstConst*
_output_shapes
:*
dtype0*
valueB"????x  2
flatten_3/Const?
flatten_3/ReshapeReshapepool_3/Squeeze:output:0flatten_3/Const:output:0*
T0*(
_output_shapes
:??????????
2
flatten_3/Reshapes
flatten_5/ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
flatten_5/Const?
flatten_5/ReshapeReshapepool_5/Squeeze:output:0flatten_5/Const:output:0*
T0*(
_output_shapes
:??????????2
flatten_5/Reshapes
flatten_7/ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
flatten_7/Const?
flatten_7/ReshapeReshapepool_7/Squeeze:output:0flatten_7/Const:output:0*
T0*(
_output_shapes
:??????????2
flatten_7/Reshapet
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate/concat/axis?
concatenate/concatConcatV2flatten_3/Reshape:output:0flatten_5/Reshape:output:0flatten_7/Reshape:output:0 concatenate/concat/axis:output:0*
N*
T0*(
_output_shapes
:??????????2
concatenate/concat?
dense_0/MatMul/ReadVariableOpReadVariableOp&dense_0_matmul_readvariableop_resource*
_output_shapes
:	?P*
dtype02
dense_0/MatMul/ReadVariableOp?
dense_0/MatMulMatMulconcatenate/concat:output:0%dense_0/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
dense_0/MatMul?
dense_0/BiasAdd/ReadVariableOpReadVariableOp'dense_0_biasadd_readvariableop_resource*
_output_shapes
:P*
dtype02 
dense_0/BiasAdd/ReadVariableOp?
dense_0/BiasAddBiasAdddense_0/MatMul:product:0&dense_0/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
dense_0/BiasAddp
dense_0/ReluReludense_0/BiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
dense_0/Relu~
drop_d0/IdentityIdentitydense_0/Relu:activations:0*
T0*'
_output_shapes
:?????????P2
drop_d0/Identityx
concatenate_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate_1/concat/axis?
concatenate_1/concatConcatV2drop_d0/Identity:output:0inputs_1"concatenate_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????Q2
concatenate_1/concat?
dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes

:QP*
dtype02
dense_1/MatMul/ReadVariableOp?
dense_1/MatMulMatMulconcatenate_1/concat:output:0%dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
dense_1/MatMul?
dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:P*
dtype02 
dense_1/BiasAdd/ReadVariableOp?
dense_1/BiasAddBiasAdddense_1/MatMul:product:0&dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
dense_1/BiasAddp
dense_1/ReluReludense_1/BiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
dense_1/Relu~
drop_d1/IdentityIdentitydense_1/Relu:activations:0*
T0*'
_output_shapes
:?????????P2
drop_d1/Identity?
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:P<*
dtype02
dense_2/MatMul/ReadVariableOp?
dense_2/MatMulMatMuldrop_d1/Identity:output:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2
dense_2/MatMul?
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:<*
dtype02 
dense_2/BiasAdd/ReadVariableOp?
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2
dense_2/BiasAddp
dense_2/ReluReludense_2/BiasAdd:output:0*
T0*'
_output_shapes
:?????????<2
dense_2/Relu~
drop_d2/IdentityIdentitydense_2/Relu:activations:0*
T0*'
_output_shapes
:?????????<2
drop_d2/Identity?
output/MatMul/ReadVariableOpReadVariableOp%output_matmul_readvariableop_resource*
_output_shapes

:<*
dtype02
output/MatMul/ReadVariableOp?
output/MatMulMatMuldrop_d2/Identity:output:0$output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
output/MatMul?
output/BiasAdd/ReadVariableOpReadVariableOp&output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
output/BiasAdd/ReadVariableOp?
output/BiasAddBiasAddoutput/MatMul:product:0%output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
output/BiasAddk
IdentityIdentityoutput/BiasAdd:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????:::::::::::::::U Q
+
_output_shapes
:?????????
"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????
"
_user_specified_name
inputs/1:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
a
C__inference_drop_d2_layer_call_and_return_conditional_losses_378319

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:?????????<2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:?????????<2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:?????????<:O K
'
_output_shapes
:?????????<
 
_user_specified_nameinputs
?
`
B__inference_drop_3_layer_call_and_return_conditional_losses_377179

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:?????????d2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:?????????d2

Identity_1"!

identity_1Identity_1:output:0**
_input_shapes
:?????????d:S O
+
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
a
E__inference_flatten_7_layer_call_and_return_conditional_losses_377229

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????(:S O
+
_output_shapes
:?????????(
 
_user_specified_nameinputs
?
b
C__inference_drop_d2_layer_call_and_return_conditional_losses_377424

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:?????????<2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:?????????<*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:?????????<2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:?????????<2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:?????????<2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:?????????<2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????<:O K
'
_output_shapes
:?????????<
 
_user_specified_nameinputs
?
C
'__inference_pool_3_layer_call_fn_377052

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*=
_output_shapes+
):'???????????????????????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_3_layer_call_and_return_conditional_losses_3770462
PartitionedCall?
IdentityIdentityPartitionedCall:output:0*
T0*=
_output_shapes+
):'???????????????????????????2

Identity"
identityIdentity:output:0*<
_input_shapes+
):'???????????????????????????:e a
=
_output_shapes+
):'???????????????????????????
 
_user_specified_nameinputs
?K
?
A__inference_model_layer_call_and_return_conditional_losses_377669

inputs
inputs_1
conv_7_377619
conv_7_377621
conv_5_377624
conv_5_377626
conv_3_377629
conv_3_377631
dense_0_377644
dense_0_377646
dense_1_377651
dense_1_377653
dense_2_377657
dense_2_377659
output_377663
output_377665
identity??conv_3/StatefulPartitionedCall?conv_5/StatefulPartitionedCall?conv_7/StatefulPartitionedCall?dense_0/StatefulPartitionedCall?dense_1/StatefulPartitionedCall?dense_2/StatefulPartitionedCall?output/StatefulPartitionedCall?
conv_7/StatefulPartitionedCallStatefulPartitionedCallinputsconv_7_377619conv_7_377621*
Tin
2*
Tout
2*+
_output_shapes
:?????????(*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_7_layer_call_and_return_conditional_losses_3770272 
conv_7/StatefulPartitionedCall?
conv_5/StatefulPartitionedCallStatefulPartitionedCallinputsconv_5_377624conv_5_377626*
Tin
2*
Tout
2*+
_output_shapes
:?????????F*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_5_layer_call_and_return_conditional_losses_3770002 
conv_5/StatefulPartitionedCall?
conv_3/StatefulPartitionedCallStatefulPartitionedCallinputsconv_3_377629conv_3_377631*
Tin
2*
Tout
2*+
_output_shapes
:?????????d*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_3_layer_call_and_return_conditional_losses_3769732 
conv_3/StatefulPartitionedCall?
drop_7/PartitionedCallPartitionedCall'conv_7/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_7_layer_call_and_return_conditional_losses_3771192
drop_7/PartitionedCall?
drop_5/PartitionedCallPartitionedCall'conv_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_5_layer_call_and_return_conditional_losses_3771492
drop_5/PartitionedCall?
drop_3/PartitionedCallPartitionedCall'conv_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_3_layer_call_and_return_conditional_losses_3771792
drop_3/PartitionedCall?
pool_7/PartitionedCallPartitionedCalldrop_7/PartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_7_layer_call_and_return_conditional_losses_3770762
pool_7/PartitionedCall?
pool_5/PartitionedCallPartitionedCalldrop_5/PartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_5_layer_call_and_return_conditional_losses_3770612
pool_5/PartitionedCall?
pool_3/PartitionedCallPartitionedCalldrop_3/PartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_3_layer_call_and_return_conditional_losses_3770462
pool_3/PartitionedCall?
flatten_3/PartitionedCallPartitionedCallpool_3/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_3_layer_call_and_return_conditional_losses_3772012
flatten_3/PartitionedCall?
flatten_5/PartitionedCallPartitionedCallpool_5/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_5_layer_call_and_return_conditional_losses_3772152
flatten_5/PartitionedCall?
flatten_7/PartitionedCallPartitionedCallpool_7/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_7_layer_call_and_return_conditional_losses_3772292
flatten_7/PartitionedCall?
concatenate/PartitionedCallPartitionedCall"flatten_3/PartitionedCall:output:0"flatten_5/PartitionedCall:output:0"flatten_7/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*P
fKRI
G__inference_concatenate_layer_call_and_return_conditional_losses_3772452
concatenate/PartitionedCall?
dense_0/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_0_377644dense_0_377646*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_0_layer_call_and_return_conditional_losses_3772662!
dense_0/StatefulPartitionedCall?
drop_d0/PartitionedCallPartitionedCall(dense_0/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d0_layer_call_and_return_conditional_losses_3772992
drop_d0/PartitionedCall?
concatenate_1/PartitionedCallPartitionedCall drop_d0/PartitionedCall:output:0inputs_1*
Tin
2*
Tout
2*'
_output_shapes
:?????????Q* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*R
fMRK
I__inference_concatenate_1_layer_call_and_return_conditional_losses_3773192
concatenate_1/PartitionedCall?
dense_1/StatefulPartitionedCallStatefulPartitionedCall&concatenate_1/PartitionedCall:output:0dense_1_377651dense_1_377653*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_1_layer_call_and_return_conditional_losses_3773392!
dense_1/StatefulPartitionedCall?
drop_d1/PartitionedCallPartitionedCall(dense_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d1_layer_call_and_return_conditional_losses_3773722
drop_d1/PartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCall drop_d1/PartitionedCall:output:0dense_2_377657dense_2_377659*
Tin
2*
Tout
2*'
_output_shapes
:?????????<*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_2_layer_call_and_return_conditional_losses_3773962!
dense_2/StatefulPartitionedCall?
drop_d2/PartitionedCallPartitionedCall(dense_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*'
_output_shapes
:?????????<* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d2_layer_call_and_return_conditional_losses_3774292
drop_d2/PartitionedCall?
output/StatefulPartitionedCallStatefulPartitionedCall drop_d2/PartitionedCall:output:0output_377663output_377665*
Tin
2*
Tout
2*'
_output_shapes
:?????????*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_output_layer_call_and_return_conditional_losses_3774522 
output/StatefulPartitionedCall?
IdentityIdentity'output/StatefulPartitionedCall:output:0^conv_3/StatefulPartitionedCall^conv_5/StatefulPartitionedCall^conv_7/StatefulPartitionedCall ^dense_0/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall ^dense_2/StatefulPartitionedCall^output/StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????::::::::::::::2@
conv_3/StatefulPartitionedCallconv_3/StatefulPartitionedCall2@
conv_5/StatefulPartitionedCallconv_5/StatefulPartitionedCall2@
conv_7/StatefulPartitionedCallconv_7/StatefulPartitionedCall2B
dense_0/StatefulPartitionedCalldense_0/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2@
output/StatefulPartitionedCalloutput/StatefulPartitionedCall:S O
+
_output_shapes
:?????????
 
_user_specified_nameinputs:OK
'
_output_shapes
:?????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?T
?
A__inference_model_layer_call_and_return_conditional_losses_377581

inputs
inputs_1
conv_7_377531
conv_7_377533
conv_5_377536
conv_5_377538
conv_3_377541
conv_3_377543
dense_0_377556
dense_0_377558
dense_1_377563
dense_1_377565
dense_2_377569
dense_2_377571
output_377575
output_377577
identity??conv_3/StatefulPartitionedCall?conv_5/StatefulPartitionedCall?conv_7/StatefulPartitionedCall?dense_0/StatefulPartitionedCall?dense_1/StatefulPartitionedCall?dense_2/StatefulPartitionedCall?drop_3/StatefulPartitionedCall?drop_5/StatefulPartitionedCall?drop_7/StatefulPartitionedCall?drop_d0/StatefulPartitionedCall?drop_d1/StatefulPartitionedCall?drop_d2/StatefulPartitionedCall?output/StatefulPartitionedCall?
conv_7/StatefulPartitionedCallStatefulPartitionedCallinputsconv_7_377531conv_7_377533*
Tin
2*
Tout
2*+
_output_shapes
:?????????(*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_7_layer_call_and_return_conditional_losses_3770272 
conv_7/StatefulPartitionedCall?
conv_5/StatefulPartitionedCallStatefulPartitionedCallinputsconv_5_377536conv_5_377538*
Tin
2*
Tout
2*+
_output_shapes
:?????????F*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_5_layer_call_and_return_conditional_losses_3770002 
conv_5/StatefulPartitionedCall?
conv_3/StatefulPartitionedCallStatefulPartitionedCallinputsconv_3_377541conv_3_377543*
Tin
2*
Tout
2*+
_output_shapes
:?????????d*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_3_layer_call_and_return_conditional_losses_3769732 
conv_3/StatefulPartitionedCall?
drop_7/StatefulPartitionedCallStatefulPartitionedCall'conv_7/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_7_layer_call_and_return_conditional_losses_3771142 
drop_7/StatefulPartitionedCall?
drop_5/StatefulPartitionedCallStatefulPartitionedCall'conv_5/StatefulPartitionedCall:output:0^drop_7/StatefulPartitionedCall*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_5_layer_call_and_return_conditional_losses_3771442 
drop_5/StatefulPartitionedCall?
drop_3/StatefulPartitionedCallStatefulPartitionedCall'conv_3/StatefulPartitionedCall:output:0^drop_5/StatefulPartitionedCall*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_3_layer_call_and_return_conditional_losses_3771742 
drop_3/StatefulPartitionedCall?
pool_7/PartitionedCallPartitionedCall'drop_7/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_7_layer_call_and_return_conditional_losses_3770762
pool_7/PartitionedCall?
pool_5/PartitionedCallPartitionedCall'drop_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_5_layer_call_and_return_conditional_losses_3770612
pool_5/PartitionedCall?
pool_3/PartitionedCallPartitionedCall'drop_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_3_layer_call_and_return_conditional_losses_3770462
pool_3/PartitionedCall?
flatten_3/PartitionedCallPartitionedCallpool_3/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_3_layer_call_and_return_conditional_losses_3772012
flatten_3/PartitionedCall?
flatten_5/PartitionedCallPartitionedCallpool_5/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_5_layer_call_and_return_conditional_losses_3772152
flatten_5/PartitionedCall?
flatten_7/PartitionedCallPartitionedCallpool_7/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_7_layer_call_and_return_conditional_losses_3772292
flatten_7/PartitionedCall?
concatenate/PartitionedCallPartitionedCall"flatten_3/PartitionedCall:output:0"flatten_5/PartitionedCall:output:0"flatten_7/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*P
fKRI
G__inference_concatenate_layer_call_and_return_conditional_losses_3772452
concatenate/PartitionedCall?
dense_0/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_0_377556dense_0_377558*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_0_layer_call_and_return_conditional_losses_3772662!
dense_0/StatefulPartitionedCall?
drop_d0/StatefulPartitionedCallStatefulPartitionedCall(dense_0/StatefulPartitionedCall:output:0^drop_3/StatefulPartitionedCall*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d0_layer_call_and_return_conditional_losses_3772942!
drop_d0/StatefulPartitionedCall?
concatenate_1/PartitionedCallPartitionedCall(drop_d0/StatefulPartitionedCall:output:0inputs_1*
Tin
2*
Tout
2*'
_output_shapes
:?????????Q* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*R
fMRK
I__inference_concatenate_1_layer_call_and_return_conditional_losses_3773192
concatenate_1/PartitionedCall?
dense_1/StatefulPartitionedCallStatefulPartitionedCall&concatenate_1/PartitionedCall:output:0dense_1_377563dense_1_377565*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_1_layer_call_and_return_conditional_losses_3773392!
dense_1/StatefulPartitionedCall?
drop_d1/StatefulPartitionedCallStatefulPartitionedCall(dense_1/StatefulPartitionedCall:output:0 ^drop_d0/StatefulPartitionedCall*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d1_layer_call_and_return_conditional_losses_3773672!
drop_d1/StatefulPartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCall(drop_d1/StatefulPartitionedCall:output:0dense_2_377569dense_2_377571*
Tin
2*
Tout
2*'
_output_shapes
:?????????<*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_2_layer_call_and_return_conditional_losses_3773962!
dense_2/StatefulPartitionedCall?
drop_d2/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0 ^drop_d1/StatefulPartitionedCall*
Tin
2*
Tout
2*'
_output_shapes
:?????????<* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d2_layer_call_and_return_conditional_losses_3774242!
drop_d2/StatefulPartitionedCall?
output/StatefulPartitionedCallStatefulPartitionedCall(drop_d2/StatefulPartitionedCall:output:0output_377575output_377577*
Tin
2*
Tout
2*'
_output_shapes
:?????????*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_output_layer_call_and_return_conditional_losses_3774522 
output/StatefulPartitionedCall?
IdentityIdentity'output/StatefulPartitionedCall:output:0^conv_3/StatefulPartitionedCall^conv_5/StatefulPartitionedCall^conv_7/StatefulPartitionedCall ^dense_0/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall ^dense_2/StatefulPartitionedCall^drop_3/StatefulPartitionedCall^drop_5/StatefulPartitionedCall^drop_7/StatefulPartitionedCall ^drop_d0/StatefulPartitionedCall ^drop_d1/StatefulPartitionedCall ^drop_d2/StatefulPartitionedCall^output/StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????::::::::::::::2@
conv_3/StatefulPartitionedCallconv_3/StatefulPartitionedCall2@
conv_5/StatefulPartitionedCallconv_5/StatefulPartitionedCall2@
conv_7/StatefulPartitionedCallconv_7/StatefulPartitionedCall2B
dense_0/StatefulPartitionedCalldense_0/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2@
drop_3/StatefulPartitionedCalldrop_3/StatefulPartitionedCall2@
drop_5/StatefulPartitionedCalldrop_5/StatefulPartitionedCall2@
drop_7/StatefulPartitionedCalldrop_7/StatefulPartitionedCall2B
drop_d0/StatefulPartitionedCalldrop_d0/StatefulPartitionedCall2B
drop_d1/StatefulPartitionedCalldrop_d1/StatefulPartitionedCall2B
drop_d2/StatefulPartitionedCalldrop_d2/StatefulPartitionedCall2@
output/StatefulPartitionedCalloutput/StatefulPartitionedCall:S O
+
_output_shapes
:?????????
 
_user_specified_nameinputs:OK
'
_output_shapes
:?????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
a
E__inference_flatten_5_layer_call_and_return_conditional_losses_378144

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????F:S O
+
_output_shapes
:?????????F
 
_user_specified_nameinputs
?
`
B__inference_drop_7_layer_call_and_return_conditional_losses_378117

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:?????????(2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:?????????(2

Identity_1"!

identity_1Identity_1:output:0**
_input_shapes
:?????????(:S O
+
_output_shapes
:?????????(
 
_user_specified_nameinputs
?
a
B__inference_drop_7_layer_call_and_return_conditional_losses_378112

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:?????????(2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:?????????(*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:?????????(2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:?????????(2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:?????????(2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:?????????(2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????(:S O
+
_output_shapes
:?????????(
 
_user_specified_nameinputs
?
a
(__inference_drop_d1_layer_call_fn_378277

inputs
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d1_layer_call_and_return_conditional_losses_3773672
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????P22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
a
C__inference_drop_d0_layer_call_and_return_conditional_losses_377299

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:?????????P2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:?????????P2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
F
*__inference_flatten_5_layer_call_fn_378149

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_5_layer_call_and_return_conditional_losses_3772152
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????F:S O
+
_output_shapes
:?????????F
 
_user_specified_nameinputs
?p
?
__inference__traced_save_378535
file_prefix,
(savev2_conv_3_kernel_read_readvariableop*
&savev2_conv_3_bias_read_readvariableop,
(savev2_conv_5_kernel_read_readvariableop*
&savev2_conv_5_bias_read_readvariableop,
(savev2_conv_7_kernel_read_readvariableop*
&savev2_conv_7_bias_read_readvariableop-
)savev2_dense_0_kernel_read_readvariableop+
'savev2_dense_0_bias_read_readvariableop-
)savev2_dense_1_kernel_read_readvariableop+
'savev2_dense_1_bias_read_readvariableop-
)savev2_dense_2_kernel_read_readvariableop+
'savev2_dense_2_bias_read_readvariableop,
(savev2_output_kernel_read_readvariableop*
&savev2_output_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop&
"savev2_total_2_read_readvariableop&
"savev2_count_2_read_readvariableop3
/savev2_adam_conv_3_kernel_m_read_readvariableop1
-savev2_adam_conv_3_bias_m_read_readvariableop3
/savev2_adam_conv_5_kernel_m_read_readvariableop1
-savev2_adam_conv_5_bias_m_read_readvariableop3
/savev2_adam_conv_7_kernel_m_read_readvariableop1
-savev2_adam_conv_7_bias_m_read_readvariableop4
0savev2_adam_dense_0_kernel_m_read_readvariableop2
.savev2_adam_dense_0_bias_m_read_readvariableop4
0savev2_adam_dense_1_kernel_m_read_readvariableop2
.savev2_adam_dense_1_bias_m_read_readvariableop4
0savev2_adam_dense_2_kernel_m_read_readvariableop2
.savev2_adam_dense_2_bias_m_read_readvariableop3
/savev2_adam_output_kernel_m_read_readvariableop1
-savev2_adam_output_bias_m_read_readvariableop3
/savev2_adam_conv_3_kernel_v_read_readvariableop1
-savev2_adam_conv_3_bias_v_read_readvariableop3
/savev2_adam_conv_5_kernel_v_read_readvariableop1
-savev2_adam_conv_5_bias_v_read_readvariableop3
/savev2_adam_conv_7_kernel_v_read_readvariableop1
-savev2_adam_conv_7_bias_v_read_readvariableop4
0savev2_adam_dense_0_kernel_v_read_readvariableop2
.savev2_adam_dense_0_bias_v_read_readvariableop4
0savev2_adam_dense_1_kernel_v_read_readvariableop2
.savev2_adam_dense_1_bias_v_read_readvariableop4
0savev2_adam_dense_2_kernel_v_read_readvariableop2
.savev2_adam_dense_2_bias_v_read_readvariableop3
/savev2_adam_output_kernel_v_read_readvariableop1
-savev2_adam_output_bias_v_read_readvariableop
savev2_1_const

identity_1??MergeV2Checkpoints?SaveV2?SaveV2_1?
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
Const?
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*<
value3B1 B+_temp_a1b4ae199ef94fefa024e1e54203a7c0/part2	
Const_1?
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shard?
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename?
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:5*
dtype0*?
value?B?5B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE2
SaveV2/tensor_names?
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:5*
dtype0*}
valuetBr5B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slices?
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0(savev2_conv_3_kernel_read_readvariableop&savev2_conv_3_bias_read_readvariableop(savev2_conv_5_kernel_read_readvariableop&savev2_conv_5_bias_read_readvariableop(savev2_conv_7_kernel_read_readvariableop&savev2_conv_7_bias_read_readvariableop)savev2_dense_0_kernel_read_readvariableop'savev2_dense_0_bias_read_readvariableop)savev2_dense_1_kernel_read_readvariableop'savev2_dense_1_bias_read_readvariableop)savev2_dense_2_kernel_read_readvariableop'savev2_dense_2_bias_read_readvariableop(savev2_output_kernel_read_readvariableop&savev2_output_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop"savev2_total_2_read_readvariableop"savev2_count_2_read_readvariableop/savev2_adam_conv_3_kernel_m_read_readvariableop-savev2_adam_conv_3_bias_m_read_readvariableop/savev2_adam_conv_5_kernel_m_read_readvariableop-savev2_adam_conv_5_bias_m_read_readvariableop/savev2_adam_conv_7_kernel_m_read_readvariableop-savev2_adam_conv_7_bias_m_read_readvariableop0savev2_adam_dense_0_kernel_m_read_readvariableop.savev2_adam_dense_0_bias_m_read_readvariableop0savev2_adam_dense_1_kernel_m_read_readvariableop.savev2_adam_dense_1_bias_m_read_readvariableop0savev2_adam_dense_2_kernel_m_read_readvariableop.savev2_adam_dense_2_bias_m_read_readvariableop/savev2_adam_output_kernel_m_read_readvariableop-savev2_adam_output_bias_m_read_readvariableop/savev2_adam_conv_3_kernel_v_read_readvariableop-savev2_adam_conv_3_bias_v_read_readvariableop/savev2_adam_conv_5_kernel_v_read_readvariableop-savev2_adam_conv_5_bias_v_read_readvariableop/savev2_adam_conv_7_kernel_v_read_readvariableop-savev2_adam_conv_7_bias_v_read_readvariableop0savev2_adam_dense_0_kernel_v_read_readvariableop.savev2_adam_dense_0_bias_v_read_readvariableop0savev2_adam_dense_1_kernel_v_read_readvariableop.savev2_adam_dense_1_bias_v_read_readvariableop0savev2_adam_dense_2_kernel_v_read_readvariableop.savev2_adam_dense_2_bias_v_read_readvariableop/savev2_adam_output_kernel_v_read_readvariableop-savev2_adam_output_bias_v_read_readvariableop"/device:CPU:0*
_output_shapes
 *C
dtypes9
725	2
SaveV2?
ShardedFilename_1/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B :2
ShardedFilename_1/shard?
ShardedFilename_1ShardedFilenameStringJoin:output:0 ShardedFilename_1/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename_1?
SaveV2_1/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2_1/tensor_names?
SaveV2_1/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*
valueB
B 2
SaveV2_1/shape_and_slices?
SaveV2_1SaveV2ShardedFilename_1:filename:0SaveV2_1/tensor_names:output:0"SaveV2_1/shape_and_slices:output:0savev2_1_const^SaveV2"/device:CPU:0*
_output_shapes
 *
dtypes
22

SaveV2_1?
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0ShardedFilename_1:filename:0^SaveV2	^SaveV2_1"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixes?
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix	^SaveV2_1"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identity?

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints^SaveV2	^SaveV2_1*
T0*
_output_shapes
: 2

Identity_1"!

identity_1Identity_1:output:0*?
_input_shapes?
?: :d:d:F:F:(:(:	?P:P:QP:P:P<:<:<:: : : : : : : : : : : :d:d:F:F:(:(:	?P:P:QP:P:P<:<:<::d:d:F:F:(:(:	?P:P:QP:P:P<:<:<:: 2(
MergeV2CheckpointsMergeV2Checkpoints2
SaveV2SaveV22
SaveV2_1SaveV2_1:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:($
"
_output_shapes
:d: 

_output_shapes
:d:($
"
_output_shapes
:F: 

_output_shapes
:F:($
"
_output_shapes
:(: 

_output_shapes
:(:%!

_output_shapes
:	?P: 

_output_shapes
:P:$	 

_output_shapes

:QP: 


_output_shapes
:P:$ 

_output_shapes

:P<: 

_output_shapes
:<:$ 

_output_shapes

:<: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :($
"
_output_shapes
:d: 

_output_shapes
:d:($
"
_output_shapes
:F: 

_output_shapes
:F:($
"
_output_shapes
:(: 

_output_shapes
:(:% !

_output_shapes
:	?P: !

_output_shapes
:P:$" 

_output_shapes

:QP: #

_output_shapes
:P:$$ 

_output_shapes

:P<: %

_output_shapes
:<:$& 

_output_shapes

:<: '

_output_shapes
::(($
"
_output_shapes
:d: )

_output_shapes
:d:(*$
"
_output_shapes
:F: +

_output_shapes
:F:(,$
"
_output_shapes
:(: -

_output_shapes
:(:%.!

_output_shapes
:	?P: /

_output_shapes
:P:$0 

_output_shapes

:QP: 1

_output_shapes
:P:$2 

_output_shapes

:P<: 3

_output_shapes
:<:$4 

_output_shapes

:<: 5

_output_shapes
::6

_output_shapes
: 
?
?
C__inference_dense_0_layer_call_and_return_conditional_losses_378186

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?P*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:P*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*/
_input_shapes
:??????????:::P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
a
C__inference_drop_d2_layer_call_and_return_conditional_losses_377429

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:?????????<2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:?????????<2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:?????????<:O K
'
_output_shapes
:?????????<
 
_user_specified_nameinputs
?
?
B__inference_conv_3_layer_call_and_return_conditional_losses_376973

inputs/
+conv1d_expanddims_1_readvariableop_resource#
biasadd_readvariableop_resource
identity?p
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
conv1d/ExpandDims/dim?
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*8
_output_shapes&
$:"??????????????????2
conv1d/ExpandDims?
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:d*
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dim?
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:d2
conv1d/ExpandDims_1?
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*8
_output_shapes&
$:"??????????????????d*
paddingVALID*
strides
2
conv1d?
conv1d/SqueezeSqueezeconv1d:output:0*
T0*4
_output_shapes"
 :??????????????????d*
squeeze_dims
2
conv1d/Squeeze?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:d*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*4
_output_shapes"
 :??????????????????d2	
BiasAdde
ReluReluBiasAdd:output:0*
T0*4
_output_shapes"
 :??????????????????d2
Relus
IdentityIdentityRelu:activations:0*
T0*4
_output_shapes"
 :??????????????????d2

Identity"
identityIdentity:output:0*;
_input_shapes*
(:??????????????????:::\ X
4
_output_shapes"
 :??????????????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
?
G__inference_concatenate_layer_call_and_return_conditional_losses_378168
inputs_0
inputs_1
inputs_2
identity\
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concat/axis?
concatConcatV2inputs_0inputs_1inputs_2concat/axis:output:0*
N*
T0*(
_output_shapes
:??????????2
concatd
IdentityIdentityconcat:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*O
_input_shapes>
<:??????????
:??????????:??????????:R N
(
_output_shapes
:??????????

"
_user_specified_name
inputs/0:RN
(
_output_shapes
:??????????
"
_user_specified_name
inputs/1:RN
(
_output_shapes
:??????????
"
_user_specified_name
inputs/2
բ
?
A__inference_model_layer_call_and_return_conditional_losses_377882
inputs_0
inputs_16
2conv_7_conv1d_expanddims_1_readvariableop_resource*
&conv_7_biasadd_readvariableop_resource6
2conv_5_conv1d_expanddims_1_readvariableop_resource*
&conv_5_biasadd_readvariableop_resource6
2conv_3_conv1d_expanddims_1_readvariableop_resource*
&conv_3_biasadd_readvariableop_resource*
&dense_0_matmul_readvariableop_resource+
'dense_0_biasadd_readvariableop_resource*
&dense_1_matmul_readvariableop_resource+
'dense_1_biasadd_readvariableop_resource*
&dense_2_matmul_readvariableop_resource+
'dense_2_biasadd_readvariableop_resource)
%output_matmul_readvariableop_resource*
&output_biasadd_readvariableop_resource
identity?~
conv_7/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
conv_7/conv1d/ExpandDims/dim?
conv_7/conv1d/ExpandDims
ExpandDimsinputs_0%conv_7/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????2
conv_7/conv1d/ExpandDims?
)conv_7/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_7_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:(*
dtype02+
)conv_7/conv1d/ExpandDims_1/ReadVariableOp?
conv_7/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
conv_7/conv1d/ExpandDims_1/dim?
conv_7/conv1d/ExpandDims_1
ExpandDims1conv_7/conv1d/ExpandDims_1/ReadVariableOp:value:0'conv_7/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:(2
conv_7/conv1d/ExpandDims_1?
conv_7/conv1dConv2D!conv_7/conv1d/ExpandDims:output:0#conv_7/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:?????????(*
paddingVALID*
strides
2
conv_7/conv1d?
conv_7/conv1d/SqueezeSqueezeconv_7/conv1d:output:0*
T0*+
_output_shapes
:?????????(*
squeeze_dims
2
conv_7/conv1d/Squeeze?
conv_7/BiasAdd/ReadVariableOpReadVariableOp&conv_7_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
conv_7/BiasAdd/ReadVariableOp?
conv_7/BiasAddBiasAddconv_7/conv1d/Squeeze:output:0%conv_7/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:?????????(2
conv_7/BiasAddq
conv_7/ReluReluconv_7/BiasAdd:output:0*
T0*+
_output_shapes
:?????????(2
conv_7/Relu~
conv_5/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
conv_5/conv1d/ExpandDims/dim?
conv_5/conv1d/ExpandDims
ExpandDimsinputs_0%conv_5/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????2
conv_5/conv1d/ExpandDims?
)conv_5/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_5_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:F*
dtype02+
)conv_5/conv1d/ExpandDims_1/ReadVariableOp?
conv_5/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
conv_5/conv1d/ExpandDims_1/dim?
conv_5/conv1d/ExpandDims_1
ExpandDims1conv_5/conv1d/ExpandDims_1/ReadVariableOp:value:0'conv_5/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:F2
conv_5/conv1d/ExpandDims_1?
conv_5/conv1dConv2D!conv_5/conv1d/ExpandDims:output:0#conv_5/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:?????????F*
paddingVALID*
strides
2
conv_5/conv1d?
conv_5/conv1d/SqueezeSqueezeconv_5/conv1d:output:0*
T0*+
_output_shapes
:?????????F*
squeeze_dims
2
conv_5/conv1d/Squeeze?
conv_5/BiasAdd/ReadVariableOpReadVariableOp&conv_5_biasadd_readvariableop_resource*
_output_shapes
:F*
dtype02
conv_5/BiasAdd/ReadVariableOp?
conv_5/BiasAddBiasAddconv_5/conv1d/Squeeze:output:0%conv_5/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:?????????F2
conv_5/BiasAddq
conv_5/ReluReluconv_5/BiasAdd:output:0*
T0*+
_output_shapes
:?????????F2
conv_5/Relu~
conv_3/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
conv_3/conv1d/ExpandDims/dim?
conv_3/conv1d/ExpandDims
ExpandDimsinputs_0%conv_3/conv1d/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????2
conv_3/conv1d/ExpandDims?
)conv_3/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv_3_conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:d*
dtype02+
)conv_3/conv1d/ExpandDims_1/ReadVariableOp?
conv_3/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
conv_3/conv1d/ExpandDims_1/dim?
conv_3/conv1d/ExpandDims_1
ExpandDims1conv_3/conv1d/ExpandDims_1/ReadVariableOp:value:0'conv_3/conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:d2
conv_3/conv1d/ExpandDims_1?
conv_3/conv1dConv2D!conv_3/conv1d/ExpandDims:output:0#conv_3/conv1d/ExpandDims_1:output:0*
T0*/
_output_shapes
:?????????d*
paddingVALID*
strides
2
conv_3/conv1d?
conv_3/conv1d/SqueezeSqueezeconv_3/conv1d:output:0*
T0*+
_output_shapes
:?????????d*
squeeze_dims
2
conv_3/conv1d/Squeeze?
conv_3/BiasAdd/ReadVariableOpReadVariableOp&conv_3_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02
conv_3/BiasAdd/ReadVariableOp?
conv_3/BiasAddBiasAddconv_3/conv1d/Squeeze:output:0%conv_3/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:?????????d2
conv_3/BiasAddq
conv_3/ReluReluconv_3/BiasAdd:output:0*
T0*+
_output_shapes
:?????????d2
conv_3/Reluq
drop_7/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
drop_7/dropout/Const?
drop_7/dropout/MulMulconv_7/Relu:activations:0drop_7/dropout/Const:output:0*
T0*+
_output_shapes
:?????????(2
drop_7/dropout/Mulu
drop_7/dropout/ShapeShapeconv_7/Relu:activations:0*
T0*
_output_shapes
:2
drop_7/dropout/Shape?
+drop_7/dropout/random_uniform/RandomUniformRandomUniformdrop_7/dropout/Shape:output:0*
T0*+
_output_shapes
:?????????(*
dtype0*
seed????2-
+drop_7/dropout/random_uniform/RandomUniform?
drop_7/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
drop_7/dropout/GreaterEqual/y?
drop_7/dropout/GreaterEqualGreaterEqual4drop_7/dropout/random_uniform/RandomUniform:output:0&drop_7/dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:?????????(2
drop_7/dropout/GreaterEqual?
drop_7/dropout/CastCastdrop_7/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:?????????(2
drop_7/dropout/Cast?
drop_7/dropout/Mul_1Muldrop_7/dropout/Mul:z:0drop_7/dropout/Cast:y:0*
T0*+
_output_shapes
:?????????(2
drop_7/dropout/Mul_1q
drop_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
drop_5/dropout/Const?
drop_5/dropout/MulMulconv_5/Relu:activations:0drop_5/dropout/Const:output:0*
T0*+
_output_shapes
:?????????F2
drop_5/dropout/Mulu
drop_5/dropout/ShapeShapeconv_5/Relu:activations:0*
T0*
_output_shapes
:2
drop_5/dropout/Shape?
+drop_5/dropout/random_uniform/RandomUniformRandomUniformdrop_5/dropout/Shape:output:0*
T0*+
_output_shapes
:?????????F*
dtype0*
seed????*
seed22-
+drop_5/dropout/random_uniform/RandomUniform?
drop_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
drop_5/dropout/GreaterEqual/y?
drop_5/dropout/GreaterEqualGreaterEqual4drop_5/dropout/random_uniform/RandomUniform:output:0&drop_5/dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:?????????F2
drop_5/dropout/GreaterEqual?
drop_5/dropout/CastCastdrop_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:?????????F2
drop_5/dropout/Cast?
drop_5/dropout/Mul_1Muldrop_5/dropout/Mul:z:0drop_5/dropout/Cast:y:0*
T0*+
_output_shapes
:?????????F2
drop_5/dropout/Mul_1q
drop_3/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
drop_3/dropout/Const?
drop_3/dropout/MulMulconv_3/Relu:activations:0drop_3/dropout/Const:output:0*
T0*+
_output_shapes
:?????????d2
drop_3/dropout/Mulu
drop_3/dropout/ShapeShapeconv_3/Relu:activations:0*
T0*
_output_shapes
:2
drop_3/dropout/Shape?
+drop_3/dropout/random_uniform/RandomUniformRandomUniformdrop_3/dropout/Shape:output:0*
T0*+
_output_shapes
:?????????d*
dtype0*
seed????*
seed22-
+drop_3/dropout/random_uniform/RandomUniform?
drop_3/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
drop_3/dropout/GreaterEqual/y?
drop_3/dropout/GreaterEqualGreaterEqual4drop_3/dropout/random_uniform/RandomUniform:output:0&drop_3/dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:?????????d2
drop_3/dropout/GreaterEqual?
drop_3/dropout/CastCastdrop_3/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:?????????d2
drop_3/dropout/Cast?
drop_3/dropout/Mul_1Muldrop_3/dropout/Mul:z:0drop_3/dropout/Cast:y:0*
T0*+
_output_shapes
:?????????d2
drop_3/dropout/Mul_1p
pool_7/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
pool_7/ExpandDims/dim?
pool_7/ExpandDims
ExpandDimsdrop_7/dropout/Mul_1:z:0pool_7/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????(2
pool_7/ExpandDims?
pool_7/AvgPoolAvgPoolpool_7/ExpandDims:output:0*
T0*/
_output_shapes
:?????????(*
ksize
*
paddingSAME*
strides
2
pool_7/AvgPool?
pool_7/SqueezeSqueezepool_7/AvgPool:output:0*
T0*+
_output_shapes
:?????????(*
squeeze_dims
2
pool_7/Squeezep
pool_5/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
pool_5/ExpandDims/dim?
pool_5/ExpandDims
ExpandDimsdrop_5/dropout/Mul_1:z:0pool_5/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????F2
pool_5/ExpandDims?
pool_5/AvgPoolAvgPoolpool_5/ExpandDims:output:0*
T0*/
_output_shapes
:?????????F*
ksize
*
paddingSAME*
strides
2
pool_5/AvgPool?
pool_5/SqueezeSqueezepool_5/AvgPool:output:0*
T0*+
_output_shapes
:?????????F*
squeeze_dims
2
pool_5/Squeezep
pool_3/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
pool_3/ExpandDims/dim?
pool_3/ExpandDims
ExpandDimsdrop_3/dropout/Mul_1:z:0pool_3/ExpandDims/dim:output:0*
T0*/
_output_shapes
:?????????d2
pool_3/ExpandDims?
pool_3/AvgPoolAvgPoolpool_3/ExpandDims:output:0*
T0*/
_output_shapes
:?????????d*
ksize
*
paddingSAME*
strides
2
pool_3/AvgPool?
pool_3/SqueezeSqueezepool_3/AvgPool:output:0*
T0*+
_output_shapes
:?????????d*
squeeze_dims
2
pool_3/Squeezes
flatten_3/ConstConst*
_output_shapes
:*
dtype0*
valueB"????x  2
flatten_3/Const?
flatten_3/ReshapeReshapepool_3/Squeeze:output:0flatten_3/Const:output:0*
T0*(
_output_shapes
:??????????
2
flatten_3/Reshapes
flatten_5/ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
flatten_5/Const?
flatten_5/ReshapeReshapepool_5/Squeeze:output:0flatten_5/Const:output:0*
T0*(
_output_shapes
:??????????2
flatten_5/Reshapes
flatten_7/ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
flatten_7/Const?
flatten_7/ReshapeReshapepool_7/Squeeze:output:0flatten_7/Const:output:0*
T0*(
_output_shapes
:??????????2
flatten_7/Reshapet
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate/concat/axis?
concatenate/concatConcatV2flatten_3/Reshape:output:0flatten_5/Reshape:output:0flatten_7/Reshape:output:0 concatenate/concat/axis:output:0*
N*
T0*(
_output_shapes
:??????????2
concatenate/concat?
dense_0/MatMul/ReadVariableOpReadVariableOp&dense_0_matmul_readvariableop_resource*
_output_shapes
:	?P*
dtype02
dense_0/MatMul/ReadVariableOp?
dense_0/MatMulMatMulconcatenate/concat:output:0%dense_0/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
dense_0/MatMul?
dense_0/BiasAdd/ReadVariableOpReadVariableOp'dense_0_biasadd_readvariableop_resource*
_output_shapes
:P*
dtype02 
dense_0/BiasAdd/ReadVariableOp?
dense_0/BiasAddBiasAdddense_0/MatMul:product:0&dense_0/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
dense_0/BiasAddp
dense_0/ReluReludense_0/BiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
dense_0/Relus
drop_d0/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
drop_d0/dropout/Const?
drop_d0/dropout/MulMuldense_0/Relu:activations:0drop_d0/dropout/Const:output:0*
T0*'
_output_shapes
:?????????P2
drop_d0/dropout/Mulx
drop_d0/dropout/ShapeShapedense_0/Relu:activations:0*
T0*
_output_shapes
:2
drop_d0/dropout/Shape?
,drop_d0/dropout/random_uniform/RandomUniformRandomUniformdrop_d0/dropout/Shape:output:0*
T0*'
_output_shapes
:?????????P*
dtype0*
seed????*
seed22.
,drop_d0/dropout/random_uniform/RandomUniform?
drop_d0/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2 
drop_d0/dropout/GreaterEqual/y?
drop_d0/dropout/GreaterEqualGreaterEqual5drop_d0/dropout/random_uniform/RandomUniform:output:0'drop_d0/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:?????????P2
drop_d0/dropout/GreaterEqual?
drop_d0/dropout/CastCast drop_d0/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:?????????P2
drop_d0/dropout/Cast?
drop_d0/dropout/Mul_1Muldrop_d0/dropout/Mul:z:0drop_d0/dropout/Cast:y:0*
T0*'
_output_shapes
:?????????P2
drop_d0/dropout/Mul_1x
concatenate_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate_1/concat/axis?
concatenate_1/concatConcatV2drop_d0/dropout/Mul_1:z:0inputs_1"concatenate_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????Q2
concatenate_1/concat?
dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes

:QP*
dtype02
dense_1/MatMul/ReadVariableOp?
dense_1/MatMulMatMulconcatenate_1/concat:output:0%dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
dense_1/MatMul?
dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:P*
dtype02 
dense_1/BiasAdd/ReadVariableOp?
dense_1/BiasAddBiasAdddense_1/MatMul:product:0&dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
dense_1/BiasAddp
dense_1/ReluReludense_1/BiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
dense_1/Relus
drop_d1/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
drop_d1/dropout/Const?
drop_d1/dropout/MulMuldense_1/Relu:activations:0drop_d1/dropout/Const:output:0*
T0*'
_output_shapes
:?????????P2
drop_d1/dropout/Mulx
drop_d1/dropout/ShapeShapedense_1/Relu:activations:0*
T0*
_output_shapes
:2
drop_d1/dropout/Shape?
,drop_d1/dropout/random_uniform/RandomUniformRandomUniformdrop_d1/dropout/Shape:output:0*
T0*'
_output_shapes
:?????????P*
dtype0*
seed????*
seed22.
,drop_d1/dropout/random_uniform/RandomUniform?
drop_d1/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2 
drop_d1/dropout/GreaterEqual/y?
drop_d1/dropout/GreaterEqualGreaterEqual5drop_d1/dropout/random_uniform/RandomUniform:output:0'drop_d1/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:?????????P2
drop_d1/dropout/GreaterEqual?
drop_d1/dropout/CastCast drop_d1/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:?????????P2
drop_d1/dropout/Cast?
drop_d1/dropout/Mul_1Muldrop_d1/dropout/Mul:z:0drop_d1/dropout/Cast:y:0*
T0*'
_output_shapes
:?????????P2
drop_d1/dropout/Mul_1?
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:P<*
dtype02
dense_2/MatMul/ReadVariableOp?
dense_2/MatMulMatMuldrop_d1/dropout/Mul_1:z:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2
dense_2/MatMul?
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:<*
dtype02 
dense_2/BiasAdd/ReadVariableOp?
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2
dense_2/BiasAddp
dense_2/ReluReludense_2/BiasAdd:output:0*
T0*'
_output_shapes
:?????????<2
dense_2/Relus
drop_d2/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
drop_d2/dropout/Const?
drop_d2/dropout/MulMuldense_2/Relu:activations:0drop_d2/dropout/Const:output:0*
T0*'
_output_shapes
:?????????<2
drop_d2/dropout/Mulx
drop_d2/dropout/ShapeShapedense_2/Relu:activations:0*
T0*
_output_shapes
:2
drop_d2/dropout/Shape?
,drop_d2/dropout/random_uniform/RandomUniformRandomUniformdrop_d2/dropout/Shape:output:0*
T0*'
_output_shapes
:?????????<*
dtype0*
seed????*
seed22.
,drop_d2/dropout/random_uniform/RandomUniform?
drop_d2/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2 
drop_d2/dropout/GreaterEqual/y?
drop_d2/dropout/GreaterEqualGreaterEqual5drop_d2/dropout/random_uniform/RandomUniform:output:0'drop_d2/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:?????????<2
drop_d2/dropout/GreaterEqual?
drop_d2/dropout/CastCast drop_d2/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:?????????<2
drop_d2/dropout/Cast?
drop_d2/dropout/Mul_1Muldrop_d2/dropout/Mul:z:0drop_d2/dropout/Cast:y:0*
T0*'
_output_shapes
:?????????<2
drop_d2/dropout/Mul_1?
output/MatMul/ReadVariableOpReadVariableOp%output_matmul_readvariableop_resource*
_output_shapes

:<*
dtype02
output/MatMul/ReadVariableOp?
output/MatMulMatMuldrop_d2/dropout/Mul_1:z:0$output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
output/MatMul?
output/BiasAdd/ReadVariableOpReadVariableOp&output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
output/BiasAdd/ReadVariableOp?
output/BiasAddBiasAddoutput/MatMul:product:0%output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
output/BiasAddk
IdentityIdentityoutput/BiasAdd:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????:::::::::::::::U Q
+
_output_shapes
:?????????
"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????
"
_user_specified_name
inputs/1:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
b
C__inference_drop_d2_layer_call_and_return_conditional_losses_378314

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:?????????<2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:?????????<*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:?????????<2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:?????????<2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:?????????<2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:?????????<2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????<:O K
'
_output_shapes
:?????????<
 
_user_specified_nameinputs
?U
?
A__inference_model_layer_call_and_return_conditional_losses_377469
input_onehot
	input_dgb
conv_7_377087
conv_7_377089
conv_5_377092
conv_5_377094
conv_3_377097
conv_3_377099
dense_0_377277
dense_0_377279
dense_1_377350
dense_1_377352
dense_2_377407
dense_2_377409
output_377463
output_377465
identity??conv_3/StatefulPartitionedCall?conv_5/StatefulPartitionedCall?conv_7/StatefulPartitionedCall?dense_0/StatefulPartitionedCall?dense_1/StatefulPartitionedCall?dense_2/StatefulPartitionedCall?drop_3/StatefulPartitionedCall?drop_5/StatefulPartitionedCall?drop_7/StatefulPartitionedCall?drop_d0/StatefulPartitionedCall?drop_d1/StatefulPartitionedCall?drop_d2/StatefulPartitionedCall?output/StatefulPartitionedCall?
conv_7/StatefulPartitionedCallStatefulPartitionedCallinput_onehotconv_7_377087conv_7_377089*
Tin
2*
Tout
2*+
_output_shapes
:?????????(*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_7_layer_call_and_return_conditional_losses_3770272 
conv_7/StatefulPartitionedCall?
conv_5/StatefulPartitionedCallStatefulPartitionedCallinput_onehotconv_5_377092conv_5_377094*
Tin
2*
Tout
2*+
_output_shapes
:?????????F*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_5_layer_call_and_return_conditional_losses_3770002 
conv_5/StatefulPartitionedCall?
conv_3/StatefulPartitionedCallStatefulPartitionedCallinput_onehotconv_3_377097conv_3_377099*
Tin
2*
Tout
2*+
_output_shapes
:?????????d*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_3_layer_call_and_return_conditional_losses_3769732 
conv_3/StatefulPartitionedCall?
drop_7/StatefulPartitionedCallStatefulPartitionedCall'conv_7/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_7_layer_call_and_return_conditional_losses_3771142 
drop_7/StatefulPartitionedCall?
drop_5/StatefulPartitionedCallStatefulPartitionedCall'conv_5/StatefulPartitionedCall:output:0^drop_7/StatefulPartitionedCall*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_5_layer_call_and_return_conditional_losses_3771442 
drop_5/StatefulPartitionedCall?
drop_3/StatefulPartitionedCallStatefulPartitionedCall'conv_3/StatefulPartitionedCall:output:0^drop_5/StatefulPartitionedCall*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_3_layer_call_and_return_conditional_losses_3771742 
drop_3/StatefulPartitionedCall?
pool_7/PartitionedCallPartitionedCall'drop_7/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_7_layer_call_and_return_conditional_losses_3770762
pool_7/PartitionedCall?
pool_5/PartitionedCallPartitionedCall'drop_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_5_layer_call_and_return_conditional_losses_3770612
pool_5/PartitionedCall?
pool_3/PartitionedCallPartitionedCall'drop_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*+
_output_shapes
:?????????d* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_3_layer_call_and_return_conditional_losses_3770462
pool_3/PartitionedCall?
flatten_3/PartitionedCallPartitionedCallpool_3/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_3_layer_call_and_return_conditional_losses_3772012
flatten_3/PartitionedCall?
flatten_5/PartitionedCallPartitionedCallpool_5/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_5_layer_call_and_return_conditional_losses_3772152
flatten_5/PartitionedCall?
flatten_7/PartitionedCallPartitionedCallpool_7/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*N
fIRG
E__inference_flatten_7_layer_call_and_return_conditional_losses_3772292
flatten_7/PartitionedCall?
concatenate/PartitionedCallPartitionedCall"flatten_3/PartitionedCall:output:0"flatten_5/PartitionedCall:output:0"flatten_7/PartitionedCall:output:0*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*P
fKRI
G__inference_concatenate_layer_call_and_return_conditional_losses_3772452
concatenate/PartitionedCall?
dense_0/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_0_377277dense_0_377279*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_0_layer_call_and_return_conditional_losses_3772662!
dense_0/StatefulPartitionedCall?
drop_d0/StatefulPartitionedCallStatefulPartitionedCall(dense_0/StatefulPartitionedCall:output:0^drop_3/StatefulPartitionedCall*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d0_layer_call_and_return_conditional_losses_3772942!
drop_d0/StatefulPartitionedCall?
concatenate_1/PartitionedCallPartitionedCall(drop_d0/StatefulPartitionedCall:output:0	input_dgb*
Tin
2*
Tout
2*'
_output_shapes
:?????????Q* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*R
fMRK
I__inference_concatenate_1_layer_call_and_return_conditional_losses_3773192
concatenate_1/PartitionedCall?
dense_1/StatefulPartitionedCallStatefulPartitionedCall&concatenate_1/PartitionedCall:output:0dense_1_377350dense_1_377352*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_1_layer_call_and_return_conditional_losses_3773392!
dense_1/StatefulPartitionedCall?
drop_d1/StatefulPartitionedCallStatefulPartitionedCall(dense_1/StatefulPartitionedCall:output:0 ^drop_d0/StatefulPartitionedCall*
Tin
2*
Tout
2*'
_output_shapes
:?????????P* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d1_layer_call_and_return_conditional_losses_3773672!
drop_d1/StatefulPartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCall(drop_d1/StatefulPartitionedCall:output:0dense_2_377407dense_2_377409*
Tin
2*
Tout
2*'
_output_shapes
:?????????<*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_2_layer_call_and_return_conditional_losses_3773962!
dense_2/StatefulPartitionedCall?
drop_d2/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0 ^drop_d1/StatefulPartitionedCall*
Tin
2*
Tout
2*'
_output_shapes
:?????????<* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d2_layer_call_and_return_conditional_losses_3774242!
drop_d2/StatefulPartitionedCall?
output/StatefulPartitionedCallStatefulPartitionedCall(drop_d2/StatefulPartitionedCall:output:0output_377463output_377465*
Tin
2*
Tout
2*'
_output_shapes
:?????????*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_output_layer_call_and_return_conditional_losses_3774522 
output/StatefulPartitionedCall?
IdentityIdentity'output/StatefulPartitionedCall:output:0^conv_3/StatefulPartitionedCall^conv_5/StatefulPartitionedCall^conv_7/StatefulPartitionedCall ^dense_0/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall ^dense_2/StatefulPartitionedCall^drop_3/StatefulPartitionedCall^drop_5/StatefulPartitionedCall^drop_7/StatefulPartitionedCall ^drop_d0/StatefulPartitionedCall ^drop_d1/StatefulPartitionedCall ^drop_d2/StatefulPartitionedCall^output/StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????::::::::::::::2@
conv_3/StatefulPartitionedCallconv_3/StatefulPartitionedCall2@
conv_5/StatefulPartitionedCallconv_5/StatefulPartitionedCall2@
conv_7/StatefulPartitionedCallconv_7/StatefulPartitionedCall2B
dense_0/StatefulPartitionedCalldense_0/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2@
drop_3/StatefulPartitionedCalldrop_3/StatefulPartitionedCall2@
drop_5/StatefulPartitionedCalldrop_5/StatefulPartitionedCall2@
drop_7/StatefulPartitionedCalldrop_7/StatefulPartitionedCall2B
drop_d0/StatefulPartitionedCalldrop_d0/StatefulPartitionedCall2B
drop_d1/StatefulPartitionedCalldrop_d1/StatefulPartitionedCall2B
drop_d2/StatefulPartitionedCalldrop_d2/StatefulPartitionedCall2@
output/StatefulPartitionedCalloutput/StatefulPartitionedCall:Y U
+
_output_shapes
:?????????
&
_user_specified_nameinput_onehot:RN
'
_output_shapes
:?????????
#
_user_specified_name	input_dGB:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
`
B__inference_drop_5_layer_call_and_return_conditional_losses_377149

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:?????????F2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:?????????F2

Identity_1"!

identity_1Identity_1:output:0**
_input_shapes
:?????????F:S O
+
_output_shapes
:?????????F
 
_user_specified_nameinputs
?
?
B__inference_output_layer_call_and_return_conditional_losses_377452

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:<*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2	
BiasAddd
IdentityIdentityBiasAdd:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????<:::O K
'
_output_shapes
:?????????<
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
C
'__inference_pool_5_layer_call_fn_377067

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*=
_output_shapes+
):'???????????????????????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_5_layer_call_and_return_conditional_losses_3770612
PartitionedCall?
IdentityIdentityPartitionedCall:output:0*
T0*=
_output_shapes+
):'???????????????????????????2

Identity"
identityIdentity:output:0*<
_input_shapes+
):'???????????????????????????:e a
=
_output_shapes+
):'???????????????????????????
 
_user_specified_nameinputs
?
D
(__inference_drop_d2_layer_call_fn_378329

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*'
_output_shapes
:?????????<* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d2_layer_call_and_return_conditional_losses_3774292
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:?????????<2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????<:O K
'
_output_shapes
:?????????<
 
_user_specified_nameinputs
?
?
C__inference_dense_2_layer_call_and_return_conditional_losses_377396

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:P<*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:<*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????<2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:?????????<2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:?????????<2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????P:::O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
C
'__inference_pool_7_layer_call_fn_377082

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*=
_output_shapes+
):'???????????????????????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_pool_7_layer_call_and_return_conditional_losses_3770762
PartitionedCall?
IdentityIdentityPartitionedCall:output:0*
T0*=
_output_shapes+
):'???????????????????????????2

Identity"
identityIdentity:output:0*<
_input_shapes+
):'???????????????????????????:e a
=
_output_shapes+
):'???????????????????????????
 
_user_specified_nameinputs
?
|
'__inference_conv_5_layer_call_fn_377010

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*4
_output_shapes"
 :??????????????????F*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_5_layer_call_and_return_conditional_losses_3770002
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*4
_output_shapes"
 :??????????????????F2

Identity"
identityIdentity:output:0*;
_input_shapes*
(:??????????????????::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
a
(__inference_drop_d2_layer_call_fn_378324

inputs
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*'
_output_shapes
:?????????<* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_drop_d2_layer_call_and_return_conditional_losses_3774242
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????<2

Identity"
identityIdentity:output:0*&
_input_shapes
:?????????<22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????<
 
_user_specified_nameinputs
?
?
B__inference_output_layer_call_and_return_conditional_losses_378339

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:<*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2	
BiasAddd
IdentityIdentityBiasAdd:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????<:::O K
'
_output_shapes
:?????????<
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
}
(__inference_dense_2_layer_call_fn_378302

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:?????????<*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_2_layer_call_and_return_conditional_losses_3773962
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????<2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????P::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
}
(__inference_dense_1_layer_call_fn_378255

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_1_layer_call_and_return_conditional_losses_3773392
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????Q::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????Q
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
f
,__inference_concatenate_layer_call_fn_378175
inputs_0
inputs_1
inputs_2
identity?
PartitionedCallPartitionedCallinputs_0inputs_1inputs_2*
Tin
2*
Tout
2*(
_output_shapes
:??????????* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*P
fKRI
G__inference_concatenate_layer_call_and_return_conditional_losses_3772452
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*O
_input_shapes>
<:??????????
:??????????:??????????:R N
(
_output_shapes
:??????????

"
_user_specified_name
inputs/0:RN
(
_output_shapes
:??????????
"
_user_specified_name
inputs/1:RN
(
_output_shapes
:??????????
"
_user_specified_name
inputs/2
?
a
E__inference_flatten_5_layer_call_and_return_conditional_losses_377215

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????F:S O
+
_output_shapes
:?????????F
 
_user_specified_nameinputs
?
|
'__inference_output_layer_call_fn_378348

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:?????????*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_output_layer_call_and_return_conditional_losses_3774522
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????<::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????<
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
C
'__inference_drop_7_layer_call_fn_378127

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*+
_output_shapes
:?????????(* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_7_layer_call_and_return_conditional_losses_3771192
PartitionedCallp
IdentityIdentityPartitionedCall:output:0*
T0*+
_output_shapes
:?????????(2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????(:S O
+
_output_shapes
:?????????(
 
_user_specified_nameinputs
?
?
B__inference_conv_5_layer_call_and_return_conditional_losses_377000

inputs/
+conv1d_expanddims_1_readvariableop_resource#
biasadd_readvariableop_resource
identity?p
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
value	B :2
conv1d/ExpandDims/dim?
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*8
_output_shapes&
$:"??????????????????2
conv1d/ExpandDims?
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*"
_output_shapes
:F*
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dim?
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*&
_output_shapes
:F2
conv1d/ExpandDims_1?
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*8
_output_shapes&
$:"??????????????????F*
paddingVALID*
strides
2
conv1d?
conv1d/SqueezeSqueezeconv1d:output:0*
T0*4
_output_shapes"
 :??????????????????F*
squeeze_dims
2
conv1d/Squeeze?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:F*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*4
_output_shapes"
 :??????????????????F2	
BiasAdde
ReluReluBiasAdd:output:0*
T0*4
_output_shapes"
 :??????????????????F2
Relus
IdentityIdentityRelu:activations:0*
T0*4
_output_shapes"
 :??????????????????F2

Identity"
identityIdentity:output:0*;
_input_shapes*
(:??????????????????:::\ X
4
_output_shapes"
 :??????????????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
?
C__inference_dense_1_layer_call_and_return_conditional_losses_378246

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:QP*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:P*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????P2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:?????????P2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????Q:::O K
'
_output_shapes
:?????????Q
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
`
B__inference_drop_7_layer_call_and_return_conditional_losses_377119

inputs

identity_1^
IdentityIdentityinputs*
T0*+
_output_shapes
:?????????(2

Identitym

Identity_1IdentityIdentity:output:0*
T0*+
_output_shapes
:?????????(2

Identity_1"!

identity_1Identity_1:output:0**
_input_shapes
:?????????(:S O
+
_output_shapes
:?????????(
 
_user_specified_nameinputs
?
s
I__inference_concatenate_1_layer_call_and_return_conditional_losses_377319

inputs
inputs_1
identity\
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concat/axis
concatConcatV2inputsinputs_1concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????Q2
concatc
IdentityIdentityconcat:output:0*
T0*'
_output_shapes
:?????????Q2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????P:?????????:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs:OK
'
_output_shapes
:?????????
 
_user_specified_nameinputs
?
`
'__inference_drop_5_layer_call_fn_378095

inputs
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*+
_output_shapes
:?????????F* 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_drop_5_layer_call_and_return_conditional_losses_3771442
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:?????????F2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????F22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:?????????F
 
_user_specified_nameinputs
?
a
B__inference_drop_7_layer_call_and_return_conditional_losses_377114

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Constw
dropout/MulMulinputsdropout/Const:output:0*
T0*+
_output_shapes
:?????????(2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*+
_output_shapes
:?????????(*
dtype0*
seed????2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*+
_output_shapes
:?????????(2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*+
_output_shapes
:?????????(2
dropout/Cast~
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*+
_output_shapes
:?????????(2
dropout/Mul_1i
IdentityIdentitydropout/Mul_1:z:0*
T0*+
_output_shapes
:?????????(2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????(:S O
+
_output_shapes
:?????????(
 
_user_specified_nameinputs
?
a
C__inference_drop_d1_layer_call_and_return_conditional_losses_377372

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:?????????P2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:?????????P2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:?????????P:O K
'
_output_shapes
:?????????P
 
_user_specified_nameinputs
?
|
'__inference_conv_7_layer_call_fn_377037

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*4
_output_shapes"
 :??????????????????(*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*K
fFRD
B__inference_conv_7_layer_call_and_return_conditional_losses_3770272
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*4
_output_shapes"
 :??????????????????(2

Identity"
identityIdentity:output:0*;
_input_shapes*
(:??????????????????::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
}
(__inference_dense_0_layer_call_fn_378195

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:?????????P*$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_dense_0_layer_call_and_return_conditional_losses_3772662
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????P2

Identity"
identityIdentity:output:0*/
_input_shapes
:??????????::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
?
?
&__inference_model_layer_call_fn_378046
inputs_0
inputs_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:?????????*0
_read_only_resource_inputs
	
**
config_proto

CPU

GPU 2J 8*J
fERC
A__inference_model_layer_call_and_return_conditional_losses_3776692
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*u
_input_shapesd
b:?????????:?????????::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:U Q
+
_output_shapes
:?????????
"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????
"
_user_specified_name
inputs/1:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
?
a
E__inference_flatten_7_layer_call_and_return_conditional_losses_378155

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"?????  2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0**
_input_shapes
:?????????(:S O
+
_output_shapes
:?????????(
 
_user_specified_nameinputs"?L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*?
serving_default?
?
	input_dGB2
serving_default_input_dGB:0?????????
I
input_onehot9
serving_default_input_onehot:0?????????:
output0
StatefulPartitionedCall:0?????????tensorflow/serving/predict:??
ǋ
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer-4
layer-5
layer-6
layer-7
	layer-8

layer-9
layer-10
layer-11
layer-12
layer-13
layer_with_weights-3
layer-14
layer-15
layer-16
layer-17
layer_with_weights-4
layer-18
layer-19
layer_with_weights-5
layer-20
layer-21
layer_with_weights-6
layer-22
	optimizer
regularization_losses
trainable_variables
	variables
	keras_api

signatures
?_default_save_signature
+?&call_and_return_all_conditional_losses
?__call__"??
_tf_keras_model??{"class_name": "Model", "name": "model", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "model", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_onehot"}, "name": "input_onehot", "inbound_nodes": []}, {"class_name": "Conv1D", "config": {"name": "conv_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "filters": 100, "kernel_size": {"class_name": "__tuple__", "items": [3]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv_3", "inbound_nodes": [[["input_onehot", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "conv_5", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "filters": 70, "kernel_size": {"class_name": "__tuple__", "items": [5]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv_5", "inbound_nodes": [[["input_onehot", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "conv_7", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "filters": 40, "kernel_size": {"class_name": "__tuple__", "items": [7]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv_7", "inbound_nodes": [[["input_onehot", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_3", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_3", "inbound_nodes": [[["conv_3", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_5", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_5", "inbound_nodes": [[["conv_5", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_7", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_7", "inbound_nodes": [[["conv_7", 0, 0, {}]]]}, {"class_name": "AveragePooling1D", "config": {"name": "pool_3", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "same", "data_format": "channels_last"}, "name": "pool_3", "inbound_nodes": [[["drop_3", 0, 0, {}]]]}, {"class_name": "AveragePooling1D", "config": {"name": "pool_5", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "same", "data_format": "channels_last"}, "name": "pool_5", "inbound_nodes": [[["drop_5", 0, 0, {}]]]}, {"class_name": "AveragePooling1D", "config": {"name": "pool_7", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "same", "data_format": "channels_last"}, "name": "pool_7", "inbound_nodes": [[["drop_7", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten_3", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten_3", "inbound_nodes": [[["pool_3", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten_5", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten_5", "inbound_nodes": [[["pool_5", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten_7", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten_7", "inbound_nodes": [[["pool_7", 0, 0, {}]]]}, {"class_name": "Concatenate", "config": {"name": "concatenate", "trainable": true, "dtype": "float32", "axis": -1}, "name": "concatenate", "inbound_nodes": [[["flatten_3", 0, 0, {}], ["flatten_5", 0, 0, {}], ["flatten_7", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_0", "trainable": true, "dtype": "float32", "units": 80, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_0", "inbound_nodes": [[["concatenate", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_d0", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_d0", "inbound_nodes": [[["dense_0", 0, 0, {}]]]}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_dGB"}, "name": "input_dGB", "inbound_nodes": []}, {"class_name": "Concatenate", "config": {"name": "concatenate_1", "trainable": true, "dtype": "float32", "axis": -1}, "name": "concatenate_1", "inbound_nodes": [[["drop_d0", 0, 0, {}], ["input_dGB", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 80, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_1", "inbound_nodes": [[["concatenate_1", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_d1", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_d1", "inbound_nodes": [[["dense_1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 60, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_2", "inbound_nodes": [[["drop_d1", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_d2", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_d2", "inbound_nodes": [[["dense_2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "output", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "output", "inbound_nodes": [[["drop_d2", 0, 0, {}]]]}], "input_layers": [["input_onehot", 0, 0], ["input_dGB", 0, 0]], "output_layers": [["output", 0, 0]]}, "build_input_shape": [{"class_name": "TensorShape", "items": [null, 30, 4]}, {"class_name": "TensorShape", "items": [null, 1]}], "is_graph_network": true, "keras_version": "2.3.0-tf", "backend": "tensorflow", "model_config": {"class_name": "Model", "config": {"name": "model", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_onehot"}, "name": "input_onehot", "inbound_nodes": []}, {"class_name": "Conv1D", "config": {"name": "conv_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "filters": 100, "kernel_size": {"class_name": "__tuple__", "items": [3]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv_3", "inbound_nodes": [[["input_onehot", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "conv_5", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "filters": 70, "kernel_size": {"class_name": "__tuple__", "items": [5]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv_5", "inbound_nodes": [[["input_onehot", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "conv_7", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "filters": 40, "kernel_size": {"class_name": "__tuple__", "items": [7]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv_7", "inbound_nodes": [[["input_onehot", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_3", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_3", "inbound_nodes": [[["conv_3", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_5", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_5", "inbound_nodes": [[["conv_5", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_7", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_7", "inbound_nodes": [[["conv_7", 0, 0, {}]]]}, {"class_name": "AveragePooling1D", "config": {"name": "pool_3", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "same", "data_format": "channels_last"}, "name": "pool_3", "inbound_nodes": [[["drop_3", 0, 0, {}]]]}, {"class_name": "AveragePooling1D", "config": {"name": "pool_5", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "same", "data_format": "channels_last"}, "name": "pool_5", "inbound_nodes": [[["drop_5", 0, 0, {}]]]}, {"class_name": "AveragePooling1D", "config": {"name": "pool_7", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "same", "data_format": "channels_last"}, "name": "pool_7", "inbound_nodes": [[["drop_7", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten_3", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten_3", "inbound_nodes": [[["pool_3", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten_5", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten_5", "inbound_nodes": [[["pool_5", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten_7", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten_7", "inbound_nodes": [[["pool_7", 0, 0, {}]]]}, {"class_name": "Concatenate", "config": {"name": "concatenate", "trainable": true, "dtype": "float32", "axis": -1}, "name": "concatenate", "inbound_nodes": [[["flatten_3", 0, 0, {}], ["flatten_5", 0, 0, {}], ["flatten_7", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_0", "trainable": true, "dtype": "float32", "units": 80, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_0", "inbound_nodes": [[["concatenate", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_d0", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_d0", "inbound_nodes": [[["dense_0", 0, 0, {}]]]}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_dGB"}, "name": "input_dGB", "inbound_nodes": []}, {"class_name": "Concatenate", "config": {"name": "concatenate_1", "trainable": true, "dtype": "float32", "axis": -1}, "name": "concatenate_1", "inbound_nodes": [[["drop_d0", 0, 0, {}], ["input_dGB", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 80, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_1", "inbound_nodes": [[["concatenate_1", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_d1", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_d1", "inbound_nodes": [[["dense_1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 60, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_2", "inbound_nodes": [[["drop_d1", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "drop_d2", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}, "name": "drop_d2", "inbound_nodes": [[["dense_2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "output", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "output", "inbound_nodes": [[["drop_d2", 0, 0, {}]]]}], "input_layers": [["input_onehot", 0, 0], ["input_dGB", 0, 0]], "output_layers": [["output", 0, 0]]}}, "training_config": {"loss": "mse", "metrics": ["mae", "mse"], "weighted_metrics": null, "loss_weights": null, "sample_weight_mode": null, "optimizer_config": {"class_name": "Adam", "config": {"name": "Adam", "learning_rate": 9.999999747378752e-05, "decay": 0.0, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "epsilon": 1e-07, "amsgrad": false}}}}
?"?
_tf_keras_input_layer?{"class_name": "InputLayer", "name": "input_onehot", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_onehot"}}
?


kernel
bias
 regularization_losses
!trainable_variables
"	variables
#	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?	
_tf_keras_layer?{"class_name": "Conv1D", "name": "conv_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "stateful": false, "config": {"name": "conv_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "filters": 100, "kernel_size": {"class_name": "__tuple__", "items": [3]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {"-1": 4}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 30, 4]}}
?


$kernel
%bias
&regularization_losses
'trainable_variables
(	variables
)	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?	
_tf_keras_layer?{"class_name": "Conv1D", "name": "conv_5", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "stateful": false, "config": {"name": "conv_5", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "filters": 70, "kernel_size": {"class_name": "__tuple__", "items": [5]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {"-1": 4}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 30, 4]}}
?


*kernel
+bias
,regularization_losses
-trainable_variables
.	variables
/	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?	
_tf_keras_layer?{"class_name": "Conv1D", "name": "conv_7", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "stateful": false, "config": {"name": "conv_7", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 30, 4]}, "dtype": "float32", "filters": 40, "kernel_size": {"class_name": "__tuple__", "items": [7]}, "strides": {"class_name": "__tuple__", "items": [1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {"-1": 4}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 30, 4]}}
?
0regularization_losses
1trainable_variables
2	variables
3	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dropout", "name": "drop_3", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "drop_3", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}
?
4regularization_losses
5trainable_variables
6	variables
7	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dropout", "name": "drop_5", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "drop_5", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}
?
8regularization_losses
9trainable_variables
:	variables
;	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dropout", "name": "drop_7", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "drop_7", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}
?
<regularization_losses
=trainable_variables
>	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "AveragePooling1D", "name": "pool_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "pool_3", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "same", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}}}
?
@regularization_losses
Atrainable_variables
B	variables
C	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "AveragePooling1D", "name": "pool_5", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "pool_5", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "same", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}}}
?
Dregularization_losses
Etrainable_variables
F	variables
G	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "AveragePooling1D", "name": "pool_7", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "pool_7", "trainable": true, "dtype": "float32", "strides": {"class_name": "__tuple__", "items": [2]}, "pool_size": {"class_name": "__tuple__", "items": [2]}, "padding": "same", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}}}
?
Hregularization_losses
Itrainable_variables
J	variables
K	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Flatten", "name": "flatten_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "flatten_3", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 1, "axes": {}}}}
?
Lregularization_losses
Mtrainable_variables
N	variables
O	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Flatten", "name": "flatten_5", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "flatten_5", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 1, "axes": {}}}}
?
Pregularization_losses
Qtrainable_variables
R	variables
S	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Flatten", "name": "flatten_7", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "flatten_7", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 1, "axes": {}}}}
?
Tregularization_losses
Utrainable_variables
V	variables
W	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Concatenate", "name": "concatenate", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "concatenate", "trainable": true, "dtype": "float32", "axis": -1}, "build_input_shape": [{"class_name": "TensorShape", "items": [null, 1400]}, {"class_name": "TensorShape", "items": [null, 910]}, {"class_name": "TensorShape", "items": [null, 480]}]}
?

Xkernel
Ybias
Zregularization_losses
[trainable_variables
\	variables
]	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dense", "name": "dense_0", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "dense_0", "trainable": true, "dtype": "float32", "units": 80, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 2790}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 2790]}}
?
^regularization_losses
_trainable_variables
`	variables
a	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dropout", "name": "drop_d0", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "drop_d0", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}
?"?
_tf_keras_input_layer?{"class_name": "InputLayer", "name": "input_dGB", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 1]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_dGB"}}
?
bregularization_losses
ctrainable_variables
d	variables
e	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Concatenate", "name": "concatenate_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "concatenate_1", "trainable": true, "dtype": "float32", "axis": -1}, "build_input_shape": [{"class_name": "TensorShape", "items": [null, 80]}, {"class_name": "TensorShape", "items": [null, 1]}]}
?

fkernel
gbias
hregularization_losses
itrainable_variables
j	variables
k	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dense", "name": "dense_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 80, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 81}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 81]}}
?
lregularization_losses
mtrainable_variables
n	variables
o	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dropout", "name": "drop_d1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "drop_d1", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}
?

pkernel
qbias
rregularization_losses
strainable_variables
t	variables
u	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dense", "name": "dense_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 60, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 80}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 80]}}
?
vregularization_losses
wtrainable_variables
x	variables
y	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dropout", "name": "drop_d2", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "drop_d2", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}
?

zkernel
{bias
|regularization_losses
}trainable_variables
~	variables
	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dense", "name": "output", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "output", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 60}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 60]}}
?
	?iter
?beta_1
?beta_2

?decay
?learning_ratem?m?$m?%m?*m?+m?Xm?Ym?fm?gm?pm?qm?zm?{m?v?v?$v?%v?*v?+v?Xv?Yv?fv?gv?pv?qv?zv?{v?"
	optimizer
 "
trackable_list_wrapper
?
0
1
$2
%3
*4
+5
X6
Y7
f8
g9
p10
q11
z12
{13"
trackable_list_wrapper
?
0
1
$2
%3
*4
+5
X6
Y7
f8
g9
p10
q11
z12
{13"
trackable_list_wrapper
?
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?layer_metrics
?layers
regularization_losses
trainable_variables
	variables
?__call__
?_default_save_signature
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
-
?serving_default"
signature_map
#:!d2conv_3/kernel
:d2conv_3/bias
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
 regularization_losses
!trainable_variables
"	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
#:!F2conv_5/kernel
:F2conv_5/bias
 "
trackable_list_wrapper
.
$0
%1"
trackable_list_wrapper
.
$0
%1"
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
&regularization_losses
'trainable_variables
(	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
#:!(2conv_7/kernel
:(2conv_7/bias
 "
trackable_list_wrapper
.
*0
+1"
trackable_list_wrapper
.
*0
+1"
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
,regularization_losses
-trainable_variables
.	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
0regularization_losses
1trainable_variables
2	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
4regularization_losses
5trainable_variables
6	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
8regularization_losses
9trainable_variables
:	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
<regularization_losses
=trainable_variables
>	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
@regularization_losses
Atrainable_variables
B	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Dregularization_losses
Etrainable_variables
F	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Hregularization_losses
Itrainable_variables
J	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Lregularization_losses
Mtrainable_variables
N	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Pregularization_losses
Qtrainable_variables
R	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Tregularization_losses
Utrainable_variables
V	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
!:	?P2dense_0/kernel
:P2dense_0/bias
 "
trackable_list_wrapper
.
X0
Y1"
trackable_list_wrapper
.
X0
Y1"
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
Zregularization_losses
[trainable_variables
\	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
^regularization_losses
_trainable_variables
`	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
bregularization_losses
ctrainable_variables
d	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 :QP2dense_1/kernel
:P2dense_1/bias
 "
trackable_list_wrapper
.
f0
g1"
trackable_list_wrapper
.
f0
g1"
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
hregularization_losses
itrainable_variables
j	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
lregularization_losses
mtrainable_variables
n	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 :P<2dense_2/kernel
:<2dense_2/bias
 "
trackable_list_wrapper
.
p0
q1"
trackable_list_wrapper
.
p0
q1"
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
rregularization_losses
strainable_variables
t	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
vregularization_losses
wtrainable_variables
x	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
:<2output/kernel
:2output/bias
 "
trackable_list_wrapper
.
z0
{1"
trackable_list_wrapper
.
z0
{1"
trackable_list_wrapper
?
?metrics
 ?layer_regularization_losses
?non_trainable_variables
?layer_metrics
?layers
|regularization_losses
}trainable_variables
~	variables
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
8
?0
?1
?2"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
?
0
1
2
3
4
5
6
7
	8

9
10
11
12
13
14
15
16
17
18
19
20
21
22"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
?

?total

?count
?	variables
?	keras_api"?
_tf_keras_metricj{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}}
?

?total

?count
?
_fn_kwargs
?	variables
?	keras_api"?
_tf_keras_metric?{"class_name": "MeanMetricWrapper", "name": "mae", "dtype": "float32", "config": {"name": "mae", "dtype": "float32", "fn": "mean_absolute_error"}}
?

?total

?count
?
_fn_kwargs
?	variables
?	keras_api"?
_tf_keras_metric?{"class_name": "MeanMetricWrapper", "name": "mse", "dtype": "float32", "config": {"name": "mse", "dtype": "float32", "fn": "mean_squared_error"}}
:  (2total
:  (2count
0
?0
?1"
trackable_list_wrapper
.
?	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
?0
?1"
trackable_list_wrapper
.
?	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
?0
?1"
trackable_list_wrapper
.
?	variables"
_generic_user_object
(:&d2Adam/conv_3/kernel/m
:d2Adam/conv_3/bias/m
(:&F2Adam/conv_5/kernel/m
:F2Adam/conv_5/bias/m
(:&(2Adam/conv_7/kernel/m
:(2Adam/conv_7/bias/m
&:$	?P2Adam/dense_0/kernel/m
:P2Adam/dense_0/bias/m
%:#QP2Adam/dense_1/kernel/m
:P2Adam/dense_1/bias/m
%:#P<2Adam/dense_2/kernel/m
:<2Adam/dense_2/bias/m
$:"<2Adam/output/kernel/m
:2Adam/output/bias/m
(:&d2Adam/conv_3/kernel/v
:d2Adam/conv_3/bias/v
(:&F2Adam/conv_5/kernel/v
:F2Adam/conv_5/bias/v
(:&(2Adam/conv_7/kernel/v
:(2Adam/conv_7/bias/v
&:$	?P2Adam/dense_0/kernel/v
:P2Adam/dense_0/bias/v
%:#QP2Adam/dense_1/kernel/v
:P2Adam/dense_1/bias/v
%:#P<2Adam/dense_2/kernel/v
:<2Adam/dense_2/bias/v
$:"<2Adam/output/kernel/v
:2Adam/output/bias/v
?2?
!__inference__wrapped_model_376956?
???
FullArgSpec
args? 
varargsjargs
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *Y?V
T?Q
*?'
input_onehot?????????
#? 
	input_dGB?????????
?2?
A__inference_model_layer_call_and_return_conditional_losses_377523
A__inference_model_layer_call_and_return_conditional_losses_377978
A__inference_model_layer_call_and_return_conditional_losses_377882
A__inference_model_layer_call_and_return_conditional_losses_377469?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
&__inference_model_layer_call_fn_378012
&__inference_model_layer_call_fn_378046
&__inference_model_layer_call_fn_377612
&__inference_model_layer_call_fn_377700?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
B__inference_conv_3_layer_call_and_return_conditional_losses_376973?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? **?'
%?"??????????????????
?2?
'__inference_conv_3_layer_call_fn_376983?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? **?'
%?"??????????????????
?2?
B__inference_conv_5_layer_call_and_return_conditional_losses_377000?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? **?'
%?"??????????????????
?2?
'__inference_conv_5_layer_call_fn_377010?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? **?'
%?"??????????????????
?2?
B__inference_conv_7_layer_call_and_return_conditional_losses_377027?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? **?'
%?"??????????????????
?2?
'__inference_conv_7_layer_call_fn_377037?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? **?'
%?"??????????????????
?2?
B__inference_drop_3_layer_call_and_return_conditional_losses_378058
B__inference_drop_3_layer_call_and_return_conditional_losses_378063?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
'__inference_drop_3_layer_call_fn_378073
'__inference_drop_3_layer_call_fn_378068?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
B__inference_drop_5_layer_call_and_return_conditional_losses_378085
B__inference_drop_5_layer_call_and_return_conditional_losses_378090?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
'__inference_drop_5_layer_call_fn_378095
'__inference_drop_5_layer_call_fn_378100?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
B__inference_drop_7_layer_call_and_return_conditional_losses_378112
B__inference_drop_7_layer_call_and_return_conditional_losses_378117?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
'__inference_drop_7_layer_call_fn_378127
'__inference_drop_7_layer_call_fn_378122?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
B__inference_pool_3_layer_call_and_return_conditional_losses_377046?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *3?0
.?+'???????????????????????????
?2?
'__inference_pool_3_layer_call_fn_377052?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *3?0
.?+'???????????????????????????
?2?
B__inference_pool_5_layer_call_and_return_conditional_losses_377061?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *3?0
.?+'???????????????????????????
?2?
'__inference_pool_5_layer_call_fn_377067?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *3?0
.?+'???????????????????????????
?2?
B__inference_pool_7_layer_call_and_return_conditional_losses_377076?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *3?0
.?+'???????????????????????????
?2?
'__inference_pool_7_layer_call_fn_377082?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *3?0
.?+'???????????????????????????
?2?
E__inference_flatten_3_layer_call_and_return_conditional_losses_378133?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
*__inference_flatten_3_layer_call_fn_378138?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
E__inference_flatten_5_layer_call_and_return_conditional_losses_378144?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
*__inference_flatten_5_layer_call_fn_378149?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
E__inference_flatten_7_layer_call_and_return_conditional_losses_378155?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
*__inference_flatten_7_layer_call_fn_378160?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
G__inference_concatenate_layer_call_and_return_conditional_losses_378168?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
,__inference_concatenate_layer_call_fn_378175?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
C__inference_dense_0_layer_call_and_return_conditional_losses_378186?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
(__inference_dense_0_layer_call_fn_378195?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
C__inference_drop_d0_layer_call_and_return_conditional_losses_378212
C__inference_drop_d0_layer_call_and_return_conditional_losses_378207?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
(__inference_drop_d0_layer_call_fn_378222
(__inference_drop_d0_layer_call_fn_378217?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
I__inference_concatenate_1_layer_call_and_return_conditional_losses_378229?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
.__inference_concatenate_1_layer_call_fn_378235?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
C__inference_dense_1_layer_call_and_return_conditional_losses_378246?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
(__inference_dense_1_layer_call_fn_378255?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
C__inference_drop_d1_layer_call_and_return_conditional_losses_378272
C__inference_drop_d1_layer_call_and_return_conditional_losses_378267?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
(__inference_drop_d1_layer_call_fn_378277
(__inference_drop_d1_layer_call_fn_378282?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
C__inference_dense_2_layer_call_and_return_conditional_losses_378293?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
(__inference_dense_2_layer_call_fn_378302?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
C__inference_drop_d2_layer_call_and_return_conditional_losses_378319
C__inference_drop_d2_layer_call_and_return_conditional_losses_378314?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
(__inference_drop_d2_layer_call_fn_378329
(__inference_drop_d2_layer_call_fn_378324?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
B__inference_output_layer_call_and_return_conditional_losses_378339?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
'__inference_output_layer_call_fn_378348?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
AB?
$__inference_signature_wrapper_377744	input_dGBinput_onehot?
!__inference__wrapped_model_376956?*+$%XYfgpqz{c?`
Y?V
T?Q
*?'
input_onehot?????????
#? 
	input_dGB?????????
? "/?,
*
output ?
output??????????
I__inference_concatenate_1_layer_call_and_return_conditional_losses_378229?Z?W
P?M
K?H
"?
inputs/0?????????P
"?
inputs/1?????????
? "%?"
?
0?????????Q
? ?
.__inference_concatenate_1_layer_call_fn_378235vZ?W
P?M
K?H
"?
inputs/0?????????P
"?
inputs/1?????????
? "??????????Q?
G__inference_concatenate_layer_call_and_return_conditional_losses_378168???~
w?t
r?o
#? 
inputs/0??????????

#? 
inputs/1??????????
#? 
inputs/2??????????
? "&?#
?
0??????????
? ?
,__inference_concatenate_layer_call_fn_378175???~
w?t
r?o
#? 
inputs/0??????????

#? 
inputs/1??????????
#? 
inputs/2??????????
? "????????????
B__inference_conv_3_layer_call_and_return_conditional_losses_376973v<?9
2?/
-?*
inputs??????????????????
? "2?/
(?%
0??????????????????d
? ?
'__inference_conv_3_layer_call_fn_376983i<?9
2?/
-?*
inputs??????????????????
? "%?"??????????????????d?
B__inference_conv_5_layer_call_and_return_conditional_losses_377000v$%<?9
2?/
-?*
inputs??????????????????
? "2?/
(?%
0??????????????????F
? ?
'__inference_conv_5_layer_call_fn_377010i$%<?9
2?/
-?*
inputs??????????????????
? "%?"??????????????????F?
B__inference_conv_7_layer_call_and_return_conditional_losses_377027v*+<?9
2?/
-?*
inputs??????????????????
? "2?/
(?%
0??????????????????(
? ?
'__inference_conv_7_layer_call_fn_377037i*+<?9
2?/
-?*
inputs??????????????????
? "%?"??????????????????(?
C__inference_dense_0_layer_call_and_return_conditional_losses_378186]XY0?-
&?#
!?
inputs??????????
? "%?"
?
0?????????P
? |
(__inference_dense_0_layer_call_fn_378195PXY0?-
&?#
!?
inputs??????????
? "??????????P?
C__inference_dense_1_layer_call_and_return_conditional_losses_378246\fg/?,
%?"
 ?
inputs?????????Q
? "%?"
?
0?????????P
? {
(__inference_dense_1_layer_call_fn_378255Ofg/?,
%?"
 ?
inputs?????????Q
? "??????????P?
C__inference_dense_2_layer_call_and_return_conditional_losses_378293\pq/?,
%?"
 ?
inputs?????????P
? "%?"
?
0?????????<
? {
(__inference_dense_2_layer_call_fn_378302Opq/?,
%?"
 ?
inputs?????????P
? "??????????<?
B__inference_drop_3_layer_call_and_return_conditional_losses_378058d7?4
-?*
$?!
inputs?????????d
p
? ")?&
?
0?????????d
? ?
B__inference_drop_3_layer_call_and_return_conditional_losses_378063d7?4
-?*
$?!
inputs?????????d
p 
? ")?&
?
0?????????d
? ?
'__inference_drop_3_layer_call_fn_378068W7?4
-?*
$?!
inputs?????????d
p
? "??????????d?
'__inference_drop_3_layer_call_fn_378073W7?4
-?*
$?!
inputs?????????d
p 
? "??????????d?
B__inference_drop_5_layer_call_and_return_conditional_losses_378085d7?4
-?*
$?!
inputs?????????F
p
? ")?&
?
0?????????F
? ?
B__inference_drop_5_layer_call_and_return_conditional_losses_378090d7?4
-?*
$?!
inputs?????????F
p 
? ")?&
?
0?????????F
? ?
'__inference_drop_5_layer_call_fn_378095W7?4
-?*
$?!
inputs?????????F
p
? "??????????F?
'__inference_drop_5_layer_call_fn_378100W7?4
-?*
$?!
inputs?????????F
p 
? "??????????F?
B__inference_drop_7_layer_call_and_return_conditional_losses_378112d7?4
-?*
$?!
inputs?????????(
p
? ")?&
?
0?????????(
? ?
B__inference_drop_7_layer_call_and_return_conditional_losses_378117d7?4
-?*
$?!
inputs?????????(
p 
? ")?&
?
0?????????(
? ?
'__inference_drop_7_layer_call_fn_378122W7?4
-?*
$?!
inputs?????????(
p
? "??????????(?
'__inference_drop_7_layer_call_fn_378127W7?4
-?*
$?!
inputs?????????(
p 
? "??????????(?
C__inference_drop_d0_layer_call_and_return_conditional_losses_378207\3?0
)?&
 ?
inputs?????????P
p
? "%?"
?
0?????????P
? ?
C__inference_drop_d0_layer_call_and_return_conditional_losses_378212\3?0
)?&
 ?
inputs?????????P
p 
? "%?"
?
0?????????P
? {
(__inference_drop_d0_layer_call_fn_378217O3?0
)?&
 ?
inputs?????????P
p
? "??????????P{
(__inference_drop_d0_layer_call_fn_378222O3?0
)?&
 ?
inputs?????????P
p 
? "??????????P?
C__inference_drop_d1_layer_call_and_return_conditional_losses_378267\3?0
)?&
 ?
inputs?????????P
p
? "%?"
?
0?????????P
? ?
C__inference_drop_d1_layer_call_and_return_conditional_losses_378272\3?0
)?&
 ?
inputs?????????P
p 
? "%?"
?
0?????????P
? {
(__inference_drop_d1_layer_call_fn_378277O3?0
)?&
 ?
inputs?????????P
p
? "??????????P{
(__inference_drop_d1_layer_call_fn_378282O3?0
)?&
 ?
inputs?????????P
p 
? "??????????P?
C__inference_drop_d2_layer_call_and_return_conditional_losses_378314\3?0
)?&
 ?
inputs?????????<
p
? "%?"
?
0?????????<
? ?
C__inference_drop_d2_layer_call_and_return_conditional_losses_378319\3?0
)?&
 ?
inputs?????????<
p 
? "%?"
?
0?????????<
? {
(__inference_drop_d2_layer_call_fn_378324O3?0
)?&
 ?
inputs?????????<
p
? "??????????<{
(__inference_drop_d2_layer_call_fn_378329O3?0
)?&
 ?
inputs?????????<
p 
? "??????????<?
E__inference_flatten_3_layer_call_and_return_conditional_losses_378133]3?0
)?&
$?!
inputs?????????d
? "&?#
?
0??????????

? ~
*__inference_flatten_3_layer_call_fn_378138P3?0
)?&
$?!
inputs?????????d
? "???????????
?
E__inference_flatten_5_layer_call_and_return_conditional_losses_378144]3?0
)?&
$?!
inputs?????????F
? "&?#
?
0??????????
? ~
*__inference_flatten_5_layer_call_fn_378149P3?0
)?&
$?!
inputs?????????F
? "????????????
E__inference_flatten_7_layer_call_and_return_conditional_losses_378155]3?0
)?&
$?!
inputs?????????(
? "&?#
?
0??????????
? ~
*__inference_flatten_7_layer_call_fn_378160P3?0
)?&
$?!
inputs?????????(
? "????????????
A__inference_model_layer_call_and_return_conditional_losses_377469?*+$%XYfgpqz{k?h
a?^
T?Q
*?'
input_onehot?????????
#? 
	input_dGB?????????
p

 
? "%?"
?
0?????????
? ?
A__inference_model_layer_call_and_return_conditional_losses_377523?*+$%XYfgpqz{k?h
a?^
T?Q
*?'
input_onehot?????????
#? 
	input_dGB?????????
p 

 
? "%?"
?
0?????????
? ?
A__inference_model_layer_call_and_return_conditional_losses_377882?*+$%XYfgpqz{f?c
\?Y
O?L
&?#
inputs/0?????????
"?
inputs/1?????????
p

 
? "%?"
?
0?????????
? ?
A__inference_model_layer_call_and_return_conditional_losses_377978?*+$%XYfgpqz{f?c
\?Y
O?L
&?#
inputs/0?????????
"?
inputs/1?????????
p 

 
? "%?"
?
0?????????
? ?
&__inference_model_layer_call_fn_377612?*+$%XYfgpqz{k?h
a?^
T?Q
*?'
input_onehot?????????
#? 
	input_dGB?????????
p

 
? "???????????
&__inference_model_layer_call_fn_377700?*+$%XYfgpqz{k?h
a?^
T?Q
*?'
input_onehot?????????
#? 
	input_dGB?????????
p 

 
? "???????????
&__inference_model_layer_call_fn_378012?*+$%XYfgpqz{f?c
\?Y
O?L
&?#
inputs/0?????????
"?
inputs/1?????????
p

 
? "???????????
&__inference_model_layer_call_fn_378046?*+$%XYfgpqz{f?c
\?Y
O?L
&?#
inputs/0?????????
"?
inputs/1?????????
p 

 
? "???????????
B__inference_output_layer_call_and_return_conditional_losses_378339\z{/?,
%?"
 ?
inputs?????????<
? "%?"
?
0?????????
? z
'__inference_output_layer_call_fn_378348Oz{/?,
%?"
 ?
inputs?????????<
? "???????????
B__inference_pool_3_layer_call_and_return_conditional_losses_377046?E?B
;?8
6?3
inputs'???????????????????????????
? ";?8
1?.
0'???????????????????????????
? ?
'__inference_pool_3_layer_call_fn_377052wE?B
;?8
6?3
inputs'???????????????????????????
? ".?+'????????????????????????????
B__inference_pool_5_layer_call_and_return_conditional_losses_377061?E?B
;?8
6?3
inputs'???????????????????????????
? ";?8
1?.
0'???????????????????????????
? ?
'__inference_pool_5_layer_call_fn_377067wE?B
;?8
6?3
inputs'???????????????????????????
? ".?+'????????????????????????????
B__inference_pool_7_layer_call_and_return_conditional_losses_377076?E?B
;?8
6?3
inputs'???????????????????????????
? ";?8
1?.
0'???????????????????????????
? ?
'__inference_pool_7_layer_call_fn_377082wE?B
;?8
6?3
inputs'???????????????????????????
? ".?+'????????????????????????????
$__inference_signature_wrapper_377744?*+$%XYfgpqz{{?x
? 
q?n
0
	input_dGB#? 
	input_dGB?????????
:
input_onehot*?'
input_onehot?????????"/?,
*
output ?
output?????????