#!/usr/bin/env python3
########################################################################
#    ------CRISPRspec and CRISPRoff scores computation Pipeline----
#
#  Copyright (c) 2018, 2021 by the contributors (see AUTHORS file)
#
#  This file is part of CRISPRoff
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  It is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this script, see file COPYING.
#  If not, see <http://www.gnu.org/licenses/>.
##########################################################################

## Generate CRISPRoff CRISPRspec scores for off-targeting assessment
__version__ = "1.1.2"

import sys, re, os, pickle, subprocess, argparse, gzip
from math import exp
from math import log10
from Bio import SeqIO

##############################################################
#################### Energy computation #########################
##############################################################

ENERGY_MODELS_PICKLE_FILE = "energy_dics.pkl"

RI_REV_NT_MAP = {'-':'', 'a':'T', 'A':'T', 'c':'G', 'C':'G', 'g':'C', 'G':'C',
              't':'A', 'T':'A', 'u':'A', 'U':'A', 'n':'N', 'N':'N'}

RI_DNA_DNA_NN={'AA':{'TT':-1.00}, 'TT':{'AA':-1.00}, 'AT':{'TA':-0.88}, 'TA':{'AT':-0.58},
            'CA':{'GT':-1.45}, 'TG':{'AC':-1.45}, 'GT':{'CA':-1.44}, 'AC':{'TG':-1.44},
            'CT':{'GA':-1.28}, 'AG':{'TC':-1.28}, 'GA':{'CT':-1.30}, 'TC':{'AG':-1.30},
            'CG':{'GC':-2.17}, 'GC':{'CG':-2.24}, 'GG':{'CC':-1.84}, 'CC':{'GG':-1.84}}

RI_MATCH_noGU = {'A':{'A':False, 'C':False, 'G':False, 'T':True},
         'C':{'A':False, 'C':False, 'G':True, 'T':False},
         'G':{'A':False, 'C':True, 'G':False, 'T':False},
         'T':{'A':True, 'C':False, 'G':False, 'T':False}}

########## READ THE 2mer 3mer 4mer energies ###############
RNA_DNA_internal_loop ={3:3.2, 4:3.555, 5:3.725, 6:3.975, 7:4.16, 8:4.33, 9:4.495, 10:4.6, 11:4.7}
RNA_DNA = None
def read_energy_parameters(ENERGY_MODELS_PICKLE_FILE = "energy_dics.pkl"):
    global RNA_DNA
    energy_reader=open(ENERGY_MODELS_PICKLE_FILE, 'rb')
    RNA_DNA = pickle.load(energy_reader)
    energy_reader.close()

####### Necessary for self-folding ######
RNAFOLD_EXE = "RNAfold"

######### positional energy contribution weights ###################
# LAST weight is filled but not used, (Left-over from some experimental option)
POS_WGH=[1.80067099242007, 1.95666668400006, 1.90472004401173, 2.13047270152512, 1.37853848098249, 1.46460783730408, 1.0, 1.387220146823, 1.51401000729362, 1.98058344620751, 1.87939168587699, 1.7222593588838, 2.02228445489326, 1.92692086621503, 2.08041972716723, 1.94496755678903, 2.14539112893591, 2.04277109036766, 2.24911493451185, 2.25]

# LAST weight is filled but not used, (Left-over from some experimental option)
DNA_POS_WGH = [1.22245576981774, 1.24561578622024, 1.37883177517399, 1.39146340276523, 1.24308180746857, 1.09598194424544, 1.0, 1.11695025382169, 1.11589045394936, 1.22243614188218, 1.21317477033274, 1.07125942316357, 1.25205871414019, 1.21445408158483, 1.20971491326295, 1.21076785001579, 1.2480898972246, 1.40301355270318, 1.41221084925493, 1.4]

######### pam correction parameters for pam-updated energy ##########
pam_ratios={"GGG":1.0, "AGG":1.0, "CGG":1.0, "TGG":1.0, "GAG":0.9, "AAG":0.9, "CAG":0.9, "TAG":0.9, "GGA":0.8, "AGA":0.8, "CGA":0.8, "TGA":0.8, "OTHERS": 0.0}
pam_ratio_count = 3

########## gRNA folding #######################
grna_folding_engs = {}
def get_rnafold_eng(seq, rid="temp_grna_id"):
    if seq not in grna_folding_engs:
        no_constraint_eng = None
        l = len(seq)
        cmd = RNAFOLD_EXE
        p = subprocess.Popen([cmd, '--noPS'], shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).communicate(input=">"+rid+"\n"+seq+"\n\n")
        if p[1]=="":
            no_constraint_eng = float(p[0].rstrip().split()[-1].replace('(','').replace(')',''))
        else:
            sys.stderr.write("#ERROR:RNAfold run went wrong: "+cmd+"; "+p[1]+"\n")
            exit()
        grna_folding_engs[seq] = no_constraint_eng
    return grna_folding_engs[seq]

# Stacking energy is distributed to positions only for interior loops
def calcRNADNAenergy(guideSeq, otSeq, GU_allowed=False):
    guideSeq = guideSeq.upper()[:-3]
    seq = ''.join([RI_REV_NT_MAP[c] for c in otSeq[:-3]])

    spos = -1
    epos = -1

    energy = [0.0]*len(guideSeq)

    MATCH = RI_MATCH_noGU

    for i in range(len(seq)):
        if MATCH[guideSeq[i]][seq[i]]:
            if spos==-1:
                spos = i
            epos = i

    i = spos
    while i<epos:
        j = i+1
        while MATCH[seq[j]][guideSeq[j]]==False:
            j = j+1

            if j > epos:
                break
        if j > epos:
            break

        loop_size = (j-i)-1
        eng_con = 0
        if loop_size<3:
            eng_con = RNA_DNA[loop_size][guideSeq[i:j+1]][seq[i:j+1]]
            # if there is a stack in the beginning or end AU GU penalty is still needed
            if loop_size == 0:
                if (i==spos and (guideSeq[i]=="T" or seq[i]=="T")) or (j==epos and (guideSeq[j]=="T" or seq[j]=="T")):
                    eng_con += 0.25
        else:
            eng_con = float(RNA_DNA_internal_loop[loop_size]) + float(RNA_DNA[0][guideSeq[i:i+2]][seq[i:i+2]]) + float(RNA_DNA[0][guideSeq[j-1:j+1]][seq[j-1:j+1]])

        for k in range(loop_size+1):
            energy[i+k] += eng_con/(loop_size+1)

        i = j

    return energy

################ DNA-DNA opening #############
def calcDNAopeningScore(otSeq):
    seq = otSeq.upper()[:-3]
    energy = [0.0]*len(seq)
    for i in range(1, len(seq)):
        energy[i] = float(RI_DNA_DNA_NN[seq[i-1]+seq[i]][RI_REV_NT_MAP[seq[i-1]]+RI_REV_NT_MAP[seq[i]]])
    return energy
#################################################################################################

REV_NT_MAP = {'-':'', 'a':'T', 'A':'T', 'c':'G', 'C':'G', 'g':'C', 'G':'C', 
              't':'A', 'T':'A', 'u':'A', 'U':'A', 'n':'N', 'N':'N'}

def rev_comp_seq(seq):
    s=''
    for c in seq:
        s =  REV_NT_MAP[c] + s
    return s

def comp_seq(seq):
    s=''
    for c in seq:
        s =  s + REV_NT_MAP[c]
    return s

############# Get interaction energy ##################
#Employ all the necessary computations on the score vector to get the final free energy
def get_eng(grna_seq, off_seq, score_func, GU_allowed=False, pos_weight=False, pam_corr=False, grna_folding=False, dna_opening=False, dna_pos_wgh=False):
    scores = score_func(grna_seq, off_seq, GU_allowed)

    if pos_weight:
        for i in range(len(scores)):
            if i < 21:
                scores[-(i+1)] = POS_WGH[-(i+1)] * scores[-(i+1)]

    off = sum(scores) + 0.0
    off = (-1.0) * off

    if grna_folding:
        off += get_rnafold_eng(grna_seq[:20])


    if dna_opening:
        dna_scores = calcDNAopeningScore(off_seq)
        if dna_pos_wgh:
            for i in range(len(dna_scores)):
                if i < 21:
                    dna_scores[-(i+1)] = DNA_POS_WGH[-(i+1)] * dna_scores[-(i+1)]
        off += sum(dna_scores)

    if pam_corr:
        if off_seq[-pam_ratio_count:] in pam_ratios.keys():
            off = off * pam_ratios[off_seq[-pam_ratio_count:]]
        else:
            off = off * pam_ratios["OTHERS"]

    return off

#################################################################
# TODO: Turn the ontarget and offSeq into class for ease of understanding the code #

## Compute POFF for the given grna and its off-targets ##
PAR_BETA = 1.000000 / (0.001987 * 310.150000)
def compute_CRISPRspec(ontarget, offSeqs, score_func, GU_allowed=False, pos_weight=False, pam_corr=False, grna_folding=False, dna_opening=False, dna_pos_wgh=False, ignored_chromosomes=set()):
    pf = 0.000000
    on = 0.000000
    CRISPRoff_scores = []

    for offSeq in offSeqs:
        offSeq_eng = get_eng(ontarget[0], offSeq[0], score_func,  GU_allowed=GU_allowed, pos_weight=pos_weight, pam_corr=pam_corr, grna_folding=grna_folding, dna_opening=dna_opening, dna_pos_wgh=dna_pos_wgh)
        if offSeq[1] not in ignored_chromosomes:
            pf = pf + exp(PAR_BETA * offSeq_eng)
        else:
            sys.stderr.write("#WARNING: This off-target sequence is ignored when computing CRISPRspec: '"+"|".join([str(x) for x in offSeq])+"'.\n")
        CRISPRoff_scores.append((offSeq, offSeq_eng))

    on_eng = get_eng(ontarget[0], ontarget[0], score_func, GU_allowed=GU_allowed, pos_weight=pos_weight, pam_corr=pam_corr, grna_folding=grna_folding, dna_opening=dna_opening, dna_pos_wgh=dna_pos_wgh)
    on = exp(PAR_BETA * on_eng)

    CRISPRoff_scores.append(((ontarget[0], ontarget[2], ontarget[3], ontarget[4], ontarget[5]), on_eng))
    return (pf/(pf+on)), CRISPRoff_scores

###############################################################
############## On-target efficiency ###########################
## Compute the on-target prediction score with Azimuth
# Input is 30nts, 4nt before the target site, 3nt after the PAM
def get_ontarget_scores_30nt(ontargets_30):
    import azimuth.model_comparison
    import numpy as np
    sequences = np.array(ontargets_30)
    try:
        predictions = azimuth.model_comparison.predict(sequences, None, None)
    except:
        sys.stderr.write("#WARNING: Elevation error: Could not compute scores for "+",".join(ontargets_30)+"\n")
        return dict(zip(ontargets_30, [0.0]*len(ontargets_30)))
    return dict(zip(ontargets_30, list(predictions)))


###############################################################
################ Input readers ################################
###############################################################
# Returns all valid Cas9 targets predicted with RIsearch2 in the genome
# Result file must have been generated with -p3 option
# grna is the sequence without the PAM addition.
# If PAM is added last parameter must be set to False
def read_risearch_results(guideSeq, ris_file, noPAM_given=True, count_mms=False, on_targets = [], offSeqs = [], off_counts = {"GG":[0]*7, "AG":[0]*7, "GA":[0]*7}):
    inf = gzip.open(ris_file, 'rt') if ris_file.endswith(".gz") else open(ris_file, 'rt')
    x = 0
    for line in inf:
        x = x + 1
        cols = line.rstrip().split("\t")
        if len(cols)>11 and len(cols[10])>3 and len(cols[11])>0:
            gid, qs, qe, tc, ts, te, tst, en, ist, iseq, pamseq, preseq= cols[:12]
            PAM = comp_seq(pamseq[:3])
            if PAM[1:3] in ["GG", "AG", "GA"]:
                offseq = comp_seq(iseq)
                if tst=="+":
                    ts = str(int(ts)-4)
                    tst="-"
                else:
                    ts = str(int(ts)-1)
                    te = str(int(te)+3)
                    tst="+"

                # determine and save the number of mismatches for this guide
                mm_count = sum ( offseq[i] != guideSeq[i] for i in range(len(offseq)) )
                if count_mms:
                    if mm_count<7:
                        off_counts[PAM[1:3]][mm_count] += 1
                if mm_count<7:
                    if (noPAM_given and (guideSeq != offseq)) or ((noPAM_given==False) and (guideSeq != offseq+PAM)):
                        offSeqs.append(((offseq+PAM), tc, ts, te, tst))
                    else:
                        on_targets.append(((offseq+PAM), (rev_comp_seq(preseq[:4])+offseq+comp_seq(pamseq[:6])), tc, ts, te, tst))
        else:
            sys.stderr.write("WARNING: RIsearch output line is missing columns, or info in columns.\nLINE: "+line)
    inf.close()
    return  offSeqs, off_counts, on_targets

# Returns all valid Cas9 targets from Cas-OFFinder result file
def read_casoff_results(guideSeq, casoff_file, count_mms=False):
    on_targets = []
    offSeqs = []
    off_counts = {"GG":[0]*7, "AG":[0]*7, "GA":[0]*7}
    inf = gzip.open(casoff_file, 'rt') if casoff_file.endswith(".gz") else open(casoff_file, 'rt')
    x = 0
    for line in inf:
        x = x + 1
        cols = line.rstrip().split("\t")
        if line[0]!="#" and len(cols)>7 and cols[0]=="X":
            bulge_type, crRNA, DNA, chromosome, position, strand, mismatches, bulge_size = cols[:8]
            PAM = DNA.upper()[-3:]
            if PAM[1:3] in ["GG", "AG", "GA"]:
                offseq = DNA.upper()
                mm_count = int(mismatches)
                # determine and save the number of mismatches for this guide
                if count_mms:
                    if mm_count<7:
                        off_counts[PAM[1:3]][mm_count] += 1

                if mm_count<7:
                    if DNA[:len(guideSeq)]==guideSeq:
                        if strand=="+":
                            on_targets.append(((offseq), DNA, chromosome, position, str(int(position)+len(offseq)), strand))
                        else:
                            on_targets.append(((offseq), DNA, chromosome, str(int(position)-len(offseq)), position, strand))
                    else:
                        if strand=="+":
                            offSeqs.append(((offseq), chromosome, position, str(int(position)+len(offseq)), strand))
                        else:
                            offSeqs.append(((offseq), chromosome, str(int(position)-len(offseq)), position, strand))
    inf.close()
    return  offSeqs, off_counts, on_targets

# Returns all valid Cas9 targets from plain new-line seperated input off-targets file
def read_standard_offtargets_input(guideSeq, in_file, count_mms=False):
    on_targets = []
    offSeqs = []
    off_counts = {"GG":[0]*7, "AG":[0]*7, "GA":[0]*7}
    inf = gzip.open(in_file, 'rt') if in_file.endswith(".gz") else open(in_file, 'rt')
    x = 0
    for line in inf:
        cols = line.rstrip().split("\t")
        if line[0]!="#" and len(cols[0])==len(guideSeq):
            x = x + 1
            offseq = cols[0].upper()
            PAM = offseq[-3:]
            if PAM[1:3] in ["GG", "AG", "GA"]:
                mm_count = sum ( offseq[i] != guideSeq[i] for i in range(len(offseq)-3) )
                # determine and save the number of mismatches for this guide
                if count_mms:
                    if mm_count<7:
                        off_counts[PAM[1:3]][mm_count] += 1
                if mm_count<7:
                    if offseq==guideSeq:
                        on_targets.append((offseq, offseq, None, None, None, None))
                    else:
                        offSeqs.append((offseq, None, None, None, None))
    inf.close()
    return  offSeqs, off_counts, on_targets

#Read given off-target file whether its risearch casoff or standard regular new-line seperated file
def read_offtargets_file(guideSeq, offtargets_file, noPAM_given=False, count_mms=False, chromosome_names=None):
    if "risearch" in offtargets_file:
        if chromosome_names==None:
            sys.stdout.write('#RUNNING: reading given risearch output file "'+offtargets_file+'" for gRNA/on-target:'+guideSeq+'.\n')
            return read_risearch_results(guideSeq, offtargets_file, noPAM_given=noPAM_given, count_mms=count_mms, on_targets = [], offSeqs = [], off_counts = {"GG":[0]*7,"AG":[0]*7,"GA":[0]*7})
        else: # READ from multiple files
            chr_inf = gzip.open(chromosome_names, 'rt') if chromosome_names.endswith(".gz") else open(chromosome_names, 'rt')
            on_targets = []
            offSeqs = []
            off_counts = {"GG":[0]*7, "AG":[0]*7, "GA":[0]*7}
            ris_dir = "/".join(offtargets_file.split("/")[:-1])
            ris_file = offtargets_file.split("/")[-1]
            for chrid in chr_inf:
                sys.stdout.write('#RUNNING: reading given risearch output file "'+"/".join([ris_dir, chrid.rstrip(), ris_file])+'".\n')
                if os.path.isfile("/".join([ris_dir, chrid.rstrip(), ris_file])):
                    offSeqs, off_counts, on_targets = read_risearch_results(guideSeq, "/".join([ris_dir, chrid.rstrip(), ris_file]), noPAM_given=noPAM_given, count_mms=count_mms, on_targets = on_targets, offSeqs = offSeqs, off_counts = off_counts)
            chr_inf.close()
            return offSeqs, off_counts, on_targets
    else:
        inf = gzip.open(offtargets_file, 'rt') if offtargets_file.endswith(".gz") else open(offtargets_file, 'rt')
        casoff = True if inf.readline().startswith("#Bulge") else False
        inf.close()
        if casoff:
            sys.stdout.write('#RUNNING: reading given Cas-OFFinder file "'+offtargets_file+'".\n')
            return read_casoff_results(guideSeq, offtargets_file, count_mms=count_mms)
        else:
            sys.stdout.write('#RUNNING: reading given offtarget file "'+offtargets_file+'" as regular off-target file.\n')
            return read_standard_offtargets_input(guideSeq, offtargets_file, count_mms=count_mms)

#check if file is fasta
def is_fasta(filename):
    with open(filename, "rt") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

# Read guide sequences from fasta file
def read_guides_fasta(fasta_file):
    ontargets={}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        if len(str(record.seq))==23:
            ontargets[record.id] = seq
        else:
            # get valid guides with PAM NGG and 23nt
            x=0
            grna_seq=""
            for i in range(len(seq)):
                if seq[i] not in ["A", "C", "G", "T", "U"]:
                    grna_seq =""
                elif len(grna_seq)<22:
                    grna_seq = grna_seq + seq[i]
                else:
                    grna_seq = grna_seq + seq[i]
                    if grna_seq[-2:] == "GG":
                        x+=1
                        ontargets[record.id+"_gRNA_"+str(x)] = grna_seq
                    grna_seq = grna_seq[1:]
            revseq = rev_comp_seq(seq)
            for i in range(len(revseq)):
                if revseq[i] not in ["A", "C", "G", "T", "U"]:
                    grna_seq =""
                elif len(grna_seq)<22:
                    grna_seq = grna_seq + revseq[i]
                else:
                    grna_seq = grna_seq + revseq[i]
                    if grna_seq[-2:] == "GG":
                        x+=1
                        ontargets[record.id+"_rev_gRNA_"+str(x)] = grna_seq
                    grna_seq = grna_seq[1:]

            if x==1:
                ontargets[record.id] = ontargets[record.id+"_gRNA_"+str(x)]
                del ontargets[record.id+"_gRNA_"+str(x)]
    return ontargets

# Read guide sequences from fasta file
def get_guides(in_file, notfile=False):
    #Read the sequences
    seq=""
    if notfile:
        seq = in_file
    else:
        with gzip.open(in_file, 'rt') if in_file.endswith(".gz") else open(in_file, 'rt') as inf:
            for line in inf:
                seq = seq + line.rstrip().upper()

    if len(seq)!=23:
        seq = seq+"N"+rev_comp_seq(seq)

    ontargets={}
    # get valid guides with PAM NGG and 23nt
    x=0
    grna_seq=""
    for i in range(len(seq)):
        if seq[i] not in ["A", "C", "G", "T", "U"]:
            grna_seq = ""
        elif len(grna_seq)<22:
            grna_seq = grna_seq + seq[i]
        else:
            grna_seq = grna_seq + seq[i]
            if grna_seq[-2:] == "GG":
                x+=1
                ontargets["gRNA_"+str(x)] = grna_seq
            grna_seq = grna_seq[1:]
    return ontargets

## Summarize the features for given on-targets
def get_energy_features_for_guides(guideSeqs):
    feature_out = {}
    for guide, guideSeq in guideSeqs.items():
        feature_out[guide]={}
        feature_out[guide]["RNA_DNA_eng"] = -1.00000 * get_eng(guideSeq, guideSeq, calcRNADNAenergy, GU_allowed=False, pos_weight=False, pam_corr=False, grna_folding=False, dna_opening=False, dna_pos_wgh=False)
        feature_out[guide]["RNA_DNA_eng_weighted"] = -1.00000 * get_eng(guideSeq, guideSeq, calcRNADNAenergy, GU_allowed=False, pos_weight=True, pam_corr=False, grna_folding=False, dna_opening=False, dna_pos_wgh=False)
        feature_out[guide]["DNA_DNA_opening"] = -1.00000 * sum(calcDNAopeningScore(guideSeq))
        feature_out[guide]["spacer_self_fold"] = get_rnafold_eng(guideSeq[:20])
        feature_out[guide]["CRISPRoff_score"] = get_eng(guideSeq, guideSeq, calcRNADNAenergy, GU_allowed=False, pos_weight=True, pam_corr=True, grna_folding=True, dna_opening=True, dna_pos_wgh=False)
    features=["RNA_DNA_eng", "RNA_DNA_eng_weighted", "DNA_DNA_opening", "spacer_self_fold", "CRISPRoff_score"]
    return features, feature_out

## Summarize the on off target info for given on_targets


# setup the argument parser
def get_parser():
    parser = argparse.ArgumentParser(description="CRISPR-OFF Webserver - Computational Cas9 off-targeting assessment / CRISPRspec and CRISPRoff scores computation pipeline v"+__version__)
    grna_input = parser.add_argument_group("# gRNA/on-target sequence(s)")
    grna_input.add_argument("--guides", metavar="<file>", type=str,
                        help="Fasta file for gRNA sequences (each gRNA is 23nt=20nt+PAM) or sequence file (fasta or unformatted) of the target to design gRNAs for.")
    grna_input.add_argument("--guide", metavar="<seq>", type=str,
                        help="single gRNA sequence to analyse (23nt, 20nt+PAM)")

    key_settings = parser.add_argument_group('# Key Settings')
    key_settings.add_argument("--duplex_energy_params", metavar='<pickle file>', type=str,
                        help="Pickled energy parameters file for nucleic acid duplexes (default: energy_dics.pkl)",
                        default='energy_dics.pkl')
    key_settings.add_argument("--rnafold_x", metavar='<command>', type=str,
                        help="RNAfold executable (default: RNAfold)",
                        default='RNAfold')

    specificity_analysis = parser.add_argument_group('# Analyze the specificity of given on-target sequences')
    specificity_analysis.add_argument("--offtargets", metavar="<off-target file>", type=str, default=None,
                         help="Targets file to read off-targets from. It could be RIsearch2 -p3 or Cas-OFFinder result files, or plain offtargets file where off-target sequences are line-seperated. Note that on-target sequence always needs to be part of the given target data.")
    specificity_analysis.add_argument("--risearch_results_folder", metavar="<path>", type=str, default=None,
                         help="Path to the folder with risearch -p3 output files (Ignores the file given by --offtargets and looks for risearch results in the given folder)")
    specificity_analysis.add_argument("--chromosome_names", metavar="<file>", type=str, default=None,
                         help="File with the name of the chromosomes (new-line seperated) to be able to locate chromosome-predicted risearch result files. When passed, risearch result files have to be located within chromosome-named sub-directories under the directory given by --risearch_results_folder")
    specificity_analysis.add_argument("--ignored_chromosomes", metavar="<file>", type=str, default=None,
                         help="Ignore chromosomes in the given list (new-line seperated) when computing the CRISPRspec specificity score. For example, certain haplotypes can be ignored for more accurate specificity computation.")
    specificity_analysis.add_argument("--no_azimuth", action="store_true",
                         help="Don't report azimuth on-target efficiency prediction scores. (Pass this argument if azimuth python package is not installed in your system.)")
    specificity_analysis.add_argument("--no_off_target_counts", action="store_true",
                         help="Don't report off-target counts")
    specificity_analysis.add_argument("--sorted_CRISPRoff_reports", action="store_true",
                         help="Sort off-targets in the report based on computed CRISPRoff score. (default: Unsorted report)")
    specificity_analysis.add_argument("--report_top", metavar="<int>", type=int, default=None,
                         help="Report only top X off-targets.")
    specificity_analysis.add_argument("--evaluate_all", action="store_true",
                         help="Pass this to evaluate all gRNAs even though their on-target is not part of the given target predictions. (default: Ignore gRNAs with no on-target)")
    specificity_analysis.add_argument("--comment_out_NAs", action="store_true",
                         help="Comment out the lines with no genomic coordinate info in CRISPRoff result files. Required for further bed-intersect within the webserver.")

    result_output = parser.add_argument_group("# Output options")
    result_output.add_argument("--specificity_report", metavar='<file>', type=str,
                        help="Output file for the on/off-target assesment of given on-target sequences. You can pass stdout or sterr for printing (default:stdout).",
                        default="stdout")
    result_output.add_argument("--CRISPRoff_scores_folder", metavar="<path>", type=str, default=None,
                        help="Path to the folder where CRISPRoff score reports are stored for each guide. (default: NOT generated.)")

    parameter_generator = parser.add_argument_group('# Compute the energy parameters for given on-target sequences')
    parameter_generator.add_argument("--guide_params_out", metavar='<file>', type=str,
                        help="Output file for computed energy parameters of given on-target sequences. You can pass 'stdout' or 'stderr' as well to print within the terminal.",
                        default=None)

    return parser

# MAIN FUNCTION here
def main():
    sys.stdout.write('#START: PIPELINE RUN HAS STARTED.\n')

    # STEP 1: Get the necessary arguments
    parser = get_parser()
    args = parser.parse_args()
    sys.stdout.write('#STEP 1: arguments parsed\n')
    for k in sorted(args.__dict__):
        if (args.__dict__[k] is not None) and args.__dict__[k] != '':
            sys.stdout.write('#ARG: ' + k + ' = ' + str(args.__dict__[k]) + '\n')

    # STEP 2: Save the energy parameters and the path to RNAfold executable
    read_energy_parameters(args.duplex_energy_params)
    sys.stdout.write('#STEP 2: energy parameters parsed from "'+args.duplex_energy_params+'"\n')
    global RNAFOLD_EXE
    RNAFOLD_EXE = args.rnafold_x

    # STEP 3: READ the guide sequences
    guideSeqs = {}
    if args.guides == None:
        if args.guide!=None:
            guideSeqs = get_guides(args.guide.upper(), notfile=True)
            sys.stdout.write('#STEP 3: "'+args.guide+'" has been read as gRNA/target seq input and we generated '+str(len(guideSeqs.keys()))+' valid gRNAs from this sequence.\n')
    elif os.path.isfile(args.guides):
        if is_fasta(args.guides):
            guideSeqs = read_guides_fasta(args.guides)
        else:
            guideSeqs = get_guides(args.guides)
        sys.stdout.write('#STEP 3: '+str(len(guideSeqs.keys()))+' guides read/generated from "'+args.guides+'"\n')
    else:
        sys.stderr.write('#ERROR: "'+args.guides+'" is not a readable input file.\n')
        exit()

    if len(guideSeqs.keys())==0:
        sys.stdout.write('#ABORTED: No guide sequence given\n')
        sys.stderr.write('#ERROR: There are no guide sequences given as input.\n')
        exit()

    # STEP 4: WRITE on-target energy features (parameters)
    if args.guide_params_out != None:
        outf = sys.stdout
        if args.guide_params_out=="stderr":
            outf = sys.stderr
        elif args.guide_params_out!="stdout":
            outf = gzip.open(args.guide_params_out, "wt") if args.guide_params_out.endswith(".gz") else open(args.guide_params_out, "wt")

        feature_names, feature_values_dic = get_energy_features_for_guides(guideSeqs)
        outf.write("\t".join(["guideID", "guideSeq"]+feature_names)+"\n")
        for guide, feature_dic in feature_values_dic.items():
            outf.write("\t".join([guide, guideSeqs[guide]]+[str(feature_dic[feature]) for feature in feature_names])+"\n")

        if args.guide_params_out not in ["stderr", "stdout"]:
            outf.close()

        sys.stdout.write('#STEP 4: on-target features have been recorded in "'+args.guide_params_out+'"\n')
    else:
        sys.stderr.write('#WARNING: Skipping STEP 4 (NO energy parameters report for the guides). \n')
        sys.stdout.write('#STEP 4: Skipping. No output chosen for reporting the energy features of the guides.\n')

    # STEP 5: off-target analysis of gRNAs
    if args.risearch_results_folder != None or args.offtargets!=None:
        ignored_chromosomes=set()
        if args.ignored_chromosomes!=None:
            if os.path.isfile(args.ignored_chromosomes):
                with open(args.ignored_chromosomes, 'rt') as inf:
                    for line in inf:
                        ignored_chromosomes.add(line.rstrip())
            else:
                sys.stdout.write('#STEP 5: Given --ignored_chromosomes file is not readable.\n')
                exit()

        #Specificity analysis results, including CRISPRoff and CRISPRspec results
        outf = sys.stdout
        if args.specificity_report=="stderr":
            outf = sys.stderr
        elif args.specificity_report!="stdout":
            outf = gzip.open(args.specificity_report, "wt") if args.specificity_report.endswith(".gz") else open(args.specificity_report, "wt")

        # Specificity report output
        specificity_report = '\t'.join(["Guide_ID", "Guide_sequence", "On_target_30nt", "Genomic_position", "CRISPRspec_specificity_score", "Azimuth_ontarget_score", "MM_counts", "MM_detailed\n"])
        # iterate over all guides
        for guideID, guideSeq in guideSeqs.items():
            offSeqs, off_counts, ontargets = None, None, None

            if args.guides!=None and args.guide!=None and guideSeq!=args.guide:
                continue
            offtargets_file = os.path.join(args.risearch_results_folder, "risearch_"+guideID+".out.gz") if args.risearch_results_folder!=None else args.offtargets
            sys.stdout.write('#STEP 5.1: off-targets are read from "'+(offtargets_file if args.risearch_results_folder==None else args.risearch_results_folder)+'".\n')
            if os.path.isfile(offtargets_file) or (args.risearch_results_folder!=None and args.chromosome_names!=None):
                offSeqs, off_counts, ontargets = read_offtargets_file(guideSeq, offtargets_file, noPAM_given=(len(guideSeq)<23), count_mms=(args.no_off_target_counts==False), chromosome_names=args.chromosome_names)

                # If there are no on-targets for this guide add it manually to target space or ignore based on --evaluate_all
                if len(ontargets)==0:
                    if args.evaluate_all:
                        sys.stderr.write("WARNING: on-target couldn't be found for "+guideID+":"+guideSeq+". Evaluating this guide by adding the guide sequence to the target space.\n")
                        ontargets.append((guideSeq, guideSeq, None, None, None, None))
                    else:
                        sys.stderr.write("WARNING: on-target couldn't be found for "+guideID+":"+guideSeq+". Skipping this guide.\n")
                        # Write CRISPRoff results
                        if args.CRISPRoff_scores_folder != None:
                            if os.path.isdir(args.CRISPRoff_scores_folder):
                                with open(os.path.join(args.CRISPRoff_scores_folder, guideSeq+".CRISPRoff.tsv"), "wt") as score_outf:
                                    score_outf.write('# No CRISPRoff and CRISPRspec score can be computed for "'+guideID+","+guideSeq+'"\n')
                                    score_outf.write('# REASON: On-target sequence do not exist within"'+offtargets_file+'".\n')
                        # Specificity report output
                        specificity_report += '\t'.join([guideID, guideSeq, guideSeq, "NA", "0", "NA", "NA", "NA"+"\n"])
                        continue

                azimuth_score_dic = {} if args.no_azimuth else get_ontarget_scores_30nt([ontarget[1] for ontarget in ontargets if len(ontarget[1])==30])
                for i in range(len(ontargets)):
                    on_target, on_target_30, tc, ts, te, tst = ontargets[i]
                    if tc in ignored_chromosomes:
                        continue

                    if guideSeq[:20]!=on_target[:20]:
                        sys.stderr.write('#WEIRD ERROR: "'+guideSeq[:20]+'" should be the same as "'+on_target[:20]+'"\n')

                    # Compute azimuth scores
                    azimuth_on = str(azimuth_score_dic[on_target_30]) if on_target_30 in azimuth_score_dic else "NA"

                    # Prepare off-target set
                    offs = [offSeq for offSeq in offSeqs]
                    for j in [k for k in range(len(ontargets)) if k!=i]:
                        off_on_target, off_on_target_30, off_tc, off_ts, off_te, off_tst = ontargets[j]
                        offs.append((off_on_target, off_tc, off_ts, off_te, off_tst))

                    #Compute CRISPRoff and CRISPRspec
                    on_prob, CRISPRoff_scored_offs = compute_CRISPRspec(ontargets[i], offs, calcRNADNAenergy, GU_allowed=False, pos_weight=True, pam_corr=True, grna_folding=True, dna_opening=True, dna_pos_wgh=False, ignored_chromosomes=ignored_chromosomes)
                    CRISPRspec = -1.0 * log10(on_prob) if on_prob>0 else "+Infinite"

                    # Write CRISPRoff results
                    if args.CRISPRoff_scores_folder != None:
                        if os.path.isdir(args.CRISPRoff_scores_folder):
                            with open(os.path.join(args.CRISPRoff_scores_folder, on_target+".CRISPRoff.tsv"), "wt") as score_outf:
                                score_outf.write('# CRISPRoff scored off-targets of "')
                                score_outf.write(",".join([on_target]+['NA' if v is None else v for v in [tc, ts, te, tst]])+'"')
                                score_outf.write(' (Off-targets read from "'+offtargets_file+'").\n')
                                score_outf.write("# Off-targets are given in bed-like format\n")
                                score_outf.write("\t".join(["# chromosome", "start", "end", "off_target_seq", "CRISPRoff_score", "strand"])+"\n")

                                report_count=0
                                for offtargets, CRISPRoff in sorted(CRISPRoff_scored_offs, key=lambda x:x[1], reverse=True) if args.sorted_CRISPRoff_reports else CRISPRoff_scored_offs:
                                    report_count += 1
                                    if args.comment_out_NAs and offtargets[1]==None:
                                        score_outf.write("#"+"\t".join(['NA' if v is None else v for v in list(offtargets[1:4])+[offtargets[0], str(CRISPRoff), offtargets[4]]])+"\n")
                                    else:
                                        score_outf.write("\t".join(['NA' if v is None else v for v in list(offtargets[1:4])+[offtargets[0], str(CRISPRoff), offtargets[4]]])+"\n")
                                    if args.report_top!=None and report_count==args.report_top:
                                        break

                    # Prepare text for mismatch count info
                    MM_counts = ",".join([ str(sum([ l[j] for l in off_counts.values()])) for j in range(7)])
                    MM_detailed = ";".join(["N"+PAM+":"+",".join([ str(off_counts[PAM][j]) for j in range(7)]) for PAM in ["GG", "AG", "GA"]])

                    # Specificity report output
                    specificity_report += '\t'.join([guideID, guideSeq, on_target_30, "|".join(['NA' if v is None else v for v in [tc, ts, te, tst]]), str(CRISPRspec), azimuth_on, MM_counts, MM_detailed+"\n"])

                    sys.stdout.write('#STEP 5.2: Scores computed and reported for given gRNA and off-targets.\n')

            else:
                sys.stderr.write('#WARNING: Skipping "'+guideID+':'+guideSeq+'" for off-target assessment (No off-targets file or risearch result folder given).\n')

        # Write the specificity report
        sys.stdout.write('#STEP 5.3: Reporting the specificity scores')
        if args.specificity_report in ["stdout", "stderr"]:
            outf.write(" below.\n\n"+specificity_report+"\n")
        else:
            sys.stdout.write(" and saving into '"+args.specificity_report+"' file.\n")
            outf.write(specificity_report)
            outf.close()
        sys.stdout.write('#STEP 5: Reporting DONE.\n')

    else:
        sys.stderr.write('#WARNING: Skipping STEP 5 (NO off-target assessment report for the guides). \n')
        sys.stdout.write('#STEP 5: Skipping. No output file specificied and/or No folder given for risearch results.\n')

    # FINISH
    sys.stdout.write('#END: FINISHED RUNNING THE PIPELINE WITH NO ERROR.\n')


if __name__== "__main__":
    main()
