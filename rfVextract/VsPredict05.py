#!/usr/bin/env python
"""
dnolivieri: (updated: 23 oct 2014)

    random forest, but using the idea of
    intervals used in the MHC exons


"""
import numpy as np
import matplotlib.pyplot as plt
import time
import os, fnmatch
import sys
import itertools
from operator import itemgetter, attrgetter
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Motif
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature


from scipy import *
import struct
import re
#from propy import PyPro
#from propy.GetProteinFromUniprot import GetProteinSequence
import json
import cPickle as pickle

rno = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

class VregionsPredict:
    def __init__(self, S, desc_method, loci_classes, mammalList):
        strand=1
        self.S = S
        self.desc_method= desc_method
        self.loci_classes = loci_classes
        self.mammalList = mammalList
        self.Lexon=320


        self.contigs = self.get_contigs(mammalList[0])
        qp = open('normalizedAA_Matrix.pkl', 'rb')
        self.normI = pickle.load(qp)
        self.predicted_seqs = []
        self.rfmodels = self.get_models()

        #self.analyze_files(mammalList)


    def get_contigs(self, mammal):
        contigs=[]
        fp = open(self.S[mammal]["contigs"], "r")
        for lk in fp:
            contigs.append(lk.strip())

        fp.close()
        print contigs
        #raw_input("The contigs")
        return contigs


    def get_models(self):
        print "Getting Models"
        rfmodels = []
        for loci in self.loci_classes:
            matfile = "trainMat_" + loci + ".pkl"
            fp = open(matfile, 'rb')
            rfmodels.append( pickle.load(fp) )

        return rfmodels


    def analyze_files(self, mammalList):
        for mammal in mammalList:
            fbar= self.S[mammal]["WGS"]
            #self.predict_seqs= self.get_VregionsRF(fbar)
            self.get_Vexon_candidates(fbar)



    def get_Vexon_candidates(self, inFile):
        print "inside get exon Vregions"
        seqbuffer=[]
        seqs_positive=[]
        seqpos_buffer=[]
        
        start_time = timeit.default_timer()
        for strand in [1,-1]:
            scnt=2
            for record in SeqIO.parse(inFile, "fasta"):
                if ( record.id.split("|")[3] not in self.contigs):
                    continue

                print record.id.split("|")[3]
                #raw_input("check if in contig")
                if strand == 1:
                    Sequence=record.seq
                else:
                    Sequence=record.seq.reverse_complement()

                ix = 0
                while ix < len(seq):
                    sbar=seq[ix: ix+20000]
                    x=[i.start()+ix for i in re.finditer("CAC", str(sbar))]
                    y=[i.start()+ix for i in re.finditer("AG", str(sbar))]
                    s=[(i,j) for i,j in itertools.product(x,y) if j>i and ( np.abs(i-j)>265  and np.abs(i-j)<285) and ((np.abs(i-j)-2)%3==0) ]

                    elapsed = timeit.default_timer() - start_time
                    print "----------", ix, "---------(", 100.*(  1.- float(len(seq) - ix)/float(len(seq))), "% )---(T=", elapsed,")---"
                    print len(s)






## ---------------MAIN ----------------------------------
if __name__ == '__main__':


    WGS_prediction = 'WGS_Prediction.json'
    json_data=open( WGS_prediction )
    S = json.load(json_data)
    json_data.close()

    classes=[ 'ighv', 'iglv', 'igkv', 'trav','trbv','trgv', 'trdv']

    mlist = []
    mlist.append(sys.argv[1])

    V = VregionsPredict(S,  desc_method='PDT',loci_classes=classes, mammalList= mlist)

