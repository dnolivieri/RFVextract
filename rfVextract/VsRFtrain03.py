#!/usr/bin/env python
"""
  dnolivieri: (started: 23 september 2014)
     - doing multiple binary trainings.

"""
import pylab as pl
import numpy as np
import sys

import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.cross_validation import cross_val_score
from numpy import genfromtxt, savetxt
import time
import itertools
import cPickle as pickle
import timeit
import pybedtools as pb

from Bio import SeqIO
from Bio.Seq import Seq

from scipy import *
import struct
import re
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence
import featureHeatmap01 as fHM
import json


AA = {1:'A', 2:'R',3:'N',4:'D',5:'C', 6:'Q', 7:'E', 8:'G', 9:'H', 10:'I',11:'L',  12:'K', 13:'M', 14:'F', 15:'P', 16:'S', 17:'T', 18:'W', 19:'Y',  20:'V',  21:'B', 22:'Z', 23:'X'}

rno = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

n_classes = 2
n_estimators = 1000
RANDOM_SEED = 13


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

class TrainingPipeline:
    def __init__(self, S, loci_classes, desc_method, bckgnd_method, pos_exists=False, neg_exists=False, do_training=False):
        self.normI=self.normalized_AAindx()


        self.S = S
        self.loci_classes  = loci_classes
        self.desc_method   = desc_method
        self.bckgnd_method = bckgnd_method
        self.pos_exists= pos_exists
        self.neg_exists= neg_exists

        self.posFile    = 'train_posSignal'
        self.bckgndFile = 'train_bckgndSignal'

        self.plot_featurevecs = False

        #self.form_pos_featurevecs()
        #self.form_neg_featurevecs()


        train_signals, train_vals =  self.form_featurevecs()
        #train_signals, train_vals =  self.form_featurevecs_fromfiles()


        npzsignals_out = "train_signals.npz"
        np.savez(npzsignals_out, dp=train_signals )

        npzvals_out = "train_vals.npz"
        np.savez(npzvals_out, dp=train_vals )


        do_training=True
        if do_training:
            A=TrainClassifier( train_signals, train_vals, 1000 )



    def form_pos_featurevecs(self):
        zP=self.get_positive_signals()

    def form_neg_featurevecs(self):
        zB=self.get_background()

    def form_featurevecs_fromfiles(self):
        dataP = np.load( self.posFile + '.npz')
        zP = dataP['dp']
        dataB = np.load( self.bckgndFile + '.npz')
        zB = dataB['dp']


        xMat=[]
        xMat.append(zP)
        xMat.append(zB)
        sigLabels = ['P','B']

        valP = np.ones( [zP.shape[0],1] )
        valB = np.zeros( [zB.shape[0],1] )
        train_signals = np.vstack([zP, zB])
        train_vals = np.vstack([valP, valB])

        return train_signals, train_vals




    def add_positive_signals(self, fname):
        D = self.descriptors_from_fasta(fname)
        return D


    def form_featurevecs(self):

        if self.neg_exists:
            zB=self.get_existing_neg_background()
        else:
            zB=self.get_background()

        npzoutfile=self.bckgndFile+".npz"
        np.savez(npzoutfile, dp=np.array(zB) )

        valB = np.zeros( [zB.shape[0],1] )

        print "-----Done with Background -----"



        if self.pos_exists:
            zP=self.get_existing_pos_signals()
        else:
            zP = np.array([])
            valP = np.array([])
            lcnt=1
            for loci in self.loci_classes:
                sig_fname = self.S['All'][loci]
                zprm=self.add_positive_signals(sig_fname)
                valprm = lcnt * np.ones( [zprm.shape[0],1] )
                if lcnt==1:
                    zP = zprm
                    valP = valprm
                else:
                    zP = np.vstack([zP, zprm])
                    valP = np.vstack([valP, valprm])
                lcnt+=1

        npzoutfile=self.posFile+".npz"
        np.savez(npzoutfile, dp=np.array(zP) )



        if self.plot_featurevecs:
            xMat=[]
            xMat.append(zP)
            xMat.append(zB)
            sigLabels = ['P','B']
            #xbarMat=np.log(xMat)
            #H = fHM.heatMap(xbarMat, sigLabels)
            H = fHM.heatMap(xMat, sigLabels)


        train_signals = np.vstack([zP, zB])
        train_vals = np.vstack([valP, valB])

        return train_signals, train_vals



    def get_existing_pos_signals(self):
        data_posSignal = self.posFile+'.npz'
        dataP = np.load( data_posSignal)
        zP = dataP['dp']

        return zP

    def get_existing_neg_background(self):
        """
        data_bkgSignal1 = './train_randbackgrnd.npz'
        data_bkgSignal2 = './train_extrctBkg.npz'
        dataB1 = np.load( data_bkgSignal1)
        dataB2 = np.load( data_bkgSignal2)
        zB1 = dataB1['dp']
        zB2 = dataB2['dp']
        zB = np.vstack([zB1, zB2])
        """
        data_bkgSignal = self.bckgndFile + '.npz'
        dataB = np.load( data_bkgSignal)
        zB = dataB['dp']
        return zB


    def descriptors_from_fasta(self, infile):
        # other methods available:
        # 'GetAAComp', 'GetAAindex1', 'GetAAindex23', 'GetALL', 'GetAPAAC', 'GetCTD',
        # 'GetDPComp', 'GetGearyAuto', 'GetGearyAutop', 'GetMoranAuto', 'GetMoranAutop',
        # 'GetMoreauBrotoAuto', 'GetMoreauBrotoAutop', 'GetPAAC', 'GetPAACp', 'GetQSO',
        # 'GetQSOp', 'GetSOCN', 'GetSOCNp', 'GetSubSeq', 'GetTPComp',

        qbar=[]
        cnt=0
        for record in SeqIO.parse(infile, "fasta"):
            descObject=PyPro.GetProDes(record.seq.tostring())
            if ('X' not in record.seq.tostring()) and ('Z' not in record.seq.tostring()) and ('B' not in record.seq.tostring()):
                if self.desc_method=='AAComp':
                    T = descObject.GetAAComp()
                elif self.desc_method=='GestQSO':
                    T = descObject.GetQSO()
                elif self.desc_method=='GetGearyAuto':
                    T = descObject.GetMoranAuto()
                elif self.desc_method=='GetCTD':
                    T=descObject.GetCTD()
                elif self.desc_method=='GetPAAC':
                    T =descObject.GetPAAC()
                elif self.desc_method=='PDT':
                    T = self.getPDT3(record.seq.tostring())
                else:
                    T=descObject.GetCTD()
                Tx = [ T[x]  for x in T.iterkeys() ]
                #print Tx
                #raw_input('press to continue')
                print cnt, Tx[0:5]
                qbar.append(Tx)
                cnt+=1
                if cnt>1e9:
                    break


        return np.array(qbar)


    def get_positive_signals(self):
        infile = './dataVgeneDB/All_mammals_Vs.fasta'
        D = self.descriptors_from_fasta(infile)
        npzoutfile=self.posFile+".npz"
        np.savez(npzoutfile, dp=np.array(D) )
        return D


    def get_background(self):
        if self.bckgnd_method=='Mine':
            infile='./dataVgeneDB/bkg.fasta'
            D = self.descriptors_from_fasta(infile)
            #npzoutfile=self.bckgndFile+".npz"
            #np.savez(npzoutfile, dp=np.array(D) )
        if self.bckgnd_method=='Mine_with_mut':
            D= self.fasta_with_mutations(infile)
            npzoutfile=self.bckgndFile+".npz"
            np.savez(npzoutfile, dp=np.array(D) )
        elif self.bckgnd_method=='Total_random':
            D=self.total_random_backgrnd()
            npzoutfile=self.bckgndFile+".npz"
            np.savez(npzoutfile, dp=np.array(D) )


        return D

    def fasta_with_mutations(self, infile):
        for record in SeqIO.parse(infile, "fasta"):
            sbar = record.seq.tostring()
            sbarL= list(sbar)
            posSeq = list(np.random.randint(1,len(sbar), int(0.95*len(sbar))) )
            print posSeq
            for p in posSeq:
                aindx=np.random.randint(1,21,1)[0]
                sbarL[p] = AA[aindx]

            rbar = ''.join(sbarL)

            descObject=PyPro.GetProDes(rbar)
            print record.seq
            print rbar
            T =descObject.GetPAAC()
            Tx = [ T[x]  for x in T.iterkeys() ]
            print Tx
            qbar.append(Tx)

        return np.array(qbar)


    def total_random_backgrnd(self):
        qbar=[]
        for i in range(50):
            indx = list(np.random.randint(1,21, 90))
            seqList = [ AA[j] for j in indx ]
            #print seqList
            seq =''
            rbar = seq.join(seqList)
            print rbar
            descObject=PyPro.GetProDes(rbar)
            if self.desc_method=='AAComp':
                T = descObject.GetAAComp()
            elif self.desc_method=='GestQSO':
                T = descObject.GetQSO()
            elif self.desc_method=='GetGearyAuto':
                T = descObject.GetMoreauBrotoAuto()
            elif self.desc_method=='GetCTD':
                T=descObject.GetCTD()
            elif self.desc_method=='GetPAAC':
                T =descObject.GetPAAC()
            elif self.desc_method=='PDT':
                T = self.getPDT3(rbar)
            else:
                T=descObject.GetCTD()
            Tx = [ T[x]  for x in T.iterkeys() ]
            print i, Tx[0:5]
            #raw_input('press to continue')
            qbar.append(Tx)

        return np.array(qbar)


    def get_AAdescriptors(self):
        pass

    def normalized_AAindx(self):
        fp = open('aaindex.txt','r')
        D=[]
        for lk in fp:
            q=[ float(i) for i in lk.split() ]
            D.append(q)

        Dvec=[]
        normI = []
        for j in  D:
            q= np.sum(np.array( j ))/20.
            denom=0.0
            for kp in j:
                denom+= (kp - q)*(kp - q)

            denom = np.sqrt(denom/20.)
            abar=[]
            for kp in j:
                abar.append( (kp - q)/denom )

            normI.append(abar)

        save_object(normI, r'normalizedAA_Matrix.pkl')
        return normI

    def getPDT3(self, seq):
        Dvec={}
        cnt=0
        for q in self.normI:
            sumDseq=0.0
            for i in range(len(seq)-3):
                sumDseq+= (q[rno[seq[i]]] - q[rno[seq[i+3]]])*(q[rno[seq[i]]] - q[rno[seq[i+3]]])

            sumDseq = sumDseq/np.float(len(seq)-3)
            Dvec.update( {str(cnt): sumDseq} )
            cnt+=1
        return Dvec


# ------------------
class TrainEachLoci:
    def __init__(self, S, loci_classes, desc_method, Nestimators):
        self.S = S
        self.desc_method   = desc_method
        self.Nestimators = Nestimators


        self.make_training_matrices()


    def get_existing_signals(self, infile):
        data = np.load( infile )
        z = data['dp']
        return z


    def make_training_matrices(self):
        # background
        zB1 = self.get_existing_signals(self.S['Bkgnd']['bkg1'] )
        print "zB1.shape=", zB1.shape
        zB2 = self.get_existing_signals(self.S['Bkgnd']['bkg2'] )
        print "zB1.shape=", zB2.shape

        zB = np.vstack([zB1, zB2])
        valB = np.zeros( [zB.shape[0],1] )

        for loci in self.S['Loci'].iterkeys():
             zP = self.get_existing_signals(self.S['Loci'][loci] )
             valP = np.ones( [zP.shape[0],1] )
             print loci,  "zP.shape=", zP.shape
             train_signals = np.vstack([zP, zB])
             train_vals = np.vstack([valP, valB])
             outfile = "trainMat_"+ loci + ".pkl"

             A=TrainClassifier( train_signals, train_vals, self.Nestimators, outfile)


# ------------------
class TrainClassifier:
    def __init__(self, train_signals, train_vals, Nestimators, action=None):
        self.train_signals = train_signals
        #self.train_vals    = train_vals
        self.Nestimators  = Nestimators
        self.outfile = "test.pkl"

        #self.train_vals=  train_vals.reshape( [train_vals.shape[0], 1])
        self.train_vals=  np.ravel(train_vals)

        do_training=False
        if do_training:
            print "do training"
            start_time = timeit.default_timer()
            rf=self.do_training()
            #save_object(rf, r'train_Matrix.pkl')
            save_object(rf, self.outfile)
            elapsed = timeit.default_timer() - start_time
            print "ELAPSED=", elapsed

        count_vectors=True
        if count_vectors:
            self.obtain_trainig_set_size()

    def do_training(self):
        rf = RandomForestClassifier(n_estimators=self.Nestimators, oob_score=True)
        rf.fit(self.train_signals, self.train_vals)

        return rf


    def obtain_trainig_set_size(self):
        print np.count_nonzero(train_vals)
        sbar=set(np.ravel(train_vals).astype(int))

        for k in list(sbar):
            print k, np.extract( train_vals==k, train_vals).size



## ---------------MAIN ----------------------------------
if __name__ == '__main__':


    Vs_Loci = 'Vs_Loci_npz.json'

    json_data=open( Vs_Loci )
    S = json.load(json_data)
    json_data.close()


    method="PDT"
    loci_classes=[ 'ighv', 'iglv', 'igkv', 'trav','trbv','trgv', 'trdv']


    """
    if method=="PDT":
        T = TrainingPipeline( S, loci_classes, desc_method='PDT', bckgnd_method='Mine', pos_exists=False, neg_exists=True, do_training=False)
    elif method=="PAAC":
        T = TrainingPipeline( S, loci_classes, desc_method='GetPAAC', bckgnd_method='Mine', pos_exists=False, neg_exists=False, do_training=True)

    """


    
    #if method=="PDT":
    #    T= TrainEachLoci(S, loci_classes, desc_method='PDT',  Nestimators=50)



    do_training=True
    if do_training:

        npzsingals_in = "train_signals.npz"
        npzvals_in = "train_vals.npz"
        signals = np.load(npzsingals_in)
        vals = np.load( npzvals_in )

        train_signals = signals['dp']
        train_vals = vals['dp']

        n_classifiers=5000
        print  "n_classifiers=", n_classifiers
        A=TrainClassifier( train_signals, train_vals, 5000 )

