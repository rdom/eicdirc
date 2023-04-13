import sys
import ROOT
import numpy as np

ROOT.gInterpreter.ProcessLine('#include "../src/PrtTools.h"')
ROOT.gSystem.Load('../build/libPrt.so')

t = ROOT.PrtTools("pik.root")

stat = int(sys.argv[1])
ne_train = stat

# x_train = np.zeros((ne_train,64,96))
# y_train = np.zeros((ne_train,1))

# x_test = np.zeros((10000,64,96))
# y_test = np.zeros((10000,1))

max_hits = 200
n_train = 8000

x_train = np.zeros((ne_train,max_hits,3))
y_train = np.zeros((ne_train,1))

x_test = np.zeros((n_train,max_hits,3))
y_test = np.zeros((n_train,1))



while t.next() and t.i() < ne_train + n_train :
    i = t.i()
    ind = 0
    for hit in t.event().getHits() :
        # tof = e.getTof()
        # tofPi = fEvent->getTofPi()
        # tofP = fEvent->getTofP()
        time = hit.getLeadTime()
        pmt = hit.getPmt()
        pix = hit.getPixel() - 1
        ch = int(hit.getChannel())
        
        if time > 70 :
            continue

        time_bin = int(5*time)
        tx = int(16 * (pmt % 6) + pix % 16);
        ty = int(16 * (pmt // 6) + pix / 16);

        if i < ne_train:
            x_train[i,ind, 1] = ch
            x_train[i,ind, 2] = time_bin
            y_train[i] = t.pid()
        else :
            x_test[i-ne_train,ind,1] = ch
            x_test[i-ne_train,ind,2] = time_bin
            y_test[i-ne_train] = t.pid()

        ind = ind + 1
        if(ind >= max_hits):
            break    
        

b = 0
for e in range(ne_train):
    if b>=32:
        b=0

    for i in range(max_hits):
        x_train[e,i, 0] = b
    b = b + 1

b = 0
for e in range(n_train):    
    if b>=32:
        b=0
    for i in range(max_hits):
        x_test[e,i, 0] = b
    b = b + 1
    
nid = "data_ind_" + str(ne_train);
np.savez(nid, x_train= x_train, y_train=y_train, x_test=x_test, y_test= y_test)
