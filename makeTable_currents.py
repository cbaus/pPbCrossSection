
import numpy as np
from numpy import sqrt, pi, exp,sign

import pickle
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

inFile = open('/afs/cern.ch/user/c/cbaus/vdmanalysis/beambeamtool/fromMonika/data/CurrentsB1_Feb2013pPb_fromScan1.pkl', 'rb')
CurrentsB1 = pickle.load(inFile)
inFile.close()

inFile = open('/afs/cern.ch/user/c/cbaus/vdmanalysis/beambeamtool/fromMonika/data/CurrentsB2_Feb2013pPb_fromScan1.pkl', 'rb')
CurrentsB2 = pickle.load(inFile)
inFile.close()

assert len(CurrentsB1) == len(CurrentsB2)

print ["%0.2f" % i for i in CurrentsB2]
CurrentsB1 =  np.multiply(CurrentsB1,[1e2]*len(CurrentsB1))
CurrentsB2 =  np.multiply(CurrentsB2,[1e2]*len(CurrentsB2))

print "Mean B1: ", "%0.1f" % np.mean(CurrentsB1), "\\times10^{9}", "this must be lead?"
print "Mean B2: ", "%0.1f" % np.mean(CurrentsB2), "\\times10^{9}"
