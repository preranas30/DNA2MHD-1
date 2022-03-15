import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from config import *
from dna_diags import *
from ev_diags import *
from spectra import *
from nlt_diags import *
#from landau_tests import *
import os

if os.path.isfile('parameters.dat'):
    read_parameters()
    if os.path.isfile(par['diagdir'][1:-1]+'/parameters.dat'):
        pass
    else:
        print "Enter directory for data files:"
        diagdir=raw_input()
        print "Directory:",diagdir
        par['diagdir']="\'"+diagdir+"\'"
else:
    change_diagdir()

#get_new_files()

