#/usr/bin/env python
# File: dna_diags.py

import re
import math
import numpy as np
import matplotlib.pyplot as plt
from config import *
import scipy as sc
import scipy.special
import matplotlib
from pylab import *
import os
import scipy.optimize
import scipy.linalg as lin
from scipy import interpolate
#from vtktools import *

def get_new_files():
    """Checks if links exist to dna output files. \n
    If links exists then these are deleted and new  \n
    links are created."""
    if os.path.isfile('energy.dat'):
        os.system('rm energy.dat')
    if os.path.isfile('energy3d.dat'):
        os.system('rm energy3d.dat')
    if os.path.isfile('fmom3d.dat'):
        os.system('rm fmom3d.dat')
    if os.path.isfile('ffm.dat'):
        os.system('rm ffm.dat')
    if os.path.isfile('nlt_shells.dat'):
        os.system('rm nlt_shells.dat')
    if os.path.isfile('energy_hermite.dat'):
        os.system('rm energy_hermite.dat')
    if os.path.isfile('eshells.dat'):
        os.system('rm eshells.dat')
    if os.path.isfile('shell_info.dat'):
        os.system('rm shell_info.dat')
    if os.path.isfile('g_out.dat'):
        os.system('rm g_out.dat')
    
    if os.path.isfile('nlt_shells_n0.dat'):
        os.system('rm nlt_shells_n0.dat')
    if os.path.isfile('nlt_shells_n1.dat'):
        os.system('rm nlt_shells_n1.dat')
    if os.path.isfile('nlt_shells_n2.dat'):
        os.system('rm nlt_shells_n2.dat')
    if os.path.isfile('nlt_shells_n3.dat'):
        os.system('rm nlt_shells_n3.dat')
    if os.path.isfile('nlt_shells_n1o8.dat'):
        os.system('rm nlt_shells_n1o8.dat')
    if os.path.isfile('nlt_shells_n2o8.dat'):
        os.system('rm nlt_shells_n2o8.dat')
    if os.path.isfile('nlt_shells_n3o8.dat'):
        os.system('rm nlt_shells_n3o8.dat')
    if os.path.isfile('nlt_shells_n4o8.dat'):
        os.system('rm nlt_shells_n4o8.dat')
    if os.path.isfile('nlt_shells_n5o8.dat'):
        os.system('rm nlt_shells_n5o8.dat')
    if os.path.isfile('nlt_shells_n6o8.dat'):
        os.system('rm nlt_shells_n6o8.dat')
    if os.path.isfile('nlt_shells_n7o8.dat'):
        os.system('rm nlt_shells_n7o8.dat')
    
    print "diagdir:",par['diagdir']
    #os.system('ls '+par['diagdir'])
    #if os.path.isfile(par['diagdir'][1:-1]+'/energy.dat'):
    #    print "We have a file!"
    #else:
    #    print "We don't have a file." 
    #    print par['diagdir'][1:-1]+'/energy.dat'
    
    if os.path.isfile(par['diagdir'][1:-1]+'/energy.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/energy.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/energy3d.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/energy3d.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/fmom3d.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/fmom3d.dat')
    else:
        print "fmom3d.dat doesn't exist!"
    if os.path.isfile(par['diagdir'][1:-1]+'/ffm.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/ffm.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/energy_hermite.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/energy_hermite.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/eshells.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/eshells.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/shell_info.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/shell_info.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/g_out.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/g_out.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells.dat')
    
    
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n0.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n0.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n1.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n1.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n2.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n2.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n3.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n3.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n1o8.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n1o8.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n2o8.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n2o8.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n3o8.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n3o8.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n4o8.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n4o8.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n5o8.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n5o8.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n6o8.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n6o8.dat')
    if os.path.isfile(par['diagdir'][1:-1]+'/nlt_shells_n7o8.dat'):
        os.system('ln -s '+par['diagdir'][1:-1]+'/nlt_shells_n7o8.dat')

def change_diagdir(diagdir='none'):
    """Changes diagdir.  The new diagdir can be entered as a keyword. \n
    If diagdir='keep' is selected, then the diagdir corresponding to \n
    the existing parameters.dat file is used. \n
    Otherwise you will be prompted for a new diagdir."""


    if(diagdir=='none'):
        print "Enter directory:"
        diagdir=raw_input()
        print "Directory:",diagdir
        if(diagdir != 'keep'):
            if os.path.isfile('parameters.dat'):
                os.system('rm parameters.dat')
        #os.system('ls -l '+diagdir)

    if(diagdir != 'keep'):
        if os.path.isfile(diagdir+'/parameters.dat'):
            os.system('ln -s '+diagdir+'/parameters.dat')
            #read_parameters()
            read_parameters()
            if(par['diagdir'][1:-1]) != diagdir:
                par['diagdir']="\'"+diagdir+"\'"
            get_new_files()
        else:
            print "Not a file:"
            stop
    else:
        read_parameters()
        if os.path.isfile(par['diagdir'][1:-1]+'/parameters.dat'):
            pass
        else:
            print "Enter directory for data files:"
            diagdir=raw_input()
            print "Directory:",diagdir
            par['diagdir']="\'"+diagdir+"\'"

def read_parameters():
    """Reads parameters from parameters.dat \n
    The parameters are in a dictionary call par \n
    and can be accessed via par['parameter_name']"""
    parfile=open('./parameters.dat','r')
    parameters_in=parfile.read()
    lines=parameters_in.split('\n')
    #    parameters={}
    #note: par comes from config.py
    num_lines=len(lines)
    print "Number of lines", num_lines
    print lines[0]
    for i in range(num_lines):
         temp=lines[i].split()    
         if temp:
              str_check_namelist=re.match("&",temp[0])
         if str_check_namelist:
              current_namelist=temp[0]
              print current_namelist
              namelists[current_namelist]=" "
         if len(temp)>2:
              #if (re.match(\d):
              str_check_sn=re.match("\d*\.?\d*[eE]-?\+?\d*",temp[2])
              str_check_int=re.match("\d*",temp[2])
              str_check_float=re.match("\d*\.\d*",temp[2])
              if (str_check_sn and str_check_sn.end()==len(temp[2])):
                   par[temp[0]]=float(temp[2])  
                   namelists[current_namelist]=namelists[current_namelist]+" "+temp[0]
              elif (str_check_float and str_check_float.end()==len(temp[2])):
                   par[temp[0]]=float(temp[2])
                   namelists[current_namelist]=namelists[current_namelist]+" "+temp[0]
              elif (str_check_int and str_check_int.end()==len(temp[2])):
                   float_temp=float(temp[2])
                   par[temp[0]]=int(float_temp)
                   namelists[current_namelist]=namelists[current_namelist]+" "+temp[0]
              else:
                   par[temp[0]]=temp[2]
                   namelists[current_namelist]=namelists[current_namelist]+" "+temp[0]

    #par['kxmax']=(par['nkx0']-1)*par['kxmin']
    #par['kymax']=(par['nky0']/2-1)*par['kymin']
    #par['kzmax']=(par['nkz0']/2-1)*par['kzmin']
    par['ky_nyq']=(par['nky0']/2)*par['kymin']
    par['kz_nyq']=(par['nkz0']/2)*par['kzmin']
    if par['etg_factor'] != 0.0:
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "Warning! field solver in dna diags not implement for ky=0 and etg_factor != 0."
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!"

def get_time_from_gout(swap_endian=False,gnl = False):
   """Returns time array taken from g_out.dat"""
   if gnl:
      file_name = 'g_out.dat'
   else:
      file_name = 'gnl_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*par['nv0']
   mem_tot=ntot*16
   time=np.empty(0)
   continue_read=1
   i=0
   while (continue_read): 
     f.seek(i*(mem_tot+8))
     i=i+1
     input=np.fromfile(f,dtype='float64',count=1)
     if swap_endian:
         input=input.newbyteorder()
     #print input
     if input==0 or input:
         time = np.append(time,input)
     else:
         continue_read=0

   f.close()
   return time

def get_time_from_energy():
   """Returns time array taken from energy3d.dat"""
   f = open('energy3d.dat','rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']
   mem_tot=4*ntot*8
   time=np.empty(0)
   continue_read=1
   i=0
   while (continue_read): 
     f.seek(i*(mem_tot+8))
     i=i+1
     input=np.fromfile(f,dtype='float64',count=1)
     #print input
     if input==0.0 or input:
         time = np.append(time,input)
     else:
         continue_read=0

   f.close()
   return time

def get_time_from_fmom3d():
   """Returns time array taken from fmom3d.dat"""
   f = open('fmom3d.dat','rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']
   mem_tot=2*ntot*16
   time=np.empty(0)
   continue_read=1
   i=0
   while (continue_read): 
     f.seek(i*(mem_tot+8))
     i=i+1
     input=np.fromfile(f,dtype='float64',count=1)
     #print input
     if input==0.0 or input:
         time = np.append(time,input)
     else:
         continue_read=0
   f.close()
   return time

def read_phi_fmom3d(which_itime):
   """Returns time time step from fmom3d.dat"""
   f = open('fmom3d.dat','rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']
   mem_tot=2*ntot*16
   phi=np.empty((par['nkx0'],par['nky0'],par['nkz0']),dtype='complex128')
   #press=np.empty((par['nkx0'],par['nky0'],par['nkz0']),dtype='complex128',count=ntot)
   f.seek(8+which_itime*(mem_tot+8))
   phi=np.fromfile(f,dtype='complex128',count=ntot)
   #f.seek(8+which_itime*(2*mem_tot+8)+mem_tot)
   #press=np.fromfile(f,dtype='float64',count=1)

   return phi

def write_parameter_file(par_out,path,suffix):
        """Writes a parameters file from current parameters. \n
        par_out:  an existing parameters dictionary \n
        path:  output path \n
        suffix: output file is parameters+suffix"""
        f=open(path+'/parameters'+suffix,'w')
        #f.write('string \n')
        for nml in namelists:
                f.write(nml+'\n')
                for var in par_out:
                        if re.search(' '+var+' ',namelists[nml]):
                                f.write(var+" = "+str(par_out[var])+'\n')
                f.write('/\n\n')
        f.close()

def plot_ffm():
    """Plots heat flux and phi squared data from ffm.dat"""
    dat=np.genfromtxt('ffm.dat') # takes a data file and generates an array with the same dimension
    plt.figure(figsize=(8.0,5.0))
    fig=plt.gcf()
    fig.subplots_adjust(bottom=0.22)
    fig.subplots_adjust(left=0.18)
    plt.plot(dat[:,0],dat[:,2])
    plt.ylabel(r'$Q(v_{ti}P_0\rho^2/R^2)$',size=13)
    plt.xlabel(r'$t(R/v_{ti})$',size=13)
    plt.show()
    plt.figure(figsize=(8.0,5.0))
    fig=plt.gcf()
    fig.subplots_adjust(bottom=0.22)
    fig.subplots_adjust(left=0.22)
    plt.plot(dat[:,0],dat[:,1])
    plt.ylabel(r'$\phi^2(e^2/T_i^2)$',size=13)
    plt.xlabel(r'$t(R/v_{ti})$',size=13)
    plt.show()

def plot_energy(start_time=-1.0,end_time=-1.0):
    """Plots data from energy.dat"""
    dat=np.genfromtxt('energy.dat') # takes a data file and generates an array with the same dimension
    plt.figure(figsize=(8.0,5.0))
    fig=plt.gcf()
    fig.subplots_adjust(bottom=0.22)
    fig.subplots_adjust(left=0.18)
    plt.plot(dat[:,0],dat[:,4],label=r'$Q$')
    plt.plot(dat[:,0],dat[:,5],label=r'$C$')
    plt.plot(dat[:,0],dat[:,6],label=r'$hyp-C$')
    plt.plot(dat[:,0],dat[:,7],label=r'$hyp-xyz$')
    plt.plot(dat[:,0],dat[:,9],label=r'$hyp-conv$')
    plt.legend(loc='upper left')
    plt.xlabel(r'$t(R/v_{ti})$',size=13)
    plt.show()
    plt.title('balance test')
    plt.plot(dat[:,0],dat[:,3],label='RHS')
    plt.plot(dat[:,0],dat[:,10],label='dE/dt')
    plt.plot(dat[:,0],np.sum(dat[:,4:10],axis=1),'x',label='sum')
    plt.legend()
    plt.xlabel(r'$t(R/v_{ti})$',size=13)
    plt.show()

    time=dat[:,0]
    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    Qavg=np.sum(dat[istart:iend,4])/ntime
    Cavg=np.sum(dat[istart:iend,5])/ntime
    HCavg=np.sum(dat[istart:iend,6])/ntime
    Hperpavg=np.sum(dat[istart:iend,7])/ntime
    Hconvavg=np.sum(dat[istart:iend,9])/ntime
    print "Q average:",Qavg
    print "C average:",Cavg
    print "Hyp coll average:",HCavg
    print "Hyp perp average:",Hperpavg
    print "Hyp conv average:",Hconvavg



def get_energy_single_k(g_1,g_2,kx_in,ky_in,kz_in,which_energy):
    """Computes the selected energy-related quantity (summed over Hermite) for the wavenumber \n
    defined by kx_in, ky_in, kz_in.  \n
    \'which_energy\' determines the energy quantity returned: \n
    which_energy=-1:  Free energy 
    which_energy=-2:  Electrostatic part 
    which_energy=-3:  Entropy part
    which_energy=0:  Total RHS of energy equation
    which_energy=1:  Collisions
    which_energy=2:  Hyper collisions
    which_energy=3:  Phase mixing
    which_energy=4:  omn term
    which_energy=5:  kz phi term
    which_energy=6:  omt (drive) term
    which_energy=7:  omt FLR term
    which_energy=8:  hyps \n
    which_energy=9:  nonlinear term \n
    This takes in two vectors so that a scalar product can be taken between to vectors."""
    #Total energy
    #-1==energy
    if(which_energy==-1):
        eop=energy_operator_single_k(g_1,kx_in,ky_in) # computes the term in the brackets on the rhs in Eq. (53)
        rhs=0.5*g_2
    #-2==electrostatic energy
    elif(which_energy==-2):
        eop=energy_operator_single_k(g_1,kx_in,ky_in,which_part=2)
        rhs=0.5*g_2
        #print "es energy",np.real(np.sum(eop*rhs))
    #-3==entropy
    elif(which_energy==-3):
        eop=energy_operator_single_k(g_1,kx_in,ky_in,which_part=1)
        rhs=0.5*g_2
    #0==nonlinear term
    elif(which_energy==9):
        #Note, must pass entire distribution function!!
        ikx=get_kindex(kx_in,par['kxmin'],par['nkx0'])
        iky=get_kindex(ky_in,par['kymin'],par['nky0'])
        ikz=get_kindex(kz_in,par['kzmin'],par['nkz0'])
        eop=energy_operator_single_k(g_1[ikx,iky,ikz,:],kx_in,ky_in)
        gshape=np.shape(g_2)
        if len(gshape) != 4:
            print "Error in get_energy_single_k!"
        rhs=get_rhs_nl_single_k(g_2,kx_in,ky_in,kz_in)
    #Some energy term from RHS
    else:
        eop=energy_operator_single_k(g_1,kx_in,ky_in)
        #0==all
        #1==Collisions
        #2==hyper collisions
        #3==phase mixing
        #4==omn term
        #5==kz phi term
        #6==omt term
        #7==omt FLR term
        #8==hyps
        rhs=get_rhs_lin_single_k(g_2,which_energy,kx_in,ky_in,kz_in)
    return np.real(np.sum(eop*rhs))

def get_henergy_single_k(g_1,g_2,kx_in,ky_in,kz_in,which_energy):
    """Computes the selected energy-related quantity for the wavenumber \n
    defined by kx_in, ky_in, kz_in.  \n
    \'which_energy\' determines the energy quantity returned: \n
    which_energy=-1:  Free energy 
    which_energy=-2:  Electrostatic part 
    which_energy=-3:  Entropy part
    which_energy=0:  Total linear RHS of energy equation
    which_energy=1:  Collisions
    which_energy=2:  Hyper collisions
    which_energy=3:  Phase mixing
    which_energy=4:  omn term
    which_energy=5:  kz phi term
    which_energy=6:  omt (drive) term
    which_energy=7:  omt FLR term
    which_energy=8:  hyps \n
    which_energy=9:  nonlinear term \n
    This takes in two vectors so that a scalar product can be taken between to vectors."""
    #Total energy
    #-1==energy
    if(which_energy==-1):
        eop=energy_operator_single_k(g_1,kx_in,ky_in) # computes the term in the brackets on the rhs in Eq. (53)
        rhs=0.5*g_2
    #-2==electrostatic energy
    elif(which_energy==-2):
        eop=energy_operator_single_k(g_1,kx_in,ky_in,which_part=2)
        rhs=0.5*g_2
        #print "es energy",np.real(np.sum(eop*rhs))
    #-3==entropy
    elif(which_energy==-3):
        eop=energy_operator_single_k(g_1,kx_in,ky_in,which_part=1)
        rhs=0.5*g_2
    #0==nonlinear term
    elif(which_energy==9):
        #Note, must have the entire distribution function!!
        ikx=get_kindex(kx_in,par['kxmin'],par['nkx0'])
        iky=get_kindex(ky_in,par['kymin'],par['nky0'])
        ikz=get_kindex(kz_in,par['kzmin'],par['nkz0'])
        eop=energy_operator_single_k(g_1[ikx,iky,ikz,:],kx_in,ky_in)
        gshape=np.shape(g_2)
        if len(gshape) != 4:
            print "Error in get_energy_single_k!"
        rhs=get_rhs_nl_single_k(g_2,kx_in,ky_in,kz_in)
    #Some energy term from RHS
    else:
        eop=energy_operator_single_k(g_1,kx_in,ky_in)
        #0==all
        #1==Collisions
        #2==hyper collisions
        #3==phase mixing
        #4==omn term
        #5==kz phi term
        #6==omt term
        #7==omt FLR term
        #8==hyps
        rhs=get_rhs_lin_single_k(g_2,which_energy,kx_in,ky_in,kz_in)
    return np.real(eop*rhs)

def energy_operator_single_k(g_in,kx_in,ky_in,which_part=0):
    """Returns the energy operator (see slab_kinetic.pdf) derived from g_in
    which_part=0: both terms
    which_part=1: entropy part
    which_part=3: phi part"""

    #which_part==0==>all
    #which_part==1==>g part
    #which_part==2==>phi part
    phi=get_phi_single_k(g_in,kx_in,ky_in)
    if which_part==0 or which_part==1:
        E_op=np.pi**0.5*np.conj(g_in)

    if which_part==0:
        E_op[0]=E_op[0]+np.pi**0.25*gyro_avg_k(kx_in,ky_in)*np.conj(phi)

    if which_part==2:
        E_op=np.zeros(len(g_in),dtype='complex128')
        E_op[0]=np.pi**0.25*gyro_avg_k(kx_in,ky_in)*np.conj(phi)
        #print "eop sum",np.sum(E_op)

    return E_op

def get_grids():
    """Returns kx,ky,kz,Hermite grids in the same form as used in the code \n
    kxgrid = 0, kxmin, . . . kxmax \n
    kygrid = 0, kymin, . . . kymax, kymax+kymin, -kymax, . . . -kymin"""
    kxgrid=np.arange((par['nkx0']))
    kxgrid=kxgrid*par['kxmin']
    kygrid=np.empty(par['nky0'])
    kzgrid=np.empty(par['nkz0'])
    herm_grid=np.arange(par['nv0'])
    herm_grid=1.0*herm_grid
    for i in range(par['nky0']/2):
        kygrid[par['nky0']-1-i]=-float(i+1)*par['kymin']
        kygrid[i]=float(i)*par['kymin']
    kygrid[par['nky0']/2]=par['nky0']/2*par['kymin']
    for i in range(par['nkz0']/2):
        kzgrid[par['nkz0']-1-i]=-float(i+1)*par['kzmin']
        kzgrid[i]=float(i)*par['kzmin']
    kzgrid[par['nkz0']/2]=par['nkz0']/2*par['kzmin']
    return kxgrid,kygrid,kzgrid,herm_grid

def read_time_step_g(which_itime,swap_endian=False,gnl = False):
   """Reads a time step from g_out.dat.  Time step determined by \'which_itime\'"""
   if gnl:
      file_name = 'gnl_out.dat'
   else:
      file_name = 'g_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*par['nv0']
   mem_tot=ntot*16
   gt0=np.empty((par['nv0'],par['nkz0'],par['nky0'],par['nkx0']))
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='complex128',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)
   f.close()
   return gt0

def read_time_step_gn(which_itime,n_in,swap_endian=False):
   """Reads a time step for a certain n from g_out.dat.  Time step determined by \'which_itime\'"""
   f = open('g_out.dat','rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']
   mem_tot=ntot*16
   #gt0=np.empty((par['nv0'],par['nkz0'],par['nky0'],par['nkx0']))
   f.seek(8+which_itime*(8+mem_tot*par['nv0'])+ntot*n_in*16)
   gt0=np.fromfile(f,dtype='complex128',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)
   f.close()
   return gt0

def get_time_from_gkfile(ikx,iky,read_nl=False):
    cikx=str(int(ikx))
    if len(cikx)==1:
        cikx='0'+cikx
    elif len(cikx) > 2:
        print cikx
        print "Error in get_time_from_gk"
        stop
    
    ciky=str(int(iky))
    if len(ciky)==1:
        ciky='0'+ciky
    elif len(ciky) > 2:
        print "Error in get_time_from_gk"
        stop

    diagdir=par['diagdir'][1:-1]
    print diagdir
        
    if read_nl:
        file_name='nl_kx'+cikx+'ky'+ciky+'.dat'
    else:
        file_name='g_kx'+cikx+'ky'+ciky+'.dat'

    print "Reading file",diagdir+'/'+file_name
    file_exists=os.path.isfile(diagdir+'/'+file_name)
    if file_exists:
        pass
    else:
        print "File does not exist:",diagdir+'/'+file_name
        stop

    f=open(diagdir+'/'+file_name,'r')
    ntot=par['nv0']*par['nkz0']
    mem_tot=ntot*16
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0

    f.close()
    return time



def read_time_step_gkfile(ikx,iky,which_itime,read_nl=False):
    cikx=str(int(ikx))
    if len(cikx)==1:
        cikx='0'+cikx
    elif len(cikx) > 2:
        print "Error in get_time_from_gk"
        stop
    
    ciky=str(int(iky))
    if len(ciky)==1:
        ciky='0'+ciky
    elif len(ciky) > 2:
        print "Error in get_time_from_gk"
        stop
    
    #cikz=str(int(ikz))
    #if len(cikz)==1:
    #    cikz='0'+cikz
    #elif len(cikz) > 2:
    #    print "Error in get_time_from_gk"
    #    stop

    diagdir=par['diagdir'][1:-1]
    #print diagdir
        
    if read_nl:
        file_name='nl_kx'+cikx+'ky'+ciky+'.dat'
    else:
        file_name='g_kx'+cikx+'ky'+ciky+'.dat'
    #print "Reading file",diagdir+'/'+file_name
    file_exists=os.path.isfile(diagdir+'/'+file_name)
    if file_exists:
        pass
    else:
        print "File does not exist:",diagdir+'/'+file_name
        stop


    f = open(diagdir+'/'+file_name,'rb')
    ntot=par['nv0']*par['nkz0']
    mem_tot=ntot*16
    #gt0=np.empty((par['nkz0'],par['nv0']))
    f.seek(8+which_itime*(8+mem_tot))
    gt0=np.fromfile(f,dtype='complex128',count=ntot)
    #print sum(gt0)
    f.close()
    return gt0

#def read_time_step_gk(kx,ky,kz,which_itime):
#    """Reads a time step from a g_k*** file.  The file is determined by kx,ky,kz, \n
#    and the time step is determined by \'which_itime\' """
#    ikx=int(kx/par['kxmin'])
#    iky=int(ky/par['kymin'])
#    if(iky<0):
#        iky=iky+par['nky0']
#    ikz=int(kz/par['kzmin'])
#    if(ikz<0):
#        ikz=ikz+par['nkz0']
#
#    #print "kx,ikx",kx,ikx
#    #print "ky,iky",ky,iky
#    #print "kz,ikz",kz,ikz
#
#    f = open('g_out.dat','rb')
#    ntot=par['nkx0']*par['nky0']*par['nkz0']*par['nv0']
#    mem_tot=ntot*16
#    kshift=(ikz*par['nky0']*par['nkx0']+iky*par['nkx0']+ikx)
#    #print 'kshift',kshift
#    gk0=np.zeros(par['nv0'],dtype='complex')
#    for i in range(par['nv0']):
#      f.seek(8+which_itime*(8+mem_tot)+kshift*16+16*i*par['nkx0']*par['nky0']*par['nkz0'])
#      a=np.fromfile(f,dtype='complex',count=1)
#      gk0[i]=a[0]
#      #print gk0[i]
#    #print sum(gt0)
#    f.close()
#    return gk0

#def test_read_time_step(kx,ky,kz,which_itime):
#    
#    ikx=int(kx/par['kxmin'])
#    iky=int(ky/par['kymin'])
#    if(iky<0):
#        iky=iky+par['nky0']
#    ikz=int(kz/par['kzmin'])
#    if(ikz<0):
#        ikz=ikz+par['nkz0']
#
#    gt0=read_time_step_g(which_itime)
#    gt1=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
#    gk0=gt1[ikx,iky,ikz,:]
#
#    gk=read_time_step_gk(kx,ky,kz,which_itime)
#
#    plt.plot(abs(gk),'-x',label='From read gk')
#    plt.plot(abs(gk0),'-+',label='From read g')
#    plt.legend()
#    plt.show()
#    
#    plt.plot(np.real(gk),'-x',label='Re From read gk')
#    plt.plot(np.real(gk0),'-+',label='Re From read g')
#    plt.plot(np.imag(gk),'-x',label='Im From read gk')
#    plt.plot(np.imag(gk0),'-+',label='Im From read g')
#    plt.legend()
#    plt.show()

def ksum_3dF(obj_in):
    """This routine returns the sum over all kx,ky,kz for a 3D object. \n
    It accounts for the implicitly deined negative-kx modes by multiplying \n
    all kx>0 modes by 2.0"""
    ksum=2.0*np.sum(np.sum(np.sum(obj_in[1:,:,:],axis=0),axis=0),axis=0)
    ksum=ksum+np.sum(np.sum(obj_in[0,:,:],axis=0),axis=0)
    #ksum=np.sum(obj_in[1:par['nkx0']-1,0:par['nky0']-2,0:par['nkz0']-2])
    #ksum=ksum*2.0
    #ksum=ksum+np.sum(obj_in[0,0:par['nky0']-2,0:par['nkz0']-2])
    return ksum 

def ksum_3d_n(obj_in):
    """This routine returns the sum over all kx,ky,kz for a 3D object. \n
    It accounts for the implicitly deined negative-kx modes by multiplying \n
    all kx>0 modes by 2.0"""
    ksum=2.0*np.sum(np.sum(np.sum(obj_in[1:,:,:,:],axis=0),axis=0),axis=0)
    ksum=ksum+np.sum(np.sum(obj_in[0,:,:,:],axis=0),axis=0)
    #ksum=np.sum(obj_in[1:par['nkx0']-1,0:par['nky0']-2,0:par['nkz0']-2])
    #ksum=ksum*2.0
    #ksum=ksum+np.sum(obj_in[0,0:par['nky0']-2,0:par['nkz0']-2])
    return ksum 


def plot_kz_vs_hermite(start_time=-1.0,end_time=-1.0,data_out=False,fit_spectra=False,swap_endian=False,which_func='power',plot_diss=False,plot_test_exp=False,test_exp=-1.5):
    """Plots hermite spectra taken from g_out.dat \n
    start_time: start time for time average (defaults to 0.0) \n
    end_time: end time for time average (defaults to end) \n
    data_out: flag for outputting data in an ascii file--useful \n
    for replotting or quick viewing later."""
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    time=get_time_from_gout(swap_endian=swap_endian)
    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    include_kz0=True

    prefactor=np.empty(par['nv0'])
    if plot_diss:
        prefactor=par['nu']*herm_grid
        plabel='C'
    else:
        prefactor[:]=1.0
        plabel='Entropy'

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    if which_func=='power':
        p0=np.empty(2)
        p0[0]=1.0
        p0[1]=-1.0
    elif which_func=='power_exp':
        #p0=np.empty(4)
        p0=np.empty(4)
        p0[0]=1.0
        p0[1]=-1.0
        p0[2]=-par['nu']/0.2
        p0[3]=1.5


    print "hypv bound:", (par['nu']/par['hyp_v'])**(1.0/float(par['hypv_order']))*par['nv0']
    entn_sum=np.zeros((par['nv0'],11),dtype='float')
    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        gt0=read_time_step_g(i,swap_endian=swap_endian)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        #Entropy
        entn_sum[:,10]=entn_sum[:,10]+get_entropy_hermite(gt0,kzind=-1,include_kz0=include_kz0)
        for k in range(10):
            kzindex=k*par['nkz0']/20
            #print 'kzindex',kzindex
            entn_sum[:,k]=entn_sum[:,k]+get_entropy_hermite(gt0,kzind=kzindex)

    entn_sum=entn_sum/float(ntime)

    plt.loglog(herm_grid,prefactor*entn_sum[:,10],basex=10,basey=10)

    temp=prefactor*entn_sum[20,10]
    temp=temp/(herm_grid**(-1))[20]
    plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
    temp=prefactor*entn_sum[20,10]
    temp=temp/(herm_grid**(-1.5))[20]
    plt.loglog(herm_grid,2.0*temp*herm_grid**(-1.5),'--',basex=10,basey=10,label=str(-1.5))

    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title(plabel+'(kz sum)')
    plt.legend(loc='lower left')
    plt.show()


    if fit_spectra:
        print "Identify fit range from plot."
    plt.loglog(herm_grid,prefactor*entn_sum[:,10],basex=10,basey=10)
    if plot_test_exp:
        temp=prefactor*entn_sum[20,10]
        temp=temp/(herm_grid**(test_exp))[20]
        plt.loglog(herm_grid,2.0*temp*herm_grid**(test_exp),'--',basex=10,basey=10,label=str(test_exp))
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title(plabel+'(kz sum)')
    plt.legend(loc='lower left')
    plt.show()
    if fit_spectra:
        startx=float(raw_input("Enter start value for fit:"))
        endx=float(raw_input("Enter end value for fit:"))
        fit=fit_function(herm_grid,entn_sum[:,10],p0,startx=startx,endx=endx,which_func=which_func)
        plt.figure(figsize=(4.5,3.0))
        fig=plt.gcf()
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(herm_grid,prefactor*entn_sum[:,10],'x-',basex=10,basey=10)
        if startx==-1:
            sx0=0
        else:
            sx0=argmin(abs(herm_grid-startx))
        #   sx0=sx0-sx0/4
        if endx==-1:
            ex0=0
        else:
            ex0=argmin(abs(herm_grid-endx))+2
        fit_x=herm_grid[sx0:ex0]
        fit_func=2.0*fit[0]*fit_x**fit[1]
        plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
        plt.annotate(r'$\propto n$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
        plt.title(plabel)
        plt.xlabel(r'$n$',size=13)
        plt.show()


    for k in range(10):
        #print prefactor
        #print prefactor*entn_sum[:,k]
        kz0=kzgrid[k*par['nkz0']/20]
        if fit_spectra:
            print "Identify fit range from plot (for power law)."
        plt.loglog(herm_grid,prefactor*entn_sum[:,k],basex=10,basey=10,label=\
                  plabel+' (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
        plt.xlabel('Hermite n')
        plt.ylabel(plabel)
        plt.legend(loc='lower left')
        plt.show()
        if fit_spectra:
            startx=float(raw_input("Enter start value for fit:"))
            endx=float(raw_input("Enter end value for fit:"))
            fit=fit_function(herm_grid,entn_sum[:,k],p0[0:2],startx=startx,endx=endx,which_func='power')

            if which_func=='power_exp':
                print "Identify fit range for exp fit."
                startx=float(raw_input("Enter start value for fit:"))
                endx=float(raw_input("Enter end value for fit:"))
                p0[0]=fit[0]
                p0[1]=fit[1]
                p0[2]=par['nu']/kz0
                p0[3]=1.5
                fit=fit_function(herm_grid,entn_sum[:,k],p0,startx=startx,endx=endx,which_func=which_func)

            plt.figure(figsize=(4.5,3.0))
            fig=plt.gcf()
            fig.subplots_adjust(left=0.2)
            fig.subplots_adjust(bottom=0.22)
            plt.loglog(herm_grid,prefactor*entn_sum[:,k],'x-',basex=10,basey=10)
            if startx==-1:
                sx0=0
            else:
                sx0=argmin(abs(herm_grid-startx))
            if endx==-1:
                ex0=0
            else:
                ex0=argmin(abs(herm_grid-endx))+2
            fit_x=herm_grid[sx0:ex0]
            if which_func=='power':
                fit_func=2.0*fit[0]*fit_x**fit[1]
            elif which_func=='power_exp':
                #fit_func=2.0*fit[0]*fit_x**fit[1]*np.e**(fit[2]*fit_x**fit[3])
                fit_func=2.0*p0[0]*fit_x**p0[1]*np.e**(fit[0]*fit_x**fit[1])

            plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
            plt.annotate(r'$\propto n$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
            plt.title(plabel)
            plt.xlabel(r'$n$',size=13)
            plt.show()
            try_again=True
            #while try_again:
            #    plt.figure(figsize=(4.5,3.0))
            #    fig=plt.gcf()
            #    fig.subplots_adjust(left=0.2)
            #    fig.subplots_adjust(bottom=0.22)
            #    plt.loglog(herm_grid,entn_sum[:,k],'x-',basex=10,basey=10)
            #    fit[0]=raw_input("a0")
            #    fit[1]=raw_input("a1")
            #    fit[2]=raw_input("a2")
            #    fit[3]=raw_input("a3")
            #    fit_func=2.0*fit[0]*fit_x**fit[1]*np.e**(fit[2]*fit_x**fit[3])
            #    plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
            #    plt.annotate(r'$\propto n$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
            #    plt.title('Entropy')
            #    plt.xlabel(r'$n$',size=13)
            #    plt.show()
            #    try_again=raw_input("Try again?")


        if data_out:
            file_name=par['diagdir'][1:-1]+'/entropy_hermite_kz'+str(kzgrid[k*par['nkz0']/20])+'.dat'
            arr_out=np.empty((par['nv0'],2)) 
            arr_out[:,0]=herm_grid
            arr_out[:,1]=entn_sum[:,k]
            np.savetxt(file_name,arr_out)



    #plt.loglog(herm_grid,entn_sum[:,10],basex=10,basey=10)
    #plt.xlabel('Hermite n')
    #plt.ylabel('Entropy')
    #plt.title('Entropy (k_z=Avg.)')
    #plt.show()

    if data_out:
        file_name=par['diagdir'][1:-1]+'/entropy_hermite_kzavg.dat'
        arr_out[:,1]=entn_sum[:,10]
        np.savetxt(file_name,arr_out)


def plot_kz_vs_hermite_lowmem(start_time=-1.0,end_time=-1.0,data_out=False,fit_spectra=False,swap_endian=False,which_func='power',plot_diss=False):
    """Plots hermite spectra taken from g_out.dat \n
    start_time: start time for time average (defaults to 0.0) \n
    end_time: end time for time average (defaults to end) \n
    data_out: flag for outputting data in an ascii file--useful \n
    for replotting or quick viewing later."""
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    time=get_time_from_gout(swap_endian=swap_endian)
    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    include_kz0=True

    prefactor=np.empty(par['nv0'])
    if plot_diss:
        prefactor=-par['nu']*herm_grid
    else:
        prefactor[:]=1.0

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    if which_func=='power':
        p0=np.empty(2)
        p0[0]=1.0
        p0[1]=-1.0
    elif which_func=='power_exp':
        #p0=np.empty(4)
        p0=np.empty(4)
        p0[0]=1.0
        p0[1]=-1.0
        p0[2]=-par['nu']/0.2
        p0[3]=1.5


    print "hypv bound:", (par['nu']/par['hyp_v'])**(1.0/float(par['hypv_order']))*par['nv0']
    entn_sum=np.zeros((par['nv0'],11),dtype='float')
    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        for n in range(par['nv0']):
            gt0=read_time_step_gn(i,n,swap_endian=swap_endian)
            gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
            #Entropy
            entn_sum[n,10]=entn_sum[n,10]+get_entropy_hermiten(gt0,kzind=-1,include_kz0=include_kz0)
            for k in range(10):
                kzindex=k*par['nkz0']/20
                #print 'kzindex',kzindex
                entn_sum[n,k]=entn_sum[n,k]+get_entropy_hermiten(gt0,kzind=kzindex)

    entn_sum=entn_sum/float(ntime)
    for k in range(10):
        kz0=kzgrid[k*par['nkz0']/20]
        if fit_spectra:
            print "Identify fit range from plot (for power law)."
        plt.loglog(herm_grid,prefactor*entn_sum[:,k],basex=10,basey=10,label=\
                  'Entropy (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
        plt.xlabel('Hermite n')
        plt.ylabel('Entropy')
        plt.legend(loc='lower left')
        plt.show()
        if fit_spectra:
            startx=float(raw_input("Enter start value for fit:"))
            endx=float(raw_input("Enter end value for fit:"))
            fit=fit_function(herm_grid,entn_sum[:,k],p0[0:2],startx=startx,endx=endx,which_func='power')

            if which_func=='power_exp':
                print "Identify fit range for exp fit."
                startx=float(raw_input("Enter start value for fit:"))
                endx=float(raw_input("Enter end value for fit:"))
                p0[0]=fit[0]
                p0[1]=fit[1]
                p0[2]=-par['nu']/kz0
                p0[3]=1.5
                fit=fit_function(herm_grid,entn_sum[:,k],p0,startx=startx,endx=endx,which_func=which_func)

            plt.figure(figsize=(4.5,3.0))
            fig=plt.gcf()
            fig.subplots_adjust(left=0.2)
            fig.subplots_adjust(bottom=0.22)
            plt.loglog(herm_grid,prefactor*entn_sum[:,k],'x-',basex=10,basey=10)
            if startx==-1:
                sx0=0
            else:
                sx0=argmin(abs(herm_grid-startx))
            if endx==-1:
                ex0=0
            else:
                ex0=argmin(abs(herm_grid-endx))+2
            fit_x=herm_grid[sx0:ex0]
            if which_func=='power':
                fit_func=2.0*fit[0]*fit_x**fit[1]
            elif which_func=='power_exp':
                #fit_func=2.0*fit[0]*fit_x**fit[1]*np.e**(fit[2]*fit_x**fit[3])
                fit_func=2.0*p0[0]*fit_x**p0[1]*np.e**(fit[0]*fit_x**fit[1])

            plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
            plt.annotate(r'$\propto n$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
            plt.title('Entropy')
            plt.xlabel(r'$n$',size=13)
            plt.show()
            try_again=True
            #while try_again:
            #    plt.figure(figsize=(4.5,3.0))
            #    fig=plt.gcf()
            #    fig.subplots_adjust(left=0.2)
            #    fig.subplots_adjust(bottom=0.22)
            #    plt.loglog(herm_grid,entn_sum[:,k],'x-',basex=10,basey=10)
            #    fit[0]=raw_input("a0")
            #    fit[1]=raw_input("a1")
            #    fit[2]=raw_input("a2")
            #    fit[3]=raw_input("a3")
            #    fit_func=2.0*fit[0]*fit_x**fit[1]*np.e**(fit[2]*fit_x**fit[3])
            #    plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
            #    plt.annotate(r'$\propto n$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
            #    plt.title('Entropy')
            #    plt.xlabel(r'$n$',size=13)
            #    plt.show()
            #    try_again=raw_input("Try again?")


        if data_out:
            file_name=par['diagdir'][1:-1]+'/entropy_hermite_kz'+str(kzgrid[k*par['nkz0']/20])+'.dat'
            arr_out=np.empty((par['nv0'],2)) 
            arr_out[:,0]=herm_grid
            arr_out[:,1]=entn_sum[:,k]
            np.savetxt(file_name,arr_out)

    if fit_spectra:
        print "Identify fit range from plot."
    plt.loglog(herm_grid,prefactor*entn_sum[:,10],basex=10,basey=10)
              #'Entropy (k_z=Avg.)')
    plt.xlabel('Hermite n')
    plt.ylabel('Entropy')
    plt.title('Entropy(kz sum)')
    plt.legend(loc='lower left')
    plt.show()
    if fit_spectra:
        startx=float(raw_input("Enter start value for fit:"))
        endx=float(raw_input("Enter end value for fit:"))
        fit=fit_function(herm_grid,entn_sum[:,10],p0,startx=startx,endx=endx,which_func=which_func)
        plt.figure(figsize=(4.5,3.0))
        fig=plt.gcf()
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(herm_grid,prefactor*entn_sum[:,10],'x-',basex=10,basey=10)
        if startx==-1:
            sx0=0
        else:
            sx0=argmin(abs(herm_grid-startx))
        #   sx0=sx0-sx0/4
        if endx==-1:
            ex0=0
        else:
            ex0=argmin(abs(herm_grid-endx))+2
        fit_x=herm_grid[sx0:ex0]
        fit_func=2.0*fit[0]*fit_x**fit[1]
        plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
        plt.annotate(r'$\propto n$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
        plt.title('Entropy')
        plt.xlabel(r'$n$',size=13)
        plt.show()



    #plt.loglog(herm_grid,entn_sum[:,10],basex=10,basey=10)
    #plt.xlabel('Hermite n')
    #plt.ylabel('Entropy')
    #plt.title('Entropy (k_z=Avg.)')
    #plt.show()

    if data_out:
        file_name=par['diagdir'][1:-1]+'/entropy_hermite_kzavg.dat'
        arr_out[:,1]=entn_sum[:,10]
        np.savetxt(file_name,arr_out)



def plot_prep2d(arr2d):
    """This routine shifts the 3D array so that the indices correspond \n
    to -kmax==>kmax instead of the FFT format defined by \'get_grids\'."""
    arr2d=np.roll(arr2d,par['nky0']/2-1,axis=1)
    arr2d0=np.zeros((2*par['nkx0']-1,par['nky0']))
    #print shape(arr2d0)
    #print shape(arr2d)
    arr2d0[par['nkx0']-1:,:]=arr2d[:,:]
    for i in range(par['nkx0']-1):
        for j in range(1,par['nky0']):
            arr2d0[i,j-1]=arr2d[par['nkx0']-1-i,par['nky0']-1-j]
    return arr2d0

def les_prep3d(arr3d,kx_red_l,kx_red_h,ky_red_ll,ky_red_hh,kz_red_ll,kz_red_hh):
    """This routine shifts the 3D array so that the indices correspond \n
    to -kmax==>kmax instead of the FFT format defined by \'get_grids\'."""
    arr3d=np.roll(arr3d,par['nky0']/2-1,axis=1)
    arr3d=np.roll(arr3d,par['nkz0']/2-1,axis=2)
    arr3d0=np.zeros((2*par['nkx0']-1,par['nky0'],par['nkz0']))
    arr3d1=np.zeros((2*par['nkx0']-1,par['nky0'],par['nkz0']))
   # print shape(arr3d0)
   # print shape(arr3d)

    arr3d0[par['nkx0']-1:,:,:]=arr3d[:,:,:]
    for i in range(par['nkx0']-1):
        for j in range(1,par['nky0']):
            arr3d0[i,j-1,:]=arr3d[par['nkx0']-1-i,par['nky0']-1-j,:]
    arr3d1 = arr3d0[kx_red_l:kx_red_h,ky_red_ll:ky_red_hh,kz_red_ll:kz_red_hh]
    return arr3d1

def les_prep1d(arr3d,nkx0_red,nky0_red,nkz0_red):    
   # arr3d_kx = np.sum(np.sum(arr3d,2),1)
   # arr3d_ky = np.sum(np.sum(arr3d,2),0)
   # arr3d_kz = np.sum(np.sum(arr3d,1),0)
 
    arr3d_kx = np.sum(np.sum(arr3d,2),1)
    arr3d_ky = np.sum(np.sum(arr3d[0:par['nkx0'],:,:],2),0)
    arr3d_kz = np.sum(np.sum(arr3d[0:par['nkx0'],:,:],1),0)
    
    
    arr3d_kx_pos = np.zeros((nkx0_red))
    arr3d_ky_pos = np.zeros((nky0_red/2))
    arr3d_kz_pos = np.zeros((nkz0_red/2))
#    print arr3d_kx_pos.shape
#    print arr3d_kx.shape
    arr3d_kx_pos[0] = arr3d_kx[nkx0_red-1]
    arr3d_ky_pos[0] = arr3d_ky[nky0_red/2-1]
    arr3d_kz_pos[0] = arr3d_kz[nkz0_red/2-1]
    for i in range(1,nkx0_red):
        arr3d_kx_pos[i] = arr3d_kx[nkx0_red -1  +i] + arr3d_kx[nkx0_red  -1 -i]
#        print i, arr3d_kx[nkx0_red -1 +i], arr3d_kx[nkx0_red -1 -i]
#        print arr3d_kx_pos[i]

    for i in range(1,nky0_red/2):
        arr3d_ky_pos[i] = arr3d_ky[nky0_red/2 -1  +i] + arr3d_ky[nky0_red/2  -1 -i]
#        print i, arr3d_ky[nky0_red/2 -1 +i], arr3d_ky[nky0_red/2 -1 -i]
#        print arr3d_ky_pos[i]

    for i in range(1,nkz0_red/2):
        arr3d_kz_pos[i] = arr3d_kz[nkz0_red/2 -1  +i] + arr3d_kz[nkz0_red/2  -1 -i]
#        print i, arr3d_kz[nkz0_red/2 -1 +i], arr3d_kz[nkz0_red/2 -1 -i]
#        print arr3d_kz_pos[i]

    return arr3d_kx_pos, arr3d_ky_pos, arr3d_kz_pos

def test_energy_equation(ikx,iky,start_time=-1.0,end_time=-1.0,data_from_gk=True,\
                 g_transform=True):

    if data_from_gk:
        time=get_time_from_gkfile(ikx,iky)
        #itime=np.arange(len(time))
    else:
        time=get_time_from_gout()
        #itime=np.arange(len(time))

    time0=np.empty(0)
    itime=np.empty(0)
    #print time
    #print itime

    tmax=-1.0
    for i in range(len(time)):
        if time[i]>tmax:
            tmax=time[i]
            time0=np.append(time0,time[i])
            itime=np.append(itime,i)
        #else:
        #    time=np.delete(time,i)
        
    time=time0

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop
        
    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    ntime=iend-istart+1
    gn=np.empty((par['nv0'],ntime),dtype='complex')
    gkz=np.empty((par['nkz0'],par['nv0'],ntime),dtype='complex')
    PM=np.empty((par['nkz0'],par['nv0'],ntime),dtype='float')
    #for i in itime[istart:iend+1]:
    for i in range(istart,iend+1):
        print i
        print 'time=',time[i],' of ',time[iend]
        it=i-istart
        if data_from_gk:
            gt0=read_time_step_gkfile(ikx,iky,itime[i])
            gt0=np.reshape(gt0,(par['nkz0'],par['nv0']),order='F')
            gkz[:,:,it]=gt0
            for k in range(par['nkz0']):
                PM[k,:,it]=get_henergy_single_k(gkz[k,:,it],gkz[k,:,it],par['kxmin']*ikx,par['kymin']*iky,kzgrid[k],3)
            if g_transform:
                for n in range(par['nv0']):
                    gkz[0,n,it]=(-1.0J)**n*gkz[0,n,it]
                    for k in range(1,par['nkz0']):
                        gkz[k,n,it]=(-1.0J*kzgrid[k]/np.abs(kzgrid[k]))**n*gkz[k,n,it]
            gn[:,it]=np.sum(gkz[:,:,it],axis=0)
        else:
            gt0=read_time_step_g(itime[i])
            gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
            gkz[:,:,it]=gt0[ikx,iky,:,:]
            if g_transform:
                for n in range(par['nv0']):
                    gkz[0,n,it]=(-1.0J)**n*gkz[0,n,it]
                    for k in range(1,par['nkz0']):
                        gkz[k,n,it]=(-1.0J*kzgrid[k]/np.abs(kzgrid[k]))**n*gkz[k,n,it]
            gn[:,it]=np.sum(gkz[:,:,it],axis=0)
        
    PM=np.sum(PM,axis=2)/float(ntime)
    Ekz=np.sum(np.conj(gkz)*gkz,axis=2)/float(ntime)
    kzn=np.empty(par['nv0'],dtype='float')

    En=np.conj(gkz)*gkz
    En=np.sum(En,axis=2)/float(ntime)
    En=np.sum(En,axis=0)

    n0p5=10.0*herm_grid**(-0.5)
    n1p0=10.0*herm_grid**(-1.0)
    n1p5=10.0*herm_grid**(-1.5)
    plt.loglog(np.arange(par['nv0']),En,label=str(i),basex=10,basey=10)
    plt.loglog(np.arange(par['nv0']),n0p5,'--',label='-0.5',basex=10,basey=10)
    plt.loglog(np.arange(par['nv0']),n1p0,':',label='-1.0',basex=10,basey=10)
    plt.loglog(np.arange(par['nv0']),n1p5,'-.',label='-1.5',basex=10,basey=10)
    plt.legend() 
    plt.title(r'$\varepsilon_n$')
    plt.show()

    J0a=get_gyroavg()
    Gam0=get_gamma0()


    akzn=np.empty((par['nkz0'],par['nv0']),dtype='float')
    delp=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    delm=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    hkzn=np.zeros((par['nkz0'],par['nv0'],ntime),dtype='complex')
    sigmap=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    sigmam=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    for n in range(par['nv0']):
        for k in range(par['nkz0']):
            akzn[k,n]=np.real(np.sqrt(np.sum(np.conj(gkz[k,n,:])*gkz[k,n,:])/float(ntime)))
            if k != par['nkz0']/2:
                hkzn[k,n,:]=gkz[k,n,:]/akzn[k,n]
    for n in range(1,par['nv0']-1):
        delp[:,n]=akzn[:,n]-akzn[:,n+1]
        delm[:,n]=-akzn[:,n]+akzn[:,n-1]
    anspect=np.sum(akzn**2,axis=0)
    plt.loglog(anspect,'x',label='an2',basex=10,basey=10)
    plt.loglog(En,label='g2',basex=10,basey=10)
    plt.legend()
    plt.show()
    for n in range(1,par['nv0']-1):
        for k in range(par['nkz0']):
            sigmap[k,n]=np.real(np.sum(np.conj(hkzn[k,n,:])*hkzn[k,n+1,:]))/float(ntime)
            sigmam[k,n]=np.real(np.sum(np.conj(hkzn[k,n,:])*hkzn[k,n-1,:]))/float(ntime)
            print "n,k,sigmap",n,k,sigmap[k,n]

    sigmap_n=np.sum(np.abs(sigmap),axis=0)/float(par['nkz0'])
    sigmam_n=np.sum(np.abs(sigmam),axis=0)/float(par['nkz0'])
    #print "sigmap_n",sigmap_n
    #print "sigmam_n",sigmam_n
    plt.plot(herm_grid,sigmap_n,label='sigma_plus')
    plt.plot(herm_grid,sigmam_n,label='sigma_minus')
    plt.legend()
    plt.show()
    plt.plot(herm_grid,sigmap_n-sigmam_n,label='sigma_plus-sigma_minus')
    plt.legend()
    plt.show()

    ######Calculating III from 4/29/13 notes
    EqIII=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    EqIV=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    EqV=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    EqVI=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    EqVII=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    test=np.zeros((par['nkz0'],par['nv0']),dtype='float')
    for k in range(par['nkz0']):
        for n in range(1,par['nv0']-1):
            EqIII[k,n]=np.pi**0.5*np.abs(kzgrid[k])*(-sigmam[k,n]*np.sqrt(herm_grid[n])*akzn[k,n]**2\
                                          -sigmam[k,n]*np.sqrt(herm_grid[n])*akzn[k,n]*delm[k,n]\
                                          +sigmap[k,n]*np.sqrt(herm_grid[n+1])*akzn[k,n]**2\
                                          -sigmap[k,n]*np.sqrt(herm_grid[n+1])*akzn[k,n]*delp[k,n])
            EqIV[k,n]=np.pi**0.5*np.abs(kzgrid[k])*(-sigmam[k,n]*np.sqrt(herm_grid[n])*akzn[k,n]**2\
                                          +sigmap[k,n]*np.sqrt(herm_grid[n+1])*akzn[k,n+1]**2)
            EqV[k,n]=np.pi**0.5*np.abs(kzgrid[k])*sigmam[k,n]*(-np.sqrt(herm_grid[n])*akzn[k,n]**2\
                                          +np.sqrt(herm_grid[n+1])*akzn[k,n+1]**2)
            EqVI[k,n]=-np.pi**0.5*np.abs(kzgrid[k])*(-np.sqrt(herm_grid[n])*akzn[k,n]**2\
                                          +np.sqrt(herm_grid[n+1])*akzn[k,n+1]**2)
            EqVII[k,n]=np.pi**0.5*np.abs(kzgrid[k])*(-np.sqrt(herm_grid[n])*akzn[k,n]**2\
                                          -np.sqrt(herm_grid[n])*akzn[k,n]*delm[k,n]\
                                          +np.sqrt(herm_grid[n+1])*akzn[k,n]**2\
                                          -np.sqrt(herm_grid[n+1])*akzn[k,n]*delp[k,n])
            test[k,n]=np.pi**0.5*np.abs(kzgrid[k])*(-sigmam[k,n]*np.sqrt(herm_grid[n])*akzn[k,n]**2\
                                          -sigmam[k,n]*np.sqrt(herm_grid[n])*akzn[k,n]*delm[k,n]\
                                          +sigmam[k,n]*np.sqrt(herm_grid[n+1])*akzn[k,n]**2\
                                          -sigmam[k,n]*np.sqrt(herm_grid[n+1])*akzn[k,n]*delp[k,n])

    plt.loglog(np.abs(np.sum(PM,axis=0)),label='PM',basex=10,basey=10)
    plt.loglog(np.abs(np.sum(EqIII,axis=0)),'-',label='EqIII',basex=10,basey=10)
    plt.loglog(np.abs(np.sum(EqIV,axis=0)),'-',label='EqIV',basex=10,basey=10)
    plt.loglog(np.abs(np.sum(EqV,axis=0)),'-',label='EqV',basex=10,basey=10)
    plt.loglog(np.abs(np.sum(EqVI,axis=0)),'-',label='EqVI',basex=10,basey=10)
    plt.loglog(np.abs(np.sum(EqVII,axis=0)),'-',label='EqVII',basex=10,basey=10)
    plt.legend()
    plt.show()

    plt.plot(np.sum(PM,axis=0),label='PM')
    plt.plot(np.sum(EqIII,axis=0),'x',label='EqIII')
    plt.legend()
    plt.show()
    plt.plot(np.sum(PM,axis=0),label='PM')
    plt.plot(np.sum(EqIV,axis=0),'x',label='EqIV')
    plt.legend()
    plt.show()
    plt.plot(np.sum(PM,axis=0),label='PM')
    plt.plot(np.sum(EqV,axis=0),'x',label='EqV')
    plt.legend()
    plt.show()
    plt.plot(np.sum(PM,axis=0),label='PM')
    plt.plot(np.sum(EqVI,axis=0),'x',label='EqVI')
    plt.legend()
    plt.show()
    plt.plot(np.sum(PM,axis=0)/np.sum(EqVI,axis=0),label='PM/EqVI')
    plt.legend()
    plt.show()




    #plt.plot(np.sum(PM,axis=0),label='PM')
    #plt.plot(np.sum(test,axis=0),'x',label='Test')
    #plt.legend()
    #plt.show()
    #plt.plot(np.sum(PM,axis=0)/np.sum(EqIII,axis=0),label='PM/EqIII')
    #plt.legend()
    #plt.show()
    #plt.title("ikz=3")
    #plt.plot(PM[3,:],label='PM')
    #plt.plot(EqIII[3,:],label='EqIII')
    #plt.legend()
    #plt.show()
    #plt.title("ikz=3")
    #plt.plot(PM[3,:]/EqIII[3,:],label='PM/EqIII')
    #plt.legend()
    #plt.show()

    ######Testing <kz> vs <kz sigma>
    kzsigma_n=np.zeros(par['nv0'],dtype='float')
    for n in range(par['nv0']):
        kzn[n]=np.sum(Ekz[:,n]*np.abs(kzgrid))/np.sum(Ekz[:,n])
        kzsigma_n[n]=np.sum(sigmam[:,n]*Ekz[:,n]*np.abs(kzgrid))/np.sum(Ekz[:,n])
            
    plt.plot(herm_grid,kzn,label='kzn')
    plt.plot(herm_grid,np.abs(kzsigma_n),label='abs(kzsigma_n)')
    plt.plot(herm_grid,kzn/np.sqrt(herm_grid),label='kzn/sqrt(n)')
    plt.plot(herm_grid,np.abs(kzsigma_n)/np.sqrt(herm_grid),label='abs(kzsigma_n/sqrt(n))')
    plt.legend()
    plt.show()

    plt.plot(herm_grid,np.abs(sigmap_n*kzn),label='sigmap*kzn')
    plt.plot(herm_grid,np.abs(sigmap_n*kzn)/np.sqrt(herm_grid),label='sigmap*kzn/sqrt(n)')
    plt.legend()
    plt.show()

def time_scales_nl(ikx,iky,start_time=-1.0,end_time=-1.0,show_plots=False):

    time=get_time_from_gkfile(ikx,iky,read_nl=True)

    time0=np.empty(0)
    itime=np.empty(0)
    #print time
    #print itime

    tmax=-1.0
    for i in range(len(time)):
        if time[i]>tmax:
            tmax=time[i]
            time0=np.append(time0,time[i])
            itime=np.append(itime,i)
        
    time=time0

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop
        
    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    ntime=iend-istart+1
    gn=np.empty((par['nv0'],ntime),dtype='complex')
    gkz=np.empty((par['nkz0'],par['nv0'],ntime),dtype='complex')
    #for i in itime[istart:iend+1]:
    for i in range(istart,iend+1):
        print i
        print 'time=',time[i],' of ',time[iend]
        it=i-istart
        gt0=read_time_step_gkfile(ikx,iky,itime[i])
        gt0=np.reshape(gt0,(par['nkz0'],par['nv0']),order='F')
        gkz[:,:,it]=gt0
        gn[:,it]=np.sum(gkz[:,:,it],axis=0)

    ptime0=time[istart:iend+1]
    corr_time=np.empty(par['nv0'],dtype='float')
    for n in range(par['nv0']):
        print "n=",n
        phicorr,tau,corr_time[n]=my_corr_func_complex(gn[n,:],gn[n,:],ptime0,show_plot=False)
   
    if show_plots:
        plt.plot(herm_grid,1/corr_time,label='1/corr time')
        plt.legend()
        plt.xlabel('n')
        plt.show()


    for n in range(10):
        omega,fs=get_frequency_spectrum(time[istart:iend+1],gn[n,:])        
        if show_plots:
            plt.plot(omega,np.conj(fs)*fs,label='frequency n='+str(n))
    if show_plots:
        plt.legend()
        plt.xlabel('omega')
        plt.show()



    #for n in range(par['nv0']):
    #    omega,fs=get_frequency_spectrum(time[istart:iend+1],gn[n,:])        
    #    plt.plot(omega,np.conj(fs)*fs,label='frequency n='+str(n))
    #    plt.legend()
    #    plt.xlabel('omega')
    #    plt.show()

    return corr_time


def time_scales_n(ikx,iky,start_time=-1.0,end_time=-1.0,show_plot_alln=False,data_from_gk=True,\
                  one_over_sqrtn=True,one_over_n=False,fit_spectra=False,\
                 plot_kz_n=False, plot_sigma_n=False,g_transform=True,compare_linear=False,\
                 use_nl=False,show_plots=True,tau_max=15.0):

    if data_from_gk:
        time=get_time_from_gkfile(ikx,iky)
        #itime=np.arange(len(time))
    else:
        time=get_time_from_gout()
        #itime=np.arange(len(time))

    time0=np.empty(0)
    itime=np.empty(0)
    #print time
    #print itime

    tmax=-1.0
    for i in range(len(time)):
        if time[i]>tmax:
            tmax=time[i]
            time0=np.append(time0,time[i])
            itime=np.append(itime,i)
        #else:
        #    time=np.delete(time,i)
        
    time=time0
    #print time
    #print itime
    #plt.plot(time)
    #plt.show()
    #plt.plot(itime)
    #plt.show()

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    #kxgrid_full=np.empty(par['nkx0']*2-1)
    #for i in range(par['nkx0']):
    #    kxgrid_full[i+par['nkx0']-1]=kxgrid[i]
    #for i in range(par['nkx0']-1):
    #    kxgrid_full[i]=-1.0*kxgrid[par['nkx0']-1-i]
    #print kxgrid_full

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop
        
    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    ntime=iend-istart+1
    gn=np.empty((par['nv0'],ntime),dtype='complex')
    gkz=np.empty((par['nkz0'],par['nv0'],ntime),dtype='complex')
    #for i in itime[istart:iend+1]:
    for i in range(istart,iend+1):
        print i
        print 'time=',time[i],' of ',time[iend]
        it=i-istart
        if data_from_gk:
            gt0=read_time_step_gkfile(ikx,iky,itime[i])
            gt0=np.reshape(gt0,(par['nkz0'],par['nv0']),order='F')
            gkz[:,:,it]=gt0
            if g_transform:
                for n in range(par['nv0']):
                    gkz[0,n,it]=(-1.0J)**n*gkz[0,n,it]
                    for k in range(1,par['nkz0']):
                        gkz[k,n,it]=(-1.0J*kzgrid[k]/np.abs(kzgrid[k]))**n*gkz[k,n,it]
            gn[:,it]=np.sum(gkz[:,:,it],axis=0)
        else:
            gt0=read_time_step_g(itime[i])
            gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
            gkz[:,:,it]=gt0[ikx,iky,:,:]
            if g_transform:
                for n in range(par['nv0']):
                    gkz[0,n,it]=(-1.0J)**n*gkz[0,n,it]
                    for k in range(1,par['nkz0']):
                        gkz[k,n,it]=(-1.0J*kzgrid[k]/np.abs(kzgrid[k]))**n*gkz[k,n,it]
            gn[:,it]=np.sum(gkz[:,:,it],axis=0)
        
    Ekz=np.real(np.sum(np.conj(gkz)*gkz,axis=2))
    ptime0=time[istart:iend+1]
    corr_time=np.empty(par['nv0'],dtype='float')
    if use_nl:
        ct_nl=time_scales_nl(ikx,iky,start_time=start_time)
    else:
        ct_nl=np.zeros(par['nv0'])
    kzn=np.empty(par['nv0'],dtype='float')
    kzsqrtn=np.empty(par['nv0'],dtype='float')
    kzon=np.empty(par['nv0'],dtype='float')

    g_nspect=np.conj(gkz)*gkz
    g_nspect=np.sum(g_nspect,axis=2)/float(ntime)
    g_nspect=np.sum(g_nspect,axis=0)

    n0p5=500.0*herm_grid**(-0.5)
    n1p0=500.0*herm_grid**(-1.0)
    n1p5=500.0*herm_grid**(-1.5)
    if show_plots:
        plt.loglog(np.arange(par['nv0']),g_nspect,label=str(i),basex=10,basey=10)
        plt.loglog(np.arange(par['nv0']),n0p5,'--',label='-0.5',basex=10,basey=10)
        plt.loglog(np.arange(par['nv0']),n1p0,':',label='-1.0',basex=10,basey=10)
        plt.loglog(np.arange(par['nv0']),n1p5,'-.',label='-1.5',basex=10,basey=10)
        plt.legend() 
        plt.title(r'$\varepsilon_n$')
        plt.show()
    np.savetxt(par['diagdir'][1:-1]+'/ts_hspect_kx'+str(ikx)+'_ky'+str(iky)+'.dat',g_nspect)

    J0a=get_gyroavg()
    Gam0=get_gamma0()
    #for i in range(par['nkx0']):
    #    for j in range(par['nky0']):
    #phi[i,j,:]=np.pi**(0.25)*J0a[i,j]*g_in[i,j,:,0]\
    #        /(par['Ti0Te']+1.0-Gam0[i,j])
    #print 'sum(phi) in get_phi',np.sum(phi)
    if not show_plots:
        show_plot_alln=False

    for n in range(par['nv0']):
        print "n=",n
        phicorr,tau,corr_time[n]=my_corr_func_complex(gn[n,:],gn[n,:],ptime0,show_plot=show_plot_alln)
        phicorr=phicorr/phicorr[0]
        if n==0:
            i_tau_max=np.argmin(abs(tau_max-tau))
            corr_out=np.empty((i_tau_max,par['nv0']))
        corr_out[:,n]=phicorr[:i_tau_max] 
        kzn[n]=np.sum(Ekz[:,n]*np.abs(kzgrid))/np.sum(Ekz[:,n])
        kzsqrtn[n]=herm_grid[n]**(-0.5)*np.sum(Ekz[:,n]*np.abs(kzgrid))/np.sum(Ekz[:,n])
        kzon[n]=herm_grid[n]**(-1)*np.sum(Ekz[:,n]*np.abs(kzgrid))/np.sum(Ekz[:,n])
        print "1/tau",1.0/corr_time[n]
        print "<kz>",kzn[n]
        if show_plots and show_plot_alln:
            plt.plot(kzgrid,np.conj(Ekz[:,n])*Ekz[:,n])
            ax=plt.axis()
            plt.vlines(kzn[n],ax[2],ax[3])
            plt.xlabel('kz')
            plt.show()

    np.savetxt(par['diagdir'][1:-1]+'/ts_tcorr_kx'+str(ikx)+'_ky'+str(iky)+'.dat',corr_out)
    np.savetxt(par['diagdir'][1:-1]+'/ts_tau_kx'+str(ikx)+'_ky'+str(iky)+'.dat',tau[:i_tau_max])
    #plt.contourf(tau[:i_tau_max],np.arange(['nv0']),corr_out,50)
    if show_plots:
        plt.contourf(corr_out,50)
        plt.plot(1.0/kzn)
        plt.show()
    
    Ekz0=np.zeros((par['nkz0']/2,par['nv0']),dtype='float')
    Ekz0[0,:]=Ekz[0,:]
    for k in range(1,par['nkz0']/2):
        print k
        Ekz0[k,:]=Ekz[k,:]+Ekz[par['nkz0']-k,:]

    print "Outputting ts_kzpect"
    np.savetxt(par['diagdir'][1:-1]+'/ts_kzspect_kx'+str(ikx)+'_ky'+str(iky)+'.dat',Ekz0)

    phi=np.empty((par['nkz0'],ntime),dtype='complex')
    phi=np.pi**(0.25)*J0a[ikx,iky]*gkz[:,0,:]\
            /(par['Ti0Te']+1.0-Gam0[ikx,iky])
    #phi0=np.sum(np.abs(np.sum(phi,axis=0)))/float(ntime)
    phisq=np.sum(np.conj(phi)*phi,axis=1)/float(ntime)
    #phibarsq=np.sum(np.e**(-kygrid[iky]**2)*np.conj(phi)*phi,axis=1)/float(ntime)
    kzphi=np.sum(phisq*np.abs(kzgrid))/np.sum(phisq)
    #kzphibar=np.sum(phibarsq*np.abs(kzgrid))/np.sum(phibarsq)
    kzphi_arr=np.empty(par['nv0'])
    #kzphibar_arr=np.empty(par['nv0'])
    kzphi_arr[:]=np.real(kzphi)
    #kzphibar_arr[:]=kzphibar

    #omega_drive=np.arange(par['nv0'],dtype='float')
    #omega_drive[:]=par['omt']*iky*par['kymin']*2.0**(-0.5)*np.pi**(-0.25)
    #omega_nl1=np.arange(par['nv0'],dtype='float')
    #omega_nl1[:]=(kygrid[iky])**2*phi0

    out_data=np.empty((par['nv0'],4),dtype='float')
    out_data[:,0]=1.0/corr_time
    out_data[:,1]=kzn
    out_data[:,2]=par['nu']*herm_grid+par['hyp_v']*(herm_grid/float(par['nv0']-1))**par['hypv_order']
    if use_nl:
        out_data[:,3]=1.0/ct_nl
    else:
        out_data[:,3]=0.0
    np.savetxt(par['diagdir'][1:-1]+'/ts_cbtest_kx'+str(ikx)+'_ky'+str(iky)+'.dat',out_data)
    if show_plots:
        plt.title(r'$k_x \rho_i = $'+str(ikx*par['kxmin'])+','+r'$k_y \rho_i = $'+str(iky*par['kymin']),size=18)
        plt.plot(1.0/corr_time,'o',markeredgewidth=0,label=r'$\tau_{NL}^{-1}(v_{t}/L)$')
        if use_nl:
            plt.plot(1.0/ct_nl,'o',markeredgewidth=0,label=r'$\tau_{NL0}^{-1}(v_{t}/L)$')
        plt.plot(kzn,'x-',label=r'$<k_z v_{t}>(v_{t}/L)$')
        #plt.plot(0.5*np.sqrt(herm_grid),'-',label=r'$0.5 \sqrt{n}$')
        #plt.plot(omega_nl1,label='kperp^2 phi0')
        if one_over_sqrtn:
            plt.plot(kzsqrtn,'x-',label=r'$<k_z>/\sqrt{n}$')
        if one_over_n:
            plt.plot(kzon,'x-',label=r'$<k_z>/n$')
        #plt.plot(kzon,'x-',label=r'$<k_z>/n$')
        plt.plot(kzphi_arr,'-',label=r'$<k_z>_\phi$')
        #plt.plot(kzphibar_arr,'-',label=r'$<k_z>_\bar{\phi}$')
        #plt.plot(omega_drive,label='omega drive')
        #plt.plot(0.5*kzn,'+-',label=r'$0.5<k_z v_{t}>(v_{t}/L)$')
        plt.plot(par['nu']*herm_grid+par['hyp_v']*(herm_grid/float(par['nv0']-1))**par['hypv_order'],'-',label=r'$\nu n$')
        ax=plt.axis()
        plt.axis((ax[0],ax[1],0,ax[3]))
        #plt.plot(0.5*kzn,label=r'$0.5 <k_z v_{t}>(v_{t}/L)$')
        plt.xlabel('n',size=18)
        plt.ylabel(r'$\omega(v_{t}/L)$',size=18)
        plt.legend(loc='upper left')
        plt.show()

    if show_plots and plot_sigma_n:
        akzn=np.empty((par['nkz0'],par['nv0']),dtype='float')
        hkzn=np.zeros((par['nkz0'],par['nv0'],ntime),dtype='complex')
        sigmap=np.zeros((par['nkz0'],par['nv0']),dtype='float')
        sigmam=np.zeros((par['nkz0'],par['nv0']),dtype='float')
        for n in range(par['nv0']):
            for k in range(par['nkz0']):
                akzn[k,n]=np.real(np.sqrt(np.sum(np.conj(gkz[k,n,:])*gkz[k,n,:])/float(ntime)))
                if k != par['nkz0']/2:
                    hkzn[k,n,:]=gkz[k,n,:]/akzn[k,n]
        anspect=np.sum(akzn**2,axis=0)
        plt.loglog(anspect,'x',label='an2',basex=10,basey=10)
        plt.loglog(g_nspect,label='g2',basex=10,basey=10)
        plt.legend()
        plt.show()
        for n in range(1,par['nv0']-1):
            for k in range(par['nkz0']):
                sigmap[k,n]=np.real(np.sum(np.conj(hkzn[k,n,:])*hkzn[k,n+1,:]))/float(ntime)
                sigmam[k,n]=np.real(np.sum(np.conj(hkzn[k,n,:])*hkzn[k,n-1,:]))/float(ntime)
                print "n,k,sigmap",n,k,sigmap[k,n]
        if compare_linear:
            ik_lin_compare=3
            k=ik_lin_compare
            plt.plot(herm_grid,sigmam[k,:],label='sigma-_ev')
            plt.plot(herm_grid,sigmap[k,:],label='sigma+_ev')
            plt.legend()
            plt.show()
            plt.plot(herm_grid,sigmap[k,:]-sigmam[k,:],label='del_sigma')
            plt.legend()
            plt.show()

        sigmap_n=np.sum(sigmap,axis=0)
        sigmam_n=np.sum(sigmam,axis=0)
        print "sigmap_n",sigmap_n
        print "sigmam_n",sigmam_n
        plt.plot(herm_grid,sigmap_n,label='sigma_plus')
        plt.plot(herm_grid,sigmam_n,label='sigma_minus')
        plt.legend()
        plt.show()
        plt.plot(herm_grid,sigmap_n-sigmam_n,label='sigma_plus-sigma_minus')
        plt.legend()
        plt.show()
        plt.plot(herm_grid,np.abs(sigmap_n*kzn),label='sigmap*kzn')
        plt.plot(herm_grid,np.abs(sigmap_n*kzn)/np.sqrt(herm_grid),label='sigmap*kzn/sqrt(n)')
        plt.legend()
        plt.show()
        #plt.plot(np.real(hn[0,:]),label='Re hn n=0')
        #lt.plot(np.imag(hn[0,:]),label='Im hn n=0')
        #lt.legend()
        #lt.show()
        if compare_linear:
            hev=np.zeros(par['nv0'],dtype='complex')
            sigmap_ev=np.zeros(par['nv0'])
            sigmam_ev=np.zeros(par['nv0'])
            mat=get_lin_matrix(par['kxmin']*ikx,par['kymin']*iky,par['kzmin']*ik_lin_compare)
            gam,ome,evec=get_ev_spectrum(mat)
            vec=evec[:,0]
            plt.plot(np.abs(vec))
            plt.title('evec')
            plt.show()
            for n in range(par['nv0']):
                hev[n]=vec[n]/np.real(np.sqrt(np.conj(vec[n])*vec[n])) 
            for n in range(1,par['nv0']-1):
                sigmam_ev[n]=np.real(np.conj(hev[n])*hev[n-1])
                sigmap_ev[n]=np.real(np.conj(hev[n])*hev[n+1])
            plt.plot(herm_grid,sigmam_ev,label='sigma-_ev')
            plt.plot(herm_grid,sigmap_ev,label='sigma+_ev')
            plt.legend()
            plt.show()
            plt.plot(herm_grid,sigmap_ev-sigmam_ev,label='del_sigma')
            plt.legend()
            plt.show()
                
                

    if show_plots and plot_kz_n:
        plt.title(r'$k_x \rho_i = $'+str(ikx*par['kxmin'])+','+r'$k_y \rho_i = $'+str(iky*par['kymin']),size=18)
        plt.loglog(kzn,'x-',label=r'$<k_z v_{t}>(v_{t}/L)$')
        ax=plt.axis()
        plt.axis((ax[0],ax[1],0,ax[3]))
        #plt.loglog(0.5*kzn,label=r'$0.5 <k_z v_{t}>(v_{t}/L)$')
        plt.xlabel('n',size=18)
        plt.ylabel(r'$\omega(v_{t}/L)$',size=18)
        plt.legend(loc='upper left')
        plt.show()

        if fit_spectra:
            p0=np.empty(2)
            p0[0]=1.0
            p0[1]=1.0
            startx=float(raw_input("Enter start value for fit:"))
            endx=float(raw_input("Enter end value for fit:"))
            fit=fit_function(herm_grid,kzn,p0,startx=startx,endx=endx,which_func='power')
            plt.figure(figsize=(4.5,3.0))
            fig=plt.gcf()
            fig.subplots_adjust(left=0.2)
            fig.subplots_adjust(bottom=0.22)
            plt.loglog(herm_grid,kzn,'x-',basex=10,basey=10)
            if startx==-1:
                sx0=0
            else:
                sx0=argmin(abs(herm_grid-startx))
            #   sx0=sx0-sx0/4
            if endx==-1:
                ex0=0
            else:
                ex0=argmin(abs(herm_grid-endx))+2
            fit_x=herm_grid[sx0:ex0]
            fit_func=1.1*fit[0]*fit_x**fit[1]
            plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
            plt.annotate(r'$\propto n$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
            plt.title('kzn')
            plt.xlabel(r'$n$',size=13)
            plt.show()
        else:
            plt.loglog(herm_grid,kzn,'x-',basex=10,basey=10)
            plt.loglog(herm_grid,0.5*kzphi*herm_grid**0.5,'--',basex=10,basey=10)
            plt.title('kzn')
            plt.xlabel(r'$n$',size=13)
            plt.show()

def get_hermite_dissipation_spectra(start_time=-1.0,end_time=-1.0,data_out=True,show_plots=False):
    """Plots time-averaged Hermite dissipation spectra from g_out.dat\n
    start_time: start time for time average (defaults to 0.0) \n
    end_time: end time for time average (defaults to end) \n
    data_out: flag for outputting data in an ascii file--useful \n
    for replotting or quick viewing later."""

    time=get_time_from_gout()
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    kxgrid_full=np.empty(par['nkx0']*2-1)
    for i in range(par['nkx0']):
        kxgrid_full[i+par['nkx0']-1]=kxgrid[i]
    for i in range(par['nkx0']-1):
        kxgrid_full[i]=-1.0*kxgrid[par['nkx0']-1-i]
    print kxgrid_full

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop
        
    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    C_sum=np.zeros(par['nv0'])
    HC_sum=np.zeros(par['nv0'])
    hypxyz_sum=np.zeros(par['nv0'])
    hypconv_sum=np.zeros(par['nv0'])
    ntime=iend-istart+1
    J0a=get_gyroavg()
    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        gt0=read_time_step_g(i)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        E_op=np.pi**0.5*np.conj(gt0)
        #Dissipation
        rhs=get_rhs_lin(gt0,1)
        etemp=E_op*rhs
        C_sum+=np.real(ksum_3d_n(etemp))
        rhs=get_rhs_lin(gt0,2)
        etemp=E_op*rhs
        HC_sum+=np.real(ksum_3d_n(etemp))
        #hyp_xyz
        rhs=get_rhs_lin(gt0,8)
        etemp=E_op*rhs
        hypxyz_sum+=np.real(ksum_3d_n(etemp))
        #hyp_conv
        rhs=get_rhs_lin(gt0,9)
        etemp=E_op*rhs
        hypconv_sum+=np.real(ksum_3d_n(etemp))

    C_sum=C_sum/float(ntime)
    HC_sum=HC_sum/float(ntime)
    hypxyz_sum=hypxyz_sum/float(ntime)
    hypconv_sum=hypconv_sum/float(ntime)

    print "C_sum",C_sum
    print "HC_sum",HC_sum
    print "hypxyzsum",hypxyz_sum
    print "hypconv_sum",hypconv_sum
    print "C_sum total",np.sum(C_sum)
    print "HC_sum total",np.sum(HC_sum)
    print "hypxyzsum total",np.sum(hypxyz_sum)
    print "hypconv_sum total",np.sum(hypconv_sum)

    if show_plots:
        plt.loglog(herm_grid,np.abs(C_sum),label='C',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(HC_sum),label='HC',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(hypxyz_sum),label='Hyp xyz',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(hypconv_sum),label='Hyp conv',basex=10,basey=10)
        plt.legend()
        plt.show()

        plt.plot(herm_grid,np.abs(C_sum),label='C')
        plt.plot(herm_grid,np.abs(HC_sum),label='HC')
        plt.plot(herm_grid,np.abs(hypxyz_sum),label='Hyp xyz')
        plt.plot(herm_grid,np.abs(hypconv_sum),label='Hyp conv')
        plt.legend()
        plt.show()

    file_name=par['diagdir'][1:-1]+'/hds_C.dat'
    np.savetxt(file_name,C_sum)
    file_name=par['diagdir'][1:-1]+'/hds_HC.dat'
    np.savetxt(file_name,HC_sum)
    file_name=par['diagdir'][1:-1]+'/hds_hypxyz.dat'
    np.savetxt(file_name,hypxyz_sum)
    file_name=par['diagdir'][1:-1]+'/hds_hypconv.dat'
    np.savetxt(file_name,hypconv_sum)


def get_energy_spectraF(start_time=-1.0,end_time=-1.0,data_out=False,which_n=-1,which_term=-1,fit_spectra=False,plot_7o3=False,which_fe='all'):
    """Plots time-averaged energy-related spectra from g_out.dat\n
    start_time: start time for time average (defaults to 0.0) \n
    end_time: end time for time average (defaults to end) \n
    which_fe: 'all'==> entropy and es, 'es'==> es, 'ent'==>entropy \n
    data_out: flag for outputting data in an ascii file--useful \n
    for replotting or quick viewing later."""

    if which_n > 0 and which_fe != 'ent':
        print "Warning! Calculating only entropy part!!"
        which_fe='ent'
        

    time=get_time_from_gout()
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    kxgrid_full=np.empty(par['nkx0']*2-1)
    for i in range(par['nkx0']):
        kxgrid_full[i+par['nkx0']-1]=kxgrid[i]
    for i in range(par['nkx0']-1):
        kxgrid_full[i]=-1.0*kxgrid[par['nkx0']-1-i]
    print kxgrid_full

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop
        
    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    diss_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    phisq_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    entk_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    hyp_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    entn_sum=np.zeros((par['nv0']))
    drive_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    ntime=iend-istart+1
    J0a=get_gyroavg()
    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        gt0=read_time_step_g(i)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        phi=get_phi(gt0)
        if which_n==-1 or which_n==0:
            phisq_sum+=np.real(np.conj(phi)*phi)
        E_op=np.pi**0.5*np.conj(gt0)
        for i in range(par['nkx0']):
            for j in range(par['nky0']):
                E_op[i,j,:,0]=E_op[i,j,:,0]+np.pi**0.25*J0a[i,j]*np.conj(phi[i,j,:])
        #Energy Drive
        rhs=get_rhs_lin(gt0,6)
        if which_n==-1:
            etemp=np.real(np.sum(E_op*rhs,axis=3))
        else:
            etemp=np.real((E_op*rhs)[:,:,:,which_n])
        drive_sum=drive_sum+etemp
        #Dissipation
        rhs=get_rhs_lin(gt0,1)
        if which_n==-1:
            etemp=np.real(np.sum(E_op*rhs,axis=3))
        else:
            etemp=np.real((E_op*rhs)[:,:,:,which_n])
        diss_sum=diss_sum+etemp
        #Entropy
        if which_n==-1:
            if which_fe=='all' or which_fe=='ent':
                entk_sum=entk_sum+get_entropy(gt0)
            if which_fe=='all' or which_fe=='es':
                entk_sum=entk_sum+get_fe_es(gt0)
    
        else:
            temp=get_entropy(gt0,sum_n=False)
            entk_sum=entk_sum+temp[:,:,:,which_n]

        if which_n==-1:
            entn_sum=entn_sum+get_entropy_hermite(gt0,kzind=-1)
        #Hyps
        rhs=get_rhs_lin(gt0,8)
        if which_n==-1:
            etemp=np.real(np.sum(E_op*rhs,axis=3))
        else:
            etemp=np.real((E_op*rhs)[:,:,:,which_n])
        hyp_sum=hyp_sum+etemp

    kx_out,ky_out,kz_out,herm_out=get_grids_shifted()

    if plot_7o3:
        spectkx_n7o3=np.arange(len(kx_out),dtype='float')
        spectky_n7o3=np.arange(len(ky_out),dtype='float')
        for i in range(len(kx_out)):
            spectkx_n7o3[i]=kx_out[i]**(-7.0/3.0)
        for i in range(len(ky_out)):
            spectky_n7o3[i]=ky_out[i]**(-7.0/3.0)
        kxp5_loc=np.argmin(np.abs(kx_out-0.6))
        print "kxp5_loc",kxp5_loc
        kyp5_loc=np.argmin(np.abs(ky_out-0.6))
        print "kyp5_loc",kyp5_loc
        kx1_loc=np.argmin(np.abs(kx_out-1.2))
        ky1_loc=np.argmin(np.abs(ky_out-1.2))
        #print "ky_out",ky_out



    if which_n==-1:
        entn_sum=entn_sum/float(ntime)
        plt.loglog(herm_out,entn_sum,basex=10,basey=10)
        plt.xlabel('Hermite n',size=18)
        plt.ylabel('Entropy',size=18)
        plt.show()
        plt.xlabel('Hermite n',size=18)
        plt.ylabel('Entropy',size=18)
        plt.semilogy(herm_out,entn_sum,basey=10)
        plt.show()

    if data_out and which_n==-1:
        file_name=par['diagdir'][1:-1]+'/entropy_hermite.dat'
        arr_out=np.empty((par['nv0'],2)) 
        arr_out[:,0]=herm_out
        arr_out[:,1]=entn_sum
        np.savetxt(file_name,arr_out)

    e2d_plot=np.empty((par['nkx0']*2-1,par['nky0']))


    if which_term== 3 or which_term==-1:
      diss=diss_sum/float(ntime)
      diss_kza=np.sum(diss,2)
      np.info(diss_kza)
      e2d_plot=plot_prep2d(diss_kza)

      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
      plt.title('Dissipation')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.ylabel(r'$k_y\rho_i$',size=13)
      time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
      plt.figtext(0.02,0.03,time_string)
      plt.colorbar()
      plt.show()


    if which_term==2 or which_term==-1:
      drive=drive_sum/float(ntime)
      drive_kza=np.sum(drive,2)
      np.info(drive_kza)
      e2d_plot=plot_prep2d(drive_kza)
      #drive_kza=np.roll(drive_kza,par['nky0']/2-1,axis=1)
      #plt.contour(ky_out,kx_out,drive_kza,200)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
      plt.title('Drive')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.ylabel(r'$k_y\rho_i$',size=13)
      plt.colorbar()
      plt.figtext(0.02,0.02,time_string)
      plt.show()

      drive_kykz=drive[0,:,:]
      drive_kykz=np.roll(drive_kykz,par['nky0']/2-1,axis=0)
      drive_kykz=np.roll(drive_kykz,par['nkz0']/2-1,axis=1)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kz_out,ky_out,drive_kykz,200)
      plt.title('Drive')
      plt.xlabel(r'$k_z R$')
      plt.ylabel(r'$k_y \rho_i$')
      plt.colorbar()
      plt.show()

    if which_term==1 or which_term==-1:
      entk_sum=np.real(entk_sum/float(ntime))
      entk_kza=np.sum(entk_sum,axis=2)
      e2d_plot=plot_prep2d(entk_kza)
      #entk_kza=np.roll(entk_kza,par['nky0']/2-1,axis=1)
      #plt.contour(ky_out,kx_out,entk_kza,200)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
      plt.title('Free Energy')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.ylabel(r'$k_y\rho_i$',size=13)
      time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
      plt.figtext(0.02,0.02,time_string)
      plt.colorbar()
      plt.show()

      entk_kykz=entk_sum[0,:,:]
      entk_kykz=np.roll(entk_kykz,par['nky0']/2-1,axis=0)
      entk_kykz=np.roll(entk_kykz,par['nkz0']/2-1,axis=1)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kz_out,ky_out,entk_kykz,200)
      plt.title('Free Energy')
      plt.xlabel(r'$k_z R$')
      plt.ylabel(r'$k_y \rho_i$')
      plt.colorbar()
      plt.show()

      if which_n==-1 or which_n==0:
          phisq_sum=phisq_sum/float(ntime)
          phisq_kza=np.sum(phisq_sum,axis=2)
          e2d_plot=plot_prep2d(phisq_kza)
          #entk_kza=np.roll(entk_kza,par['nky0']/2-1,axis=1)
          #plt.contour(ky_out,kx_out,entk_kza,200)
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.2)
          plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
          plt.title(r'$\phi^2$')
          plt.xlabel(r'$k_x\rho_i$',size=13)
          plt.ylabel(r'$k_y\rho_i$',size=13)
          time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
          plt.figtext(0.02,0.02,time_string)
          plt.colorbar()
          plt.show()

          phisq_kykz=phisq_sum[0,:,:]
          phisq_kykz=np.roll(phisq_kykz,par['nky0']/2-1,axis=0)
          phisq_kykz=np.roll(phisq_kykz,par['nkz0']/2-1,axis=1)
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.2)
          plt.contourf(kz_out,ky_out,phisq_kykz,200)
          plt.title(r'$\phi^2$')
          plt.xlabel(r'$k_z R$')
          plt.ylabel(r'$k_y \rho_i$')
          plt.colorbar()
          plt.show()





    if which_term==4 or which_term==-1:
      hyp_sum=hyp_sum/float(ntime)
      hyp_kza=np.sum(hyp_sum,axis=2)
      e2d_plot=plot_prep2d(hyp_kza)
      #entk_kza=np.roll(entk_kza,par['nky0']/2-1,axis=1)
      #plt.contour(ky_out,kx_out,entk_kza,200)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
      plt.title('Hyps')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.ylabel(r'$k_y\rho_i$',size=13)
      time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
      plt.figtext(0.02,0.02,time_string)
      plt.colorbar()
      plt.show()



    if which_term==3 or which_term==-1:
      diss_kx=np.sum(diss_kza,1)
      diss_kx[1:]=2.0*diss_kx[1:]
      diss_ky=np.sum(diss_kza[:,:],0)
      for i in range(1,par['nky0']/2):
          diss_ky[i]=diss_ky[i]+diss_ky[par['nky0']-i]
          diss_ky[par['nky0']-i]=diss_ky[i]
      diss_ky=np.roll(diss_ky,par['nky0']/2-1,axis=0)
      diss_kya=np.sum(diss,1)
      diss_kz=np.sum(diss_kya[:,:],0)
      for i in range(1,par['nkz0']/2):
          diss_kz[i]=diss_kz[i]+diss_kz[par['nkz0']-i]
          diss_kz[par['nkz0']-i]=diss_kz[i]
      diss_kz=np.roll(diss_kz,par['nkz0']/2-1,axis=0)

      #Linear plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kx_out,diss_kx,'x-')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.title('Dissipation')
      plt.xlim((0,par['kxmax0']))
      plt.show()



      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(ky_out,diss_ky,'x-')
      plt.title('Dissipation')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kymax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kz_out,diss_kz,'x-')
      plt.title('Dissipation')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kzmax0']))
      plt.show()


      if data_out:
          file_name=par['diagdir'][1:-1]+'/diss0_kx.dat'
          arr_out=np.empty((par['nkx0'],2)) 
          arr_out[:,0]=kx_out
          arr_out[:,1]=diss_kx
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/diss0_ky.dat'
          arr_out=np.empty((par['nky0'],2)) 
          arr_out[:,0]=ky_out
          arr_out[:,1]=diss_ky
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/diss0_kz.dat'
          arr_out=np.empty((par['nkz0'],2)) 
          arr_out[:,0]=kz_out
          arr_out[:,1]=diss_kz
          np.savetxt(file_name,arr_out)

      #log-log plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kx_out,-1.0*diss_kx,'x-',basex=10,basey=10)
      plt.title('Dissipation')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(ky_out,-1.0*diss_ky,'x-',basex=10,basey=10)
      plt.title('Dissipation')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kz_out,-1.0*diss_kz,'x-',basex=10,basey=10)
      plt.title('Dissipation')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()

    if which_term==2 or which_term==-1:
      drive_kx=np.sum(drive_kza,1)
      drive_kx[1:]=2.0*drive_kx[1:]
      drive_ky=np.sum(drive_kza[:,:],0)
      for i in range(1,par['nky0']/2):
          drive_ky[i]=drive_ky[i]+drive_ky[par['nky0']-i]
          drive_ky[par['nky0']-i]=drive_ky[i]
      drive_ky=np.roll(drive_ky,par['nky0']/2-1,axis=0)
      drive_kya=np.sum(drive,1)
      drive_kz=np.sum(drive_kya[:,:],0)
      for i in range(1,par['nkz0']/2):
          drive_kz[i]=drive_kz[i]+drive_kz[par['nkz0']-i]
          drive_kz[par['nkz0']-i]=drive_kz[i]
      drive_kz=np.roll(drive_kz,par['nkz0']/2-1,axis=0)

      if data_out:
          file_name=par['diagdir'][1:-1]+'/drive0_kx.dat'
          arr_out=np.empty((par['nkx0'],2)) 
          arr_out[:,0]=kx_out
          arr_out[:,1]=drive_kx
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/drive0_ky.dat'
          arr_out=np.empty((par['nky0'],2)) 
          arr_out[:,0]=ky_out
          arr_out[:,1]=drive_ky
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/drive0_kz.dat'
          arr_out=np.empty((par['nkz0'],2)) 
          arr_out[:,0]=kz_out
          arr_out[:,1]=drive_kz
          np.savetxt(file_name,arr_out)



    #drive_kx=np.sum(drive_kza,1)
    #drive_ky=drive_kza[0,:]+2.0*np.sum(drive_kza[1:,:],0)
    ##drive_ky=np.roll(drive_ky,par['nky0']/2-1,axis=0)
    #drive_kya=np.sum(drive,1)
    #drive_kz=drive_kya[0,:]+2.0*np.sum(drive_kya[1:,:],0)
    #drive_kz=np.roll(drive_kz,par['nkz0']/2-1,axis=0)

      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kx_out,drive_kx,'x-')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.title('Drive')
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kxmax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(ky_out,drive_ky,'x-')
      plt.title('Drive')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kymax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kz_out,drive_kz,'x-')
      plt.title('Drive')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kzmax0']))
      plt.show()


      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kx_out,drive_kx,'x-',basex=10,basey=10)
      plt.title('Drive')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(ky_out,drive_ky,'x-',basex=10,basey=10)
      plt.title('Drive')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kz_out,drive_kz,'x-',basex=10,basey=10)
      plt.title('Drive')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()


    if which_term==1 or which_term==-1:

      if which_n==-1 or which_n==0:
          phisq_kx=np.sum(phisq_kza,1)
          phisq_kx[1:]=2.0*phisq_kx[1:]
          phisq_ky=np.sum(phisq_kza[:,:],0)
          for i in range(1,par['nky0']/2):
              phisq_ky[i]=phisq_ky[i]+phisq_ky[par['nky0']-i]
              phisq_ky[par['nky0']-i]=phisq_ky[i]
          phisq_ky=np.roll(phisq_ky,par['nky0']/2-1,axis=0)
          phisq_kya=np.sum(phisq_sum,1)
          phisq_kz=np.sum(phisq_kya[:,:],0)
          for i in range(1,par['nkz0']/2):
              phisq_kz[i]=phisq_kz[i]+phisq_kz[par['nkz0']-i]
              phisq_kz[par['nkz0']-i]=phisq_kz[i]
          phisq_kz=np.roll(phisq_kz,par['nkz0']/2-1,axis=0)


          #log-log plot for phisq
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.22)
          if fit_spectra:
              print "Identify fit range from plot."
          plt.loglog(kx_out,phisq_kx,'x-',basex=10,basey=10)
          plt.title(r'$\phi^2$')
          plt.xlabel(r'$k_x\rho_i$',size=13)
          plt.figtext(0.02,0.02,time_string)
          plt.show()
          if fit_spectra:
              startx=float(raw_input("Enter start value for fit:"))
              endx=float(raw_input("Enter end value for fit:"))
    
              fit=fit_function(kx_out,phisq_kx,(1.0,-1.0),startx=startx,endx=endx,which_func='power')
              plt.figure(figsize=(4.5,3.0))
              fig=plt.gcf()
              fig.subplots_adjust(left=0.2)
              fig.subplots_adjust(bottom=0.22)
              plt.loglog(kx_out,phisq_kx,'x-',basex=10,basey=10)
              if startx==-1:
                  sx0=0
              else:
                  sx0=argmin(abs(kx_out-startx))
             #     sx0=sx0-sx0/4
              if endx==-1:
                  ex0=0
              else:
                  ex0=argmin(abs(kx_out-endx))
              fit_x=kx_out[sx0:ex0]
              fit_func=2.0*fit[0]*fit_x**fit[1]
              plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
              plt.annotate(r'$\propto k$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
              plt.title(r'$\phi^2$')
              plt.xlabel(r'$k_x\rho_i$',size=13)
              plt.figtext(0.02,0.02,time_string)
              plt.show()
               
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.22)
          if fit_spectra:
              print "Identify fit range from plot."
          plt.loglog(ky_out,phisq_ky,'x-',basex=10,basey=10)
          plt.title(r'$\phi^2$')
          plt.xlabel(r'$k_y\rho_i$',size=13)
          plt.figtext(0.02,0.02,time_string)
          plt.show()
          if fit_spectra:
              startx=float(raw_input("Enter start value for fit:"))
              endx=float(raw_input("Enter end value for fit:"))
              fit=fit_function(ky_out,phisq_ky,(1.0,-1.0),startx=startx,endx=endx,which_func='power')
              plt.figure(figsize=(4.5,3.0))
              fig=plt.gcf()
              fig.subplots_adjust(left=0.2)
              fig.subplots_adjust(bottom=0.22)
              plt.loglog(ky_out,phisq_ky,'x-',basex=10,basey=10)
              if startx==-1:
                  sx0=0
              else:
                  sx0=argmin(abs(ky_out-startx))
              if endx==-1:
                  ex0=0
              else:
                  ex0=argmin(abs(ky_out-endx))
              #print "sx0,ex0",sx0,ex0
              #print len(ky_out)
              #print "ky_out",ky_out
              fit_x=ky_out[sx0:ex0]
              fit_func=2.0*fit[0]*fit_x**fit[1]
              #print fit_x
              #print fit_func
              plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
              plt.annotate(r'$\propto k$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
              plt.title(r'$\phi^2$')
              plt.xlabel(r'$k_y\rho_i$',size=13)
              plt.figtext(0.02,0.02,time_string)
              plt.show()

          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.22)
          if fit_spectra:
              print "Identify fit range from plot."
          plt.loglog(kz_out,phisq_kz,'x-',basex=10,basey=10)
          plt.title(r'$\phi^2$')
          plt.xlabel(r'$k_z R$',size=13)
          plt.figtext(0.02,0.02,time_string)
          plt.show()
          if fit_spectra:
              startx=float(raw_input("Enter start value for fit:"))
              endx=float(raw_input("Enter end value for fit:"))
              fit=fit_function(kz_out,phisq_kz,(1.0,-1.0),startx=startx,endx=endx,which_func='power')
              plt.figure(figsize=(4.5,3.0))
              fig=plt.gcf()
              fig.subplots_adjust(left=0.2)
              fig.subplots_adjust(bottom=0.22)
              plt.loglog(kz_out,phisq_kz,'x-',basex=10,basey=10)
              if startx==-1:
                  sx0=0
              else:
                  sx0=argmin(abs(kz_out-startx))
              if endx==-1:
                  ex0=0
              else:
                  ex0=argmin(abs(kz_out-endx))
              #print "sx0,ex0",sx0,ex0
              #print len(kz_out)
              #print "kz_out",kz_out
              fit_x=kz_out[sx0:ex0]
              fit_func=2.0*fit[0]*fit_x**fit[1]
              #print fit_x
              #print fit_func
              plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
              plt.annotate(r'$\propto k$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
              plt.title(r'$\phi^2$')
              plt.xlabel(r'$k_z R$',size=13)
              plt.figtext(0.02,0.02,time_string)
              plt.show()
          if data_out:
              file_name=par['diagdir'][1:-1]+'/phisq_kx.dat'
              arr_out=np.empty((par['nkx0'],2)) 
              arr_out[:,0]=kx_out
              arr_out[:,1]=phisq_kx
              np.savetxt(file_name,arr_out)
              file_name=par['diagdir'][1:-1]+'/phisq_ky.dat'
              arr_out=np.empty((par['nky0'],2)) 
              arr_out[:,0]=ky_out
              arr_out[:,1]=phisq_ky
              np.savetxt(file_name,arr_out)
              file_name=par['diagdir'][1:-1]+'/phisq_kz.dat'
              arr_out=np.empty((par['nkz0'],2)) 
              arr_out[:,0]=kz_out
              arr_out[:,1]=phisq_kz
              np.savetxt(file_name,arr_out)




      entk_kx=np.sum(entk_kza,1)
      print "type entk_kx", type(entk_kx[0])
      entk_kx[1:]=2.0*entk_kx[1:]
      entk_ky=np.sum(entk_kza[:,:],0)
      for i in range(1,par['nky0']/2):
          entk_ky[i]=entk_ky[i]+entk_ky[par['nky0']-i]
          entk_ky[par['nky0']-i]=entk_ky[i]
      entk_ky=np.roll(entk_ky,par['nky0']/2-1,axis=0)
      entk_kya=np.sum(entk_sum,1)
      entk_kz=np.sum(entk_kya[:,:],0)
      for i in range(1,par['nkz0']/2):
          entk_kz[i]=entk_kz[i]+entk_kz[par['nkz0']-i]
          entk_kz[par['nkz0']-i]=entk_kz[i]
      entk_kz=np.roll(entk_kz,par['nkz0']/2-1,axis=0)

      if data_out:
          file_name=par['diagdir'][1:-1]+'/entk0_kx.dat'
          arr_out=np.empty((par['nkx0'],2)) 
          arr_out[:,0]=kx_out
          arr_out[:,1]=entk_kx
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/entk0_ky.dat'
          arr_out=np.empty((par['nky0'],2)) 
          arr_out[:,0]=ky_out
          arr_out[:,1]=entk_ky
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/entk0_kz.dat'
          arr_out=np.empty((par['nkz0'],2)) 
          arr_out[:,0]=kz_out
          arr_out[:,1]=entk_kz
          np.savetxt(file_name,arr_out)



      #entk_kx=np.sum(entk_kza,1)
      #entk_ky=entk_kza[0,:]+2.0*np.sum(entk_kza[1:,:],0)
      ##diss_ky=np.roll(diss_ky,par['nky0']/2-1,axis=0)
      #entk_kya=np.sum(entk_sum,1)
      #entk_kz=entk_kya[0,:]+2.0*np.sum(entk_kya[1:,:],0)
      #entk_kz=np.roll(entk_kz,par['nkz0']/2-1,axis=0)

      #Linear plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kx_out,entk_kx,'x-')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.title('Free Energy')
      plt.xlim((0,par['kxmax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(ky_out,entk_ky,'x-')
      plt.title('Free Energy')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kymax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kz_out,entk_kz,'x-')
      plt.title('Free Energy')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kzmax0']))
      plt.show()

      #log-log plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      if fit_spectra:
          print "Identify fit range from plot."
      plt.loglog(kx_out,entk_kx,'x-',basex=10,basey=10)
      if plot_7o3:
          spectkx_n7o3/=0.45*spectkx_n7o3[kxp5_loc]/entk_kx[kxp5_loc]
          plt.loglog(kx_out[kxp5_loc:kx1_loc],spectkx_n7o3[kxp5_loc:kx1_loc],':',color='black',basex=10,basey=10,label=r'$\propto k^{-7/3}$')
          plt.legend(loc='lower left')
      plt.title('Free Energy')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      if fit_spectra:
          startx=float(raw_input("Enter start value for fit:"))
          endx=float(raw_input("Enter end value for fit:"))

          fit=fit_function(kx_out,entk_kx,(1.0,-1.0),startx=startx,endx=endx,which_func='power')
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.22)
          plt.loglog(kx_out,entk_kx,'x-',basex=10,basey=10)
          if startx==-1:
              sx0=0
          else:
              sx0=argmin(abs(kx_out-startx))
         #     sx0=sx0-sx0/4
          if endx==-1:
              ex0=0
          else:
              ex0=argmin(abs(kx_out-endx))
          fit_x=kx_out[sx0:ex0]
          fit_func=2.0*fit[0]*fit_x**fit[1]
          plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
          plt.annotate(r'$\propto k$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
          plt.title('Free Energy')
          plt.xlabel(r'$k_x\rho_i$',size=13)
          plt.figtext(0.02,0.02,time_string)
          plt.show()
           
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      if fit_spectra:
          print "Identify fit range from plot."
      plt.loglog(ky_out,entk_ky,'x-',basex=10,basey=10)
      if plot_7o3:
          spectky_n7o3/=0.45*spectky_n7o3[kyp5_loc]/entk_ky[kyp5_loc]
          plt.loglog(ky_out[kyp5_loc:ky1_loc],spectky_n7o3[kyp5_loc:ky1_loc],':',color='black',basex=10,basey=10,label=r'$\propto k^{-7/3}$')
          plt.legend(loc='lower left')
      plt.title('Free Energy')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      if fit_spectra:
          startx=float(raw_input("Enter start value for fit:"))
          endx=float(raw_input("Enter end value for fit:"))
          fit=fit_function(ky_out,entk_ky,(1.0,-1.0),startx=startx,endx=endx,which_func='power')
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.22)
          plt.loglog(ky_out,entk_ky,'x-',basex=10,basey=10)
          if startx==-1:
              sx0=0
          else:
              sx0=argmin(abs(ky_out-startx))
          if endx==-1:
              ex0=0
          else:
              ex0=argmin(abs(ky_out-endx))
          #print "sx0,ex0",sx0,ex0
          #print len(ky_out)
          #print "ky_out",ky_out
          fit_x=ky_out[sx0:ex0]
          fit_func=2.0*fit[0]*fit_x**fit[1]
          #print fit_x
          #print fit_func
          plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
          plt.annotate(r'$\propto k$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
          plt.title('Free Energy')
          plt.xlabel(r'$k_y\rho_i$',size=13)
          plt.figtext(0.02,0.02,time_string)
          plt.show()

      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      if fit_spectra:
          print "Identify fit range from plot."
      plt.loglog(kz_out,entk_kz,'x-',basex=10,basey=10)
      plt.title('Free Energy')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      if fit_spectra:
          startx=float(raw_input("Enter start value for fit:"))
          endx=float(raw_input("Enter end value for fit:"))
          fit=fit_function(kz_out,entk_kz,(1.0,-1.0),startx=startx,endx=endx,which_func='power')
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.22)
          plt.loglog(kz_out,entk_kz,'x-',basex=10,basey=10)
          if startx==-1:
              sx0=0
          else:
              sx0=argmin(abs(kz_out-startx))
          if endx==-1:
              ex0=0
          else:
              ex0=argmin(abs(kz_out-endx))
          #print "sx0,ex0",sx0,ex0
          #print len(kz_out)
          #print "kz_out",kz_out
          fit_x=kz_out[sx0:ex0]
          fit_func=2.0*fit[0]*fit_x**fit[1]
          #print fit_x
          #print fit_func
          plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
          plt.annotate(r'$\propto k$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
          plt.title('Free Energy')
          plt.xlabel(r'$k_y\rho_i$',size=13)
          plt.figtext(0.02,0.02,time_string)
          plt.show()


    if which_term==4 or which_term==-1:
      hyp_kx=np.sum(hyp_kza,1)
      hyp_kx[1:]=2.0*hyp_kx[1:]
      hyp_ky=np.sum(hyp_kza[:,:],0)
      for i in range(1,par['nky0']/2):
          hyp_ky[i]=hyp_ky[i]+hyp_ky[par['nky0']-i]
          hyp_ky[par['nky0']-i]=hyp_ky[i]
      hyp_ky=np.roll(hyp_ky,par['nky0']/2-1,axis=0)
      hyp_kya=np.sum(hyp_sum,1)
      hyp_kz=np.sum(hyp_kya[:,:],0)
      for i in range(1,par['nkz0']/2):
          hyp_kz[i]=hyp_kz[i]+hyp_kz[par['nkz0']-i]
          hyp_kz[par['nkz0']-i]=hyp_kz[i]
      hyp_kz=np.roll(hyp_kz,par['nkz0']/2-1,axis=0)
  
      if data_out:
          file_name=par['diagdir'][1:-1]+'/hyp0_kx.dat'
          arr_out=np.empty((par['nkx0'],2)) 
          arr_out[:,0]=kx_out
          arr_out[:,1]=hyp_kx
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/hyp0_ky.dat'
          arr_out=np.empty((par['nky0'],2)) 
          arr_out[:,0]=ky_out
          arr_out[:,1]=hyp_ky
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/hyp0_kz.dat'
          arr_out=np.empty((par['nkz0'],2)) 
          arr_out[:,0]=kz_out
          arr_out[:,1]=hyp_kz
          np.savetxt(file_name,arr_out)

      #Linear plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kx_out,hyp_kx,'x-')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.title('Hyps')
      plt.xlim((0,par['kxmax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(ky_out,hyp_ky,'x-')
      plt.title('Hyps')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kymax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kz_out,hyp_kz,'x-')
      plt.title('Hyps')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kzmax0']))
      plt.show()

      #log-log plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kx_out,-hyp_kx,'x-',basex=10,basey=10)
      plt.title('Hyps')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(ky_out,-hyp_ky,'x-',basex=10,basey=10)
      plt.title('Hyps')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kz_out,-hyp_kz,'x-',basex=10,basey=10)
      plt.title('Hyps')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()


    #diss=diss_sum/float(ntime)
    #diss_kza=np.sum(diss,2)
    #np.info(diss_kza)
    #e2d_plot=plot_prep2d(diss_kza)
    ##plt.contour(ky_out,kxgrid_full,e2d_plot,200)
    #plt.figure(figsize=(4.5,3.0))
    #fig=plt.gcf()
    #fig.subplots_adjust(left=0.2)
    #fig.subplots_adjust(bottom=0.2)
    #plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
    #plt.title('Dissipation')
    #plt.xlabel(r'$k_x\rho_i$',size=13)
    #plt.ylabel(r'$k_y\rho_i$',size=13)
    #time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
    #plt.figtext(0.02,0.02,time_string)
    #plt.colorbar()
    #plt.show()

    #drive=drive_sum/float(ntime)
    #drive_kza=np.sum(drive,2)
    #np.info(drive_kza)
    #e2d_plot=plot_prep2d(drive_kza)
    ##drive_kza=np.roll(drive_kza,par['nky0']/2-1,axis=1)
    ##plt.contour(ky_out,kx_out,drive_kza,200)
    #plt.figure(figsize=(4.5,3.0))
    #fig=plt.gcf()
    #fig.subplots_adjust(left=0.2)
    #fig.subplots_adjust(bottom=0.2)
    #plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
    #plt.title('Drive')
    #plt.xlabel(r'$k_x\rho_i$',size=13)
    #plt.ylabel(r'$k_y\rho_i$',size=13)
    #plt.colorbar()
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()

    #entk_sum=entk_sum/float(ntime)
    #entk_kza=np.sum(entk_sum,axis=2)
    #e2d_plot=plot_prep2d(entk_kza)
    ##entk_kza=np.roll(entk_kza,par['nky0']/2-1,axis=1)
    ##plt.contour(ky_out,kx_out,entk_kza,200)
    #plt.figure(figsize=(4.5,3.0))
    #fig=plt.gcf()
    #fig.subplots_adjust(left=0.2)
    #fig.subplots_adjust(bottom=0.2)
    #plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
    #plt.title('Free Energy')
    #plt.xlabel(r'$k_x\rho_i$',size=13)
    #plt.ylabel(r'$k_y\rho_i$',size=13)
    #time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
    #plt.figtext(0.02,0.02,time_string)
    #plt.colorbar()
    #plt.show()

    #hyp_sum=hyp_sum/float(ntime)
    #hyp_kza=np.sum(hyp_sum,axis=2)
    #e2d_plot=plot_prep2d(hyp_kza)
    #plt.figure(figsize=(4.5,3.0))
    #fig=plt.gcf()
    #fig.subplots_adjust(left=0.2)
    #fig.subplots_adjust(bottom=0.2)
    #plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
    #plt.title('Hyps')
    #plt.xlabel(r'$k_x\rho_i$',size=13)
    #plt.ylabel(r'$k_y\rho_i$',size=13)
    #time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
    #plt.figtext(0.02,0.02,time_string)
    #plt.colorbar()
    #plt.show()

    #diss_kx=np.sum(diss_kza,1)
    #diss_kx[1:]=2.0*diss_kx[1:]
    #diss_ky=np.sum(diss_kza[:,:],0)
    #for i in range(1,par['nky0']/2):
    #    diss_ky[i]=diss_ky[i]+diss_ky[par['nky0']-i]
    #    diss_ky[par['nky0']-i]=diss_ky[i]
    #diss_ky=np.roll(diss_ky,par['nky0']/2-1,axis=0)
    #diss_kya=np.sum(diss,1)
    #diss_kz=np.sum(diss_kya[:,:],0)
    #for i in range(1,par['nkz0']/2):
    #    diss_kz[i]=diss_kz[i]+diss_kz[par['nkz0']-i]
    #    diss_kz[par['nkz0']-i]=diss_kz[i]
    #diss_kz=np.roll(diss_kz,par['nkz0']/2-1,axis=0)

    ##Linear plot
    #plt.plot(kx_out,diss_kx)
    #plt.xlabel(r'$k_x\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.title('Dissipation')
    #plt.xlim((0,par['kxmax0']))
    #plt.show()
    #plt.plot(ky_out,diss_ky)
    #plt.title('Dissipation')
    #plt.xlabel(r'$k_y\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.xlim((0,par['kymax0']))
    #plt.show()
    #plt.plot(kz_out,diss_kz)
    #plt.title('Dissipation')
    #plt.xlabel(r'$k_z\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.xlim((0,par['kzmax0']))
    #plt.show()

    #if data_out:
    #    file_name=par['diagdir'][1:-1]+'/diss_kx.dat'
    #    arr_out=np.empty((par['nkx0'],2)) 
    #    arr_out[:,0]=kx_out
    #    arr_out[:,1]=diss_kx
    #    np.savetxt(file_name,arr_out)
    #    file_name=par['diagdir'][1:-1]+'/diss_ky.dat'
    #    arr_out=np.empty((par['nky0'],2)) 
    #    arr_out[:,0]=ky_out
    #    arr_out[:,1]=diss_ky
    #    np.savetxt(file_name,arr_out)
    #    file_name=par['diagdir'][1:-1]+'/diss_kz.dat'
    #    arr_out=np.empty((par['nkz0'],2)) 
    #    arr_out[:,0]=kz_out
    #    arr_out[:,1]=diss_kz
    #    np.savetxt(file_name,arr_out)

    ##log-log plot
    #plt.loglog(kx_out,-1.0*diss_kx,basex=10,basey=10)
    #plt.title('Dissipation')
    #plt.xlabel(r'$k_x\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()
    #plt.loglog(ky_out,-1.0*diss_ky,basex=10,basey=10)
    #plt.title('Dissipation')
    #plt.xlabel(r'$k_y\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()
    #plt.loglog(kz_out,-1.0*diss_kz,basex=10,basey=10)
    #plt.title('Dissipation')
    #plt.xlabel(r'$k_z\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()

    #drive_kx=np.sum(drive_kza,1)
    #drive_kx[1:]=2.0*drive_kx[1:]
    #drive_ky=np.sum(drive_kza[:,:],0)
    #for i in range(1,par['nky0']/2):
    #    drive_ky[i]=drive_ky[i]+drive_ky[par['nky0']-i]
    #    drive_ky[par['nky0']-i]=drive_ky[i]
    #drive_ky=np.roll(drive_ky,par['nky0']/2-1,axis=0)
    #drive_kya=np.sum(drive,1)
    #drive_kz=np.sum(drive_kya[:,:],0)
    #for i in range(1,par['nkz0']/2):
    #    drive_kz[i]=drive_kz[i]+drive_kz[par['nkz0']-i]
    #    drive_kz[par['nkz0']-i]=drive_kz[i]
    #drive_kz=np.roll(drive_kz,par['nkz0']/2-1,axis=0)

    #if data_out:
    #    file_name=par['diagdir'][1:-1]+'/drive_kx.dat'
    #    arr_out=np.empty((par['nkx0'],2)) 
    #    arr_out[:,0]=kx_out
    #    arr_out[:,1]=drive_kx
    #    np.savetxt(file_name,arr_out)
    #    file_name=par['diagdir'][1:-1]+'/drive_ky.dat'
    #    arr_out=np.empty((par['nky0'],2)) 
    #    arr_out[:,0]=ky_out
    #    arr_out[:,1]=drive_ky
    #    np.savetxt(file_name,arr_out)
    #    file_name=par['diagdir'][1:-1]+'/drive_kz.dat'
    #    arr_out=np.empty((par['nkz0'],2)) 
    #    arr_out[:,0]=kz_out
    #    arr_out[:,1]=drive_kz
    #    np.savetxt(file_name,arr_out)

    #plt.plot(kx_out,drive_kx)
    #plt.xlabel(r'$k_x\rho_i$',size=18)
    #plt.title('Drive')
    #plt.figtext(0.02,0.02,time_string)
    #plt.xlim((0,par['kxmax0']))
    #plt.show()
    #plt.plot(ky_out,drive_ky)
    #plt.title('Drive')
    #plt.xlabel(r'$k_y\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.xlim((0,par['kymax0']))
    #plt.show()
    #plt.plot(kz_out,drive_kz)
    #plt.title('Drive')
    #plt.xlabel(r'$k_z\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.xlim((0,par['kzmax0']))
    #plt.show()


    #plt.loglog(kx_out,drive_kx,basex=10,basey=10)
    #plt.title('Drive')
    #plt.xlabel(r'$k_x\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()
    #plt.loglog(ky_out,drive_ky,basex=10,basey=10)
    #plt.title('Drive')
    #plt.xlabel(r'$k_y\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()
    #plt.loglog(kz_out,drive_kz,basex=10,basey=10)
    #plt.title('Drive')
    #plt.xlabel(r'$k_z\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()


    #entk_kx=np.sum(entk_kza,1)
    #entk_kx[1:]=2.0*entk_kx[1:]
    #entk_ky=np.sum(entk_kza[:,:],0)
    #for i in range(1,par['nky0']/2):
    #    entk_ky[i]=entk_ky[i]+entk_ky[par['nky0']-i]
    #    entk_ky[par['nky0']-i]=entk_ky[i]
    #entk_ky=np.roll(entk_ky,par['nky0']/2-1,axis=0)
    #entk_kya=np.sum(entk_sum,1)
    #entk_kz=np.sum(entk_kya[:,:],0)
    #for i in range(1,par['nkz0']/2):
    #    entk_kz[i]=entk_kz[i]+entk_kz[par['nkz0']-i]
    #    entk_kz[par['nkz0']-i]=entk_kz[i]
    #entk_kz=np.roll(entk_kz,par['nkz0']/2-1,axis=0)

    #if data_out:
    #    file_name=par['diagdir'][1:-1]+'/entk_kx.dat'
    #    arr_out=np.empty((par['nkx0'],2)) 
    #    arr_out[:,0]=kx_out
    #    arr_out[:,1]=entk_kx
    #    np.savetxt(file_name,arr_out)
    #    file_name=par['diagdir'][1:-1]+'/entk_ky.dat'
    #    arr_out=np.empty((par['nky0'],2)) 
    #    arr_out[:,0]=ky_out
    #    arr_out[:,1]=entk_ky
    #    np.savetxt(file_name,arr_out)
    #    file_name=par['diagdir'][1:-1]+'/entk_kz.dat'
    #    arr_out=np.empty((par['nkz0'],2)) 
    #    arr_out[:,0]=kz_out
    #    arr_out[:,1]=entk_kz
    #    np.savetxt(file_name,arr_out)

    ##Linear plot
    #plt.plot(kx_out,entk_kx)
    #plt.xlabel(r'$k_x\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.title('Entropy')
    #plt.xlim((0,par['kxmax0']))
    #plt.show()
    #plt.plot(ky_out,entk_ky)
    #plt.title('Entropy')
    #plt.xlabel(r'$k_y\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.xlim((0,par['kymax0']))
    #plt.show()
    #plt.plot(kz_out,entk_kz)
    #plt.title('Entropy')
    #plt.xlabel(r'$k_z\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.xlim((0,par['kzmax0']))
    #plt.show()

    ##log-log plot
    #plt.loglog(kx_out,entk_kx,basex=10,basey=10)
    #plt.title('Entropy')
    #plt.xlabel(r'$k_x\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()
    #plt.loglog(ky_out,entk_ky,basex=10,basey=10)
    #plt.title('Entropy')
    #plt.xlabel(r'$k_y\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()
    #plt.loglog(kz_out,entk_kz,basex=10,basey=10)
    #plt.title('Entropy')
    #plt.xlabel(r'$k_z\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()

    #hyp_kx=np.sum(hyp_kza,1)
    #hyp_kx[1:]=2.0*hyp_kx[1:]
    #hyp_ky=np.sum(hyp_kza[:,:],0)
    #for i in range(1,par['nky0']/2):
    #    hyp_ky[i]=hyp_ky[i]+hyp_ky[par['nky0']-i]
    #    hyp_ky[par['nky0']-i]=hyp_ky[i]
    #hyp_ky=np.roll(hyp_ky,par['nky0']/2-1,axis=0)
    #hyp_kya=np.sum(hyp_sum,1)
    #hyp_kz=np.sum(hyp_kya[:,:],0)
    #for i in range(1,par['nkz0']/2):
    #    hyp_kz[i]=hyp_kz[i]+hyp_kz[par['nkz0']-i]
    #    hyp_kz[par['nkz0']-i]=hyp_kz[i]
    #hyp_kz=np.roll(hyp_kz,par['nkz0']/2-1,axis=0)

    #if data_out:
    #    file_name=par['diagdir'][1:-1]+'/hyp_kx.dat'
    #    arr_out=np.empty((par['nkx0'],2)) 
    #    arr_out[:,0]=kx_out
    #    arr_out[:,1]=hyp_kx
    #    np.savetxt(file_name,arr_out)
    #    file_name=par['diagdir'][1:-1]+'/hyp_ky.dat'
    #    arr_out=np.empty((par['nky0'],2)) 
    #    arr_out[:,0]=ky_out
    #    arr_out[:,1]=hyp_ky
    #    np.savetxt(file_name,arr_out)
    #    file_name=par['diagdir'][1:-1]+'/hyp_kz.dat'
    #    arr_out=np.empty((par['nkz0'],2)) 
    #    arr_out[:,0]=kz_out
    #    arr_out[:,1]=hyp_kz
    #    np.savetxt(file_name,arr_out)

    ##Linear plot
    #plt.plot(kx_out,hyp_kx)
    #plt.xlabel(r'$k_x\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.title('Hyps')
    #plt.xlim((0,par['kxmax0']))
    #plt.show()
    #plt.plot(ky_out,hyp_ky)
    #plt.title('Hyps')
    #plt.xlabel(r'$k_y\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.xlim((0,par['kymax0']))
    #plt.show()
    #plt.plot(kz_out,hyp_kz)
    #plt.title('Hyps')
    #plt.xlabel(r'$k_z\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.xlim((0,par['kzmax0']))
    #plt.show()

    ##log-log plot
    #plt.loglog(kx_out,-hyp_kx,basex=10,basey=10)
    #plt.title('Hyps')
    #plt.xlabel(r'$k_x\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()
    #plt.loglog(ky_out,-hyp_ky,basex=10,basey=10)
    #plt.title('Hyps')
    #plt.xlabel(r'$k_y\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()
    #plt.loglog(kz_out,-hyp_kz,basex=10,basey=10)
    #plt.title('Hyps')
    #plt.xlabel(r'$k_z\rho_i$',size=18)
    #plt.figtext(0.02,0.02,time_string)
    #plt.show()

def get_energy_spectra_nk(start_time=-1.0,end_time=-1.0,data_out=False,delta_kperp = 0.1,which_term=2):
    """Plots time-averaged n,k_perp resolved spectra \n
    start_time: start time for time average (defaults to 0.0) \n
    end_time: end time for time average (defaults to end) \n
    data_out: flag for outputting data in an ascii file--useful \n
    delta_kperp: width of kperp bins \n
    which_term: 1=free energy, 2=coll+hcoll, 3=hyp_kperp \n
    for replotting or quick viewing later."""

    time=get_time_from_gout()
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    kxgrid_full=np.empty(par['nkx0']*2-1)
    for i in range(par['nkx0']):
        kxgrid_full[i+par['nkx0']-1]=kxgrid[i]
    for i in range(par['nkx0']-1):
        kxgrid_full[i]=-1.0*kxgrid[par['nkx0']-1-i]
    print kxgrid_full

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop
        
    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    eterm_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0'],par['nv0']))
    #hypk_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0'],par['nv0']))
    ntime=iend-istart+1
    J0a=get_gyroavg()

    if which_term==1:
        #free energy
        print "Calculating free energy spectra."
        file_base='energy'
    elif which_term==2:
        #collisions and hypercollisions
        print "Calculating collisional and hypercollisional dissipation."
        file_base='coll'
    elif which_term==3:
        #kperp hyperdiffusion
        print "Calculating kperp and kz hyperdiffusion."
        file_base='hypk'
    else:
        print "Error!"
        stop

    kpmax=np.sqrt(par['kxmax0']**2+par['kymax0']**2)
    print "kpmax",kpmax
    nkp=ceil(kpmax/delta_kperp)
    kpmax_top=nkp*delta_kperp
    print "nkp",nkp
    kpgrid=np.arange(nkp)/float(nkp-1)*(kpmax_top-delta_kperp)+delta_kperp
    print "kpgrid",kpgrid
    print "kzgrid",kzgrid
    

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        gt0=read_time_step_g(i)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        phi=get_phi(gt0)
        E_op=np.pi**0.5*np.conj(gt0)
        for i in range(par['nkx0']):
            for j in range(par['nky0']):
                E_op[i,j,:,0]=E_op[i,j,:,0]+np.pi**0.25*J0a[i,j]*np.conj(phi[i,j,:])
        #Dissipation
        if which_term==1:
            #free energy
            eterm_sum+=np.real((E_op*gt0)[:,:,:,:])
        elif which_term==2:
            #collisions and hypercollisions
            rhs=get_rhs_lin(gt0,1)
            rhs+=get_rhs_lin(gt0,2)
            eterm_sum+=np.real((E_op*rhs)[:,:,:,:])
        elif which_term==3:
            #kperp hyperdiffusion
            rhs=get_rhs_lin(gt0,8)
            eterm_sum+=np.real((E_op*rhs)[:,:,:,:])
        else:
            print "Error!"
            stop
        #Entropy
        #if which_n==-1:
        #    if which_fe=='all' or which_fe=='ent':
        #        entk_sum=entk_sum+get_entropy(gt0)
        #    if which_fe=='all' or which_fe=='es':
        #        entk_sum=entk_sum+get_fe_es(gt0)
        #else:
        #    temp=get_entropy(gt0,sum_n=False)
        #    entk_sum=entk_sum+temp[:,:,:,which_n]

        #if which_n==-1:
        #    entn_sum=entn_sum+get_entropy_hermite(gt0,kzind=-1)
        #Hyps
        #rhs=get_rhs_lin(gt0,8)
        #if which_n==-1:
        #    etemp=np.real(np.sum(E_op*rhs,axis=3))
        #else:
        #    etemp=np.real((E_op*rhs)[:,:,:,which_n])
        #hyp_sum=hyp_sum+etemp

    eterm_sum=eterm_sum/float(ntime)
    eterm_kpkzn=np.zeros((nkp,par['nkz0'],par['nv0']),dtype='float')
    eterm_kpn=np.zeros((nkp,par['nv0']),dtype='float')

    for i in range(par['nkx0']):
        for j in range(par['nky0']):
          if j != par['nky0']/2:
            kperp=np.sqrt(kxgrid[i]**2+kygrid[j]**2)    
            kpind=np.argmin(np.abs(kpgrid-kperp))
            if kpgrid[kpind] < kperp:
                kpind+=1
            #print "kperp,kpind",kperp,kpind
            #print "kpgrid[kpind]",kpgrid[kpind]
            eterm_kpkzn[kpind,:,:]+=eterm_sum[i,j,:,:]

    eterm_kpn=np.sum(eterm_kpkzn,axis=1)
    eterm_kzn=np.sum(eterm_kpkzn,axis=0)
    print np.shape(eterm_kzn)

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    fig.subplots_adjust(left=0.2)
    fig.subplots_adjust(bottom=0.2)
    plt.contourf(kpgrid,herm_grid,np.transpose(eterm_kpn),100)
    plt.title(file_base)
    plt.xlabel(r'$k_{\perp}\rho_i$',size=13)
    plt.ylabel(r'Hermite n',size=13)
    time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
    #plt.figtext(0.02,0.03,time_string)
    plt.colorbar()
    plt.show()

    #eterm_kzn_temp=np.roll(eterm_kzn,par['nkz0']/2-1,axis=0)
    #kzgrid_temp=np.roll(kzgrid,par['nkz0']/2-1)
    eterm_kzn_out=np.zeros((par['nkz0']/2-1,par['nv0']),dtype='float')
    eterm_kzn_out[0,:]=eterm_kzn[0,:]
    for i in range(1,par['nkz0']/2-1):
        eterm_kzn_out[i,:]=eterm_kzn[i,:]+eterm_kzn[par['nkz0']-i,:] 
    
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    fig.subplots_adjust(left=0.2)
    fig.subplots_adjust(bottom=0.2)
    plt.contourf(kzgrid[:par['nkz0']/2-1],herm_grid,np.transpose(eterm_kzn_out),100)
    plt.title(file_base)
    plt.xlabel(r'$k_z L_n$',size=13)
    plt.ylabel(r'Hermite n',size=13)
    time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
    #plt.figtext(0.02,0.03,time_string)
    plt.colorbar()
    plt.show()

    file_name=par['diagdir'][1:-1]+'/'+file_base+'_kpn.dat'
    np.savetxt(file_name,eterm_kpn)
    file_name=par['diagdir'][1:-1]+'/'+file_base+'_kpn_grid.dat'
    np.savetxt(file_name,kpgrid)

    file_name=par['diagdir'][1:-1]+'/'+file_base+'_kzn.dat'
    np.savetxt(file_name,eterm_kzn_out)
    file_name=par['diagdir'][1:-1]+'/'+file_base+'_kzn_grid.dat'
    np.savetxt(file_name,kzgrid[0:par['nkz0']/2-1])


def read_time_step_energy3d(which_time,swap_endian=False):
   """Reads a time step from energy3d.dat.  Time step determined by \'which_itime\'"""
   f = open('energy3d.dat','rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*4
   #In file:
   #1: FE
   #2: Drive
   #3: Coll and Hcoll
   #4: Other dissipation
   mem_tot=ntot*16
   en0=np.empty((par['nkz0'],par['nky0'],par['nkx0'],4))
   f.seek(8+which_itime*(8+mem_tot))
   en0=np.fromfile(f,dtype='float64',count=ntot)
   if swap_endian:
       en0=en0.newbyteorder()
   #print sum(gt0)
   f.close()
   return en0


def get_energy_spectra0(start_time=-1.0,end_time=-1.0,data_out=False,which_term=-1,fit_spectra=False,plot_7o3=False):
    """Plots time-averaged energy-related spectra from energy3d.dat\n
    start_time: start time for time average (defaults to 0.0) \n
    end_time: end time for time average (defaults to end) \n
    data_out: flag for outputting data in an ascii file--useful \n
    for replotting or quick viewing later. \n
    which_term=-1 (default) ==> all
    which_term=1  ==> Free energy
    which_term=2  ==> Drive
    which_term=3  ==> Dissipation 
    which_term=4  ==> hyper-diffusion """

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    #print "kxgrid",kxgrid
    #print "kygrid",kygrid
    #print "kzgrid",kzgrid
    time=get_time_from_energy()

    kxgrid_full=np.empty(par['nkx0']*2-1)
    for i in range(par['nkx0']):
        kxgrid_full[i+par['nkx0']-1]=kxgrid[i]
    for i in range(par['nkx0']-1):
        kxgrid_full[i]=-1.0*kxgrid[par['nkx0']-1-i]
    print kxgrid_full

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    diss_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    entk_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    drive_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    hyp_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    #en_in=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    ntime=iend-istart+1
    J0a=get_gyroavg()
    ntot=par['nkx0']*par['nky0']*par['nkz0']
    mem_tot=ntot*8
    #Note in energy3d.dat:
    #Energy
    #Drive
    #Coll+Hyp Coll
    #Hyp's
    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('energy3d.dat','rb')
        f.seek(i*(8+4*mem_tot))
        ttemp=np.fromfile(f,dtype='float64',count=1)
        print ttemp
        f.seek(8+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        entk_sum=entk_sum+en_in
        f.seek(8+mem_tot+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        drive_sum=drive_sum+en_in
        f.seek(8+2*mem_tot+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        diss_sum=diss_sum+en_in
        f.seek(8+3*mem_tot+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        hyp_sum=hyp_sum+en_in
        f.close()

    kx_out,ky_out,kz_out,herm_out=get_grids_shifted()
    print "ky_out",ky_out

    if plot_7o3:
        spectkx_n7o3=np.arange(len(kx_out),dtype='float')
        spectky_n7o3=np.arange(len(ky_out),dtype='float')
        for i in range(len(kx_out)):
            spectkx_n7o3[i]=kx_out[i]**(-7.0/3.0)
            #spectkx_n7o3[i]=kx_out[i]**(-7.0/3.0)*np.e**(-2.0/3.0*kx_out[i]**2)
        for i in range(len(ky_out)):
            spectky_n7o3[i]=ky_out[i]**(-7.0/3.0)
            #spectky_n7o3[i]=ky_out[i]**(-7.0/3.0)*np.e**(-2.0/3.0*ky_out[i]**2)
        kxp5_loc=np.argmin(np.abs(kx_out-0.6))
        print "kxp5_loc",kxp5_loc
        kyp5_loc=np.argmin(np.abs(ky_out-0.6))
        print "kyp5_loc",kyp5_loc
        kx1_loc=np.argmin(np.abs(kx_out-1.2))
        ky1_loc=np.argmin(np.abs(ky_out-1.2))
        #print "ky_out",ky_out

    e2d_plot=np.empty((par['nkx0']*2-1,par['nky0']))
    
    if which_term== 3 or which_term==-1:
      diss=diss_sum/float(ntime)
      diss_kza=np.sum(diss,2)
      np.info(diss_kza)
      e2d_plot=plot_prep2d(diss_kza)

      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
      plt.title('Dissipation')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.ylabel(r'$k_y\rho_i$',size=13)
      time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
      plt.figtext(0.02,0.03,time_string)
      plt.colorbar()
      plt.show()
      if data_out:
          file_name=par['diagdir'][1:-1]+'/es_diss0_kxky.dat'
          np.savetxt(file_name,e2d_plot)


    if which_term==2 or which_term==-1:
      drive=drive_sum/float(ntime)
      drive_kza=np.sum(drive,2)
      np.info(drive_kza)
      e2d_plot=plot_prep2d(drive_kza)
      #drive_kza=np.roll(drive_kza,par['nky0']/2-1,axis=1)
      #plt.contour(ky_out,kx_out,drive_kza,200)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
      plt.title('Drive')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.ylabel(r'$k_y\rho_i$',size=13)
      plt.colorbar()
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      if data_out:
          file_name=par['diagdir'][1:-1]+'/es_drive0_kxky.dat'
          np.savetxt(file_name,e2d_plot)

      drive_kykz=drive[0,:,:]
      drive_kykz=np.roll(drive_kykz,par['nky0']/2-1,axis=0)
      drive_kykz=np.roll(drive_kykz,par['nkz0']/2-1,axis=1)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kz_out,ky_out,drive_kykz,200)
      plt.title('Drive')
      plt.xlabel(r'$k_z R$')
      plt.ylabel(r'$k_y \rho_i$')
      plt.colorbar()
      plt.show()

    if which_term==1 or which_term==-1:
      entk_sum=entk_sum/float(ntime)
      entk_kza=np.sum(entk_sum,axis=2)
      e2d_plot=plot_prep2d(entk_kza)
      #entk_kza=np.roll(entk_kza,par['nky0']/2-1,axis=1)
      #plt.contour(ky_out,kx_out,entk_kza,200)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
      plt.title('Free Energy')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.ylabel(r'$k_y\rho_i$',size=13)
      time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
      plt.figtext(0.02,0.02,time_string)
      plt.colorbar()
      plt.show()
      if data_out:
          file_name=par['diagdir'][1:-1]+'/es_fenergy0_kxky.dat'
          np.savetxt(file_name,e2d_plot)

      entk_kykz=entk_sum[0,:,:]
      entk_kykz=np.roll(entk_kykz,par['nky0']/2-1,axis=0)
      entk_kykz=np.roll(entk_kykz,par['nkz0']/2-1,axis=1)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kz_out,ky_out,entk_kykz,200)
      plt.title('Free Energy')
      plt.xlabel(r'$k_z R$')
      plt.ylabel(r'$k_y \rho_i$')
      plt.colorbar()
      plt.show()

      kz_eff = np.empty(len(ky_out))
      for i in range(len(ky_out)):
          kz_eff[i]=np.sum(np.abs(kz_out[:])*entk_kykz[i,:])/np.sum(entk_kykz[i,:])
      plt.figure(figsize=(4.0,3.0))
      plt.subplots_adjust(bottom=0.15)
      plt.subplots_adjust(left=0.2)
      plt.plot(ky_out,kz_eff,'x-',linewidth=2)
      plt.xlabel(r'$k_y \rho_i$')
      plt.ylabel(r'$k_z^{eff}L$')
      ax=plt.axis()
      plt.axis((0,1.0,0,1.0))
      plt.show()

      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kz_out,ky_out,np.log(entk_kykz),200)
      plt.title('Free Energy')
      plt.xlabel(r'$k_z R$')
      plt.ylabel(r'$k_y \rho_i$')
      plt.colorbar()
      plt.show()


    if which_term==4 or which_term==-1:
      hyp_sum=hyp_sum/float(ntime)
      hyp_kza=np.sum(hyp_sum,axis=2)
      e2d_plot=plot_prep2d(hyp_kza)
      #entk_kza=np.roll(entk_kza,par['nky0']/2-1,axis=1)
      #plt.contour(ky_out,kx_out,entk_kza,200)
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
      plt.title('Hyps')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.ylabel(r'$k_y\rho_i$',size=13)
      time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
      plt.figtext(0.02,0.02,time_string)
      plt.colorbar()
      plt.show()
      if data_out:
          file_name=par['diagdir'][1:-1]+'/es_hyp0_kxky.dat'
          np.savetxt(file_name,e2d_plot)



    if which_term==3 or which_term==-1:
      diss_kx=np.sum(diss_kza,1)
      diss_kx[1:]=2.0*diss_kx[1:]
      diss_ky=np.sum(diss_kza[:,:],0)
      for i in range(1,par['nky0']/2):
          diss_ky[i]=diss_ky[i]+diss_ky[par['nky0']-i]
          diss_ky[par['nky0']-i]=diss_ky[i]
      diss_ky=np.roll(diss_ky,par['nky0']/2-1,axis=0)
      diss_kya=np.sum(diss,1)
      diss_kz=np.sum(diss_kya[:,:],0)
      for i in range(1,par['nkz0']/2):
          diss_kz[i]=diss_kz[i]+diss_kz[par['nkz0']-i]
          diss_kz[par['nkz0']-i]=diss_kz[i]
      diss_kz=np.roll(diss_kz,par['nkz0']/2-1,axis=0)

      #Linear plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kx_out,diss_kx,'x-')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.title('Dissipation')
      plt.xlim((0,par['kxmax0']))
      plt.show()



      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(ky_out,diss_ky,'x-')
      plt.title('Dissipation')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kymax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kz_out,diss_kz,'x-')
      plt.title('Dissipation')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kzmax0']))
      plt.show()


      if data_out:
          file_name=par['diagdir'][1:-1]+'/es_diss0_kx.dat'
          arr_out=np.empty((par['nkx0'],2)) 
          arr_out[:,0]=kx_out
          arr_out[:,1]=diss_kx
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/es_diss0_ky.dat'
          arr_out=np.empty((par['nky0'],2)) 
          arr_out[:,0]=ky_out
          arr_out[:,1]=diss_ky
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/es_diss0_kz.dat'
          arr_out=np.empty((par['nkz0'],2)) 
          arr_out[:,0]=kz_out
          arr_out[:,1]=diss_kz
          np.savetxt(file_name,arr_out)

      #log-log plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kx_out,-1.0*diss_kx,'x-',basex=10,basey=10)
      plt.title('Dissipation')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(ky_out,-1.0*diss_ky,'x-',basex=10,basey=10)
      plt.title('Dissipation')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kz_out,-1.0*diss_kz,'x-',basex=10,basey=10)
      plt.title('Dissipation')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()

    if which_term==2 or which_term==-1:
      drive_kx=np.sum(drive_kza,1)
      drive_kx[1:]=2.0*drive_kx[1:]
      drive_ky=np.sum(drive_kza[:,:],0)
      for i in range(1,par['nky0']/2):
          drive_ky[i]=drive_ky[i]+drive_ky[par['nky0']-i]
          drive_ky[par['nky0']-i]=drive_ky[i]
      drive_ky=np.roll(drive_ky,par['nky0']/2-1,axis=0)
      drive_kya=np.sum(drive,1)
      drive_kz=np.sum(drive_kya[:,:],0)
      for i in range(1,par['nkz0']/2):
          drive_kz[i]=drive_kz[i]+drive_kz[par['nkz0']-i]
          drive_kz[par['nkz0']-i]=drive_kz[i]
      drive_kz=np.roll(drive_kz,par['nkz0']/2-1,axis=0)

      if data_out:
          file_name=par['diagdir'][1:-1]+'/es_drive0_kx.dat'
          arr_out=np.empty((par['nkx0'],2)) 
          arr_out[:,0]=kx_out
          arr_out[:,1]=drive_kx
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/es_drive0_ky.dat'
          arr_out=np.empty((par['nky0'],2)) 
          arr_out[:,0]=ky_out
          arr_out[:,1]=drive_ky
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/es_drive0_kz.dat'
          arr_out=np.empty((par['nkz0'],2)) 
          arr_out[:,0]=kz_out
          arr_out[:,1]=drive_kz
          np.savetxt(file_name,arr_out)



    #drive_kx=np.sum(drive_kza,1)
    #drive_ky=drive_kza[0,:]+2.0*np.sum(drive_kza[1:,:],0)
    ##drive_ky=np.roll(drive_ky,par['nky0']/2-1,axis=0)
    #drive_kya=np.sum(drive,1)
    #drive_kz=drive_kya[0,:]+2.0*np.sum(drive_kya[1:,:],0)
    #drive_kz=np.roll(drive_kz,par['nkz0']/2-1,axis=0)

      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kx_out,drive_kx,'x-')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.title('Drive')
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kxmax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(ky_out,drive_ky,'x-')
      plt.title('Drive')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kymax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kz_out,drive_kz,'x-')
      plt.title('Drive')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kzmax0']))
      plt.show()


      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kx_out,drive_kx,'x-',basex=10,basey=10)
      plt.title('Drive')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(ky_out,drive_ky,'x-',basex=10,basey=10)
      plt.title('Drive')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kz_out,drive_kz,'x-',basex=10,basey=10)
      plt.title('Drive')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()


    if which_term==1 or which_term==-1:
      entk_kx=np.sum(entk_kza,1)
      entk_kx[1:]=2.0*entk_kx[1:]
      entk_ky=np.sum(entk_kza[:,:],0)
      for i in range(1,par['nky0']/2):
          entk_ky[i]=entk_ky[i]+entk_ky[par['nky0']-i]
          entk_ky[par['nky0']-i]=entk_ky[i]
      entk_ky=np.roll(entk_ky,par['nky0']/2-1,axis=0)
      entk_kya=np.sum(entk_sum,1)
      entk_kz=np.sum(entk_kya[:,:],0)
      for i in range(1,par['nkz0']/2):
          entk_kz[i]=entk_kz[i]+entk_kz[par['nkz0']-i]
          entk_kz[par['nkz0']-i]=entk_kz[i]
      entk_kz=np.roll(entk_kz,par['nkz0']/2-1,axis=0)

      if data_out:
          file_name=par['diagdir'][1:-1]+'/es_fenergy_kx.dat'
          arr_out=np.empty((par['nkx0'],2)) 
          arr_out[:,0]=kx_out
          arr_out[:,1]=entk_kx
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/es_fenergy_ky.dat'
          arr_out=np.empty((par['nky0'],2)) 
          arr_out[:,0]=ky_out
          arr_out[:,1]=entk_ky
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/es_fenergy_kz.dat'
          arr_out=np.empty((par['nkz0'],2)) 
          arr_out[:,0]=kz_out
          arr_out[:,1]=entk_kz
          np.savetxt(file_name,arr_out)



      #entk_kx=np.sum(entk_kza,1)
      #entk_ky=entk_kza[0,:]+2.0*np.sum(entk_kza[1:,:],0)
      ##diss_ky=np.roll(diss_ky,par['nky0']/2-1,axis=0)
      #entk_kya=np.sum(entk_sum,1)
      #entk_kz=entk_kya[0,:]+2.0*np.sum(entk_kya[1:,:],0)
      #entk_kz=np.roll(entk_kz,par['nkz0']/2-1,axis=0)

      #Linear plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kx_out,entk_kx,'x-')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.title('Free Energy')
      plt.xlim((0,par['kxmax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(ky_out,entk_ky,'x-')
      plt.title('Free Energy')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kymax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kz_out,entk_kz,'x-')
      plt.title('Free Energy')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kzmax0']))
      plt.show()

      #log-log plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      if fit_spectra:
          print "Identify fit range from plot."
      plt.loglog(kx_out,entk_kx,'x-',basex=10,basey=10)
      if plot_7o3:
          spectkx_n7o3/=0.45*spectkx_n7o3[kxp5_loc]/entk_kx[kxp5_loc]
          plt.loglog(kx_out[kxp5_loc:kx1_loc],spectkx_n7o3[kxp5_loc:kx1_loc],':',color='black',basex=10,basey=10,label=r'$\propto k^{-7/3}$')
          plt.legend(loc='lower left')
      plt.title('Free Energy')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      if fit_spectra:
          startx=float(raw_input("Enter start value for fit:"))
          endx=float(raw_input("Enter end value for fit:"))

          fit=fit_function(kx_out,entk_kx,(1.0,-1.0),startx=startx,endx=endx,which_func='power')
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.22)
          plt.loglog(kx_out,entk_kx,'x-',basex=10,basey=10)
          if startx==-1:
              sx0=0
          else:
              sx0=argmin(abs(kx_out-startx))
         #     sx0=sx0-sx0/4
          if endx==-1:
              ex0=0
          else:
              ex0=argmin(abs(kx_out-endx))
          fit_x=kx_out[sx0:ex0]
          fit_func=2.0*fit[0]*fit_x**fit[1]
          plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
          plt.annotate(r'$\propto k$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
          plt.title('Free Energy')
          plt.xlabel(r'$k_x\rho_i$',size=13)
          plt.figtext(0.02,0.02,time_string)
          plt.show()
           
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      if fit_spectra:
          print "Identify fit range from plot."
      plt.loglog(ky_out,entk_ky,'x-',basex=10,basey=10)
      if plot_7o3:
          spectky_n7o3/=0.45*spectky_n7o3[kyp5_loc]/entk_ky[kyp5_loc]
          plt.loglog(ky_out[kyp5_loc:ky1_loc],spectky_n7o3[kyp5_loc:ky1_loc],':',color='black',basex=10,basey=10,label=r'$\propto k^{-7/3}$')
          plt.legend(loc='lower left')
      plt.title('Free Energy')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      if fit_spectra:
          startx=float(raw_input("Enter start value for fit:"))
          endx=float(raw_input("Enter end value for fit:"))
          fit=fit_function(ky_out,entk_ky,(1.0,-1.0),startx=startx,endx=endx,which_func='power')
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.22)
          plt.loglog(ky_out,entk_ky,'x-',basex=10,basey=10)
          if startx==-1:
              sx0=0
          else:
              sx0=argmin(abs(ky_out-startx))
          if endx==-1:
              ex0=0
          else:
              ex0=argmin(abs(ky_out-endx))
          #print "sx0,ex0",sx0,ex0
          #print len(ky_out)
          #print "ky_out",ky_out
          fit_x=ky_out[sx0:ex0]
          fit_func=2.0*fit[0]*fit_x**fit[1]
          #print fit_x
          #print fit_func
          plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
          plt.annotate(r'$\propto k$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
          plt.title('Free Energy')
          plt.xlabel(r'$k_y\rho_i$',size=13)
          plt.figtext(0.02,0.02,time_string)
          plt.show()

      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      if fit_spectra:
          print "Identify fit range from plot."
      plt.loglog(kz_out,entk_kz,'x-',basex=10,basey=10)
      plt.title('Free Energy')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      if fit_spectra:
          startx=float(raw_input("Enter start value for fit:"))
          endx=float(raw_input("Enter end value for fit:"))
          fit=fit_function(kz_out,entk_kz,(1.0,-1.0),startx=startx,endx=endx,which_func='power')
          plt.figure(figsize=(4.5,3.0))
          fig=plt.gcf()
          fig.subplots_adjust(left=0.2)
          fig.subplots_adjust(bottom=0.22)
          plt.loglog(kz_out,entk_kz,'x-',basex=10,basey=10)
          if startx==-1:
              sx0=0
          else:
              sx0=argmin(abs(kz_out-startx))
          if endx==-1:
              ex0=0
          else:
              ex0=argmin(abs(kz_out-endx))
          #print "sx0,ex0",sx0,ex0
          #print len(kz_out)
          #print "kz_out",kz_out
          fit_x=kz_out[sx0:ex0]
          fit_func=2.0*fit[0]*fit_x**fit[1]
          #print fit_x
          #print fit_func
          plt.loglog(fit_x,fit_func,'--',basex=10,basey=10)
          plt.annotate(r'$\propto k$'+'^'+str(fit[1])[:5],(fit_x[len(fit_x)/2],fit_func[len(fit_x)/2]))
          plt.title('Free Energy')
          plt.xlabel(r'$k_y\rho_i$',size=13)
          plt.figtext(0.02,0.02,time_string)
          plt.show()


    if which_term==4 or which_term==-1:
      hyp_kx=np.sum(hyp_kza,1)
      hyp_kx[1:]=2.0*hyp_kx[1:]
      hyp_ky=np.sum(hyp_kza[:,:],0)
      for i in range(1,par['nky0']/2):
          hyp_ky[i]=hyp_ky[i]+hyp_ky[par['nky0']-i]
          hyp_ky[par['nky0']-i]=hyp_ky[i]
      hyp_ky=np.roll(hyp_ky,par['nky0']/2-1,axis=0)
      hyp_kya=np.sum(hyp_sum,1)
      hyp_kz=np.sum(hyp_kya[:,:],0)
      for i in range(1,par['nkz0']/2):
          hyp_kz[i]=hyp_kz[i]+hyp_kz[par['nkz0']-i]
          hyp_kz[par['nkz0']-i]=hyp_kz[i]
      hyp_kz=np.roll(hyp_kz,par['nkz0']/2-1,axis=0)
  
      if data_out:
          file_name=par['diagdir'][1:-1]+'/es_hyp0_kx.dat'
          arr_out=np.empty((par['nkx0'],2)) 
          arr_out[:,0]=kx_out
          arr_out[:,1]=hyp_kx
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/es_hyp0_ky.dat'
          arr_out=np.empty((par['nky0'],2)) 
          arr_out[:,0]=ky_out
          arr_out[:,1]=hyp_ky
          np.savetxt(file_name,arr_out)
          file_name=par['diagdir'][1:-1]+'/es_hyp0_kz.dat'
          arr_out=np.empty((par['nkz0'],2)) 
          arr_out[:,0]=kz_out
          arr_out[:,1]=hyp_kz
          np.savetxt(file_name,arr_out)

      #Linear plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kx_out,hyp_kx,'x-')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.title('Hyps')
      plt.xlim((0,par['kxmax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(ky_out,hyp_ky,'x-')
      plt.title('Hyps')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kymax0']))
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.2)
      plt.plot(kz_out,hyp_kz,'x-')
      plt.title('Hyps')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.xlim((0,par['kzmax0']))
      plt.show()

      #log-log plot
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kx_out,-hyp_kx,'x-',basex=10,basey=10)
      plt.title('Hyps')
      plt.xlabel(r'$k_x\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(ky_out,-hyp_ky,'x-',basex=10,basey=10)
      plt.title('Hyps')
      plt.xlabel(r'$k_y\rho_i$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()
      plt.figure(figsize=(4.5,3.0))
      fig=plt.gcf()
      fig.subplots_adjust(left=0.2)
      fig.subplots_adjust(bottom=0.22)
      plt.loglog(kz_out,-hyp_kz,'x-',basex=10,basey=10)
      plt.title('Hyps')
      plt.xlabel(r'$k_z R$',size=13)
      plt.figtext(0.02,0.02,time_string)
      plt.show()

def get_energy_spectra0_red(start_time=-1.0,end_time=-1.0,nkx0_red = -1, \
                    nky0_red = -1, nkz0_red = -1,data_out=False):
    """Take a DNS and reduce its box to a LES domain """

    #For reducing the size of the array
    if nkx0_red==-1.0:
       nkx0_red=par['nkx0']

    if nky0_red==-1.0:
       nky0_red=par['nky0']

    if nkz0_red==-1.0:
       nkz0_red=par['nkz0']

    kx_red_l = par['nkx0']    -nkx0_red 
    kx_red_h = par['nkx0']    + nkx0_red - 1
    ky_red_l = par['nky0']/2  - nky0_red/2
    ky_red_h = par['nky0']/2  + nky0_red/2 -1
    kz_red_l = par['nkz0']/2  - nkz0_red/2
    kz_red_h = par['nkz0']/2  + nkz0_red/2 -1


    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    
    kxgrid_full=np.empty(par['nkx0']*2-1)
    for i in range(par['nkx0']):
        kxgrid_full[i+par['nkx0']-1]=kxgrid[i]
    for i in range(par['nkx0']-1):
        kxgrid_full[i]=-1.0*kxgrid[par['nkx0']-1-i]
    print kxgrid_full


    time=get_time_from_energy()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    diss_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    entk_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    drive_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    hyp_sum=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    #en_in=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    ntime=iend-istart+1
    J0a=get_gyroavg()
    ntot=par['nkx0']*par['nky0']*par['nkz0']
    mem_tot=ntot*8
    #Note in energy3d.dat:
    #Energy
    #Drive
    #Coll+Hyp Coll
    #Hyp's
    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('energy3d.dat','rb')
        f.seek(i*(8+4*mem_tot))
        ttemp=np.fromfile(f,dtype='float64',count=1)
        print ttemp
        f.seek(8+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        entk_sum=entk_sum+en_in
        f.seek(8+mem_tot+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        drive_sum=drive_sum+en_in
        f.seek(8+2*mem_tot+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        diss_sum=diss_sum+en_in
        f.seek(8+3*mem_tot+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        hyp_sum=hyp_sum+en_in
        f.close()

    time_string='Time=['+str(time[istart])[:5]+','+str(time[iend])[:5]+']'
    
    kx_out,ky_out,kz_out,herm_out=get_grids_shifted()
    print "ky_out",ky_out
    print "kz_out",kz_out
    print "kx_out",kx_out
    
    #Take only the positive part (for plotting purposes)
    kx_out_pos = kxgrid_full[par['nkx0']-1:kx_red_h]
    ky_out_pos = ky_out[par['nky0']/2-1:ky_red_h]
    kz_out_pos = kz_out[par['nkz0']/2-1:kz_red_h]


    diss=diss_sum/float(ntime)
    drive=drive_sum/float(ntime)
    entk_sum=entk_sum/float(ntime)
    hyp_sum=hyp_sum/float(ntime)

    diss_3d = les_prep3d(diss,kx_red_l,kx_red_h,ky_red_l,ky_red_h+1,kz_red_l,kz_red_h+1)
    diss_kx_pos,diss_ky_pos,diss_kz_pos= les_prep1d(diss_3d,nkx0_red,nky0_red,nkz0_red)

#    #Linear plot
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(kx_out_pos,diss_kx_pos,'r-')
#    plt.xlabel(r'$k_x\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.title('Dissipation')
# #   plt.xlim((0,par['kxmax0']))
#    plt.show()
#

#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(ky_out_pos,diss_ky_pos,'r-')
#    plt.title('Dissipation')
#    plt.xlabel(r'$k_y\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
##    plt.xlim((0,par['kymax0']))
#    plt.show()
#
#
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(kz_out_pos,diss_kz_pos,'r-')
#    plt.title('Dissipation')
#    plt.xlabel(r'$k_z R$',size=13)
#    plt.figtext(0.02,0.02,time_string)
##    plt.xlim((0,par['kzmax0']))
#    plt.show()
#

    if data_out:
        file_name=par['diagdir'][1:-1]+'/diss0_kx_red.dat'
        arr_out=np.empty((nkx0_red,2)) 
        arr_out[:,0]=kx_out_pos
        arr_out[:,1]=diss_kx_pos
        np.savetxt(file_name,arr_out)
        file_name=par['diagdir'][1:-1]+'/diss0_ky_red.dat'
        arr_out=np.empty((nky0_red/2,2)) 
        arr_out[:,0]=ky_out_pos
        arr_out[:,1]=diss_ky_pos
        np.savetxt(file_name,arr_out)
        file_name=par['diagdir'][1:-1]+'/diss0_kz_red.dat'
        arr_out=np.empty((nkz0_red/2,2)) 
        arr_out[:,0]=kz_out_pos
        arr_out[:,1]=diss_kz_pos
        np.savetxt(file_name,arr_out)

##    #log-log plot
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(kx_out_pos,-1.0*diss_kx_pos,basex=10,basey=10)
#    plt.title('Dissipation')
#    plt.xlabel(r'$k_x\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(ky_out_pos,-1.0*diss_ky_pos,basex=10,basey=10)
#    plt.title('Dissipation')
#    plt.xlabel(r'$k_y\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(kz_out_pos,-1.0*diss_kz_pos,basex=10,basey=10)
#    plt.title('Dissipation')
#    plt.xlabel(r'$k_z R$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#
#
    drive_3d = les_prep3d(drive,kx_red_l,kx_red_h,ky_red_l,ky_red_h+1,kz_red_l,kz_red_h+1)
    drive_kx_pos,drive_ky_pos,drive_kz_pos= les_prep1d(drive_3d,nkx0_red,nky0_red,nkz0_red)

##    #Linear plot
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(kx_out_pos,drive_kx_pos,'r-')
#    plt.xlabel(r'$k_x\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.title('Drive')
# #   plt.xlim((0,par['kxmax0']))
#    plt.show()
#
#
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(ky_out_pos,drive_ky_pos,'r-')
#    plt.title('Drive')
#    plt.xlabel(r'$k_y\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
##    plt.xlim((0,par['kymax0']))
#    plt.show()
#
#
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(kz_out_pos,drive_kz_pos,'r-')
#    plt.title('Drive')
#    plt.xlabel(r'$k_z R$',size=13)
#    plt.figtext(0.02,0.02,time_string)
##    plt.xlim((0,par['kzmax0']))
#    plt.show()
#

    if data_out:
        file_name=par['diagdir'][1:-1]+'/drive0_kx_red.dat'
        arr_out=np.empty((nkx0_red,2)) 
        arr_out[:,0]=kx_out_pos
        arr_out[:,1]=drive_kx_pos
        np.savetxt(file_name,arr_out)
        file_name=par['diagdir'][1:-1]+'/drive0_ky_red.dat'
        arr_out=np.empty((nky0_red/2,2)) 
        arr_out[:,0]=ky_out_pos
        arr_out[:,1]=drive_ky_pos
        np.savetxt(file_name,arr_out)
        file_name=par['diagdir'][1:-1]+'/drive0_kz_red.dat'
        arr_out=np.empty((nkz0_red/2,2)) 
        arr_out[:,0]=kz_out_pos
        arr_out[:,1]=drive_kz_pos
        np.savetxt(file_name,arr_out)


   
#    #log-log plot
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(kx_out_pos,drive_kx_pos,basex=10,basey=10)
#    plt.title('Drive')
#    plt.xlabel(r'$k_x\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(ky_out_pos,drive_ky_pos,basex=10,basey=10)
#    plt.title('Drive')
#    plt.xlabel(r'$k_y\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(kz_out_pos,drive_kz_pos,basex=10,basey=10)
#    plt.title('Drive')
#    plt.xlabel(r'$k_z R$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#
    entk_3d = les_prep3d(entk_sum,kx_red_l,kx_red_h,ky_red_l,ky_red_h+1,kz_red_l,kz_red_h+1)
    entk_kx_pos,entk_ky_pos,entk_kz_pos= les_prep1d(entk_3d,nkx0_red,nky0_red,nkz0_red)

#    #Linear plot
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(kx_out_pos,entk_kx_pos,'r-')
#    plt.xlabel(r'$k_x\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.title('Free energy')
# #   plt.xlim((0,par['kxmax0']))
#    plt.show()
#
#
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(ky_out_pos,entk_ky_pos,'r-')
#    plt.title('Free energy')
#    plt.xlabel(r'$k_y\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
##    plt.xlim((0,par['kymax0']))
#    plt.show()
#
#
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(kz_out_pos,entk_kz_pos,'r-')
#    plt.title('Free energy')
#    plt.xlabel(r'$k_z R$',size=13)
#    plt.figtext(0.02,0.02,time_string)
##    plt.xlim((0,par['kzmax0']))
#    plt.show()
#

    if data_out:
        file_name=par['diagdir'][1:-1]+'/entk0_kx_red.dat'
        arr_out=np.empty((nkx0_red,2)) 
        arr_out[:,0]=kx_out_pos
        arr_out[:,1]=entk_kx_pos
        np.savetxt(file_name,arr_out)
        file_name=par['diagdir'][1:-1]+'/entk0_ky_red.dat'
        arr_out=np.empty((nky0_red/2,2)) 
        arr_out[:,0]=ky_out_pos
        arr_out[:,1]=entk_ky_pos
        np.savetxt(file_name,arr_out)
        file_name=par['diagdir'][1:-1]+'/entk0_kz_red.dat'
        arr_out=np.empty((nkz0_red/2,2)) 
        arr_out[:,0]=kz_out_pos
        arr_out[:,1]=entk_kz_pos
        np.savetxt(file_name,arr_out)

#    #log-log plot
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(kx_out_pos,entk_kx_pos,basex=10,basey=10)
#    plt.title('Free energy')
#    plt.xlabel(r'$k_x\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(ky_out_pos,entk_ky_pos,basex=10,basey=10)
#    plt.title('Free energy')
#    plt.xlabel(r'$k_y\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(kz_out_pos,entk_kz_pos,basex=10,basey=10)
#    plt.title('Free energy')
#    plt.xlabel(r'$k_z R$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#
    hyp_3d = les_prep3d(hyp_sum,kx_red_l,kx_red_h,ky_red_l,ky_red_h+1,kz_red_l,kz_red_h+1)
    hyp_kx_pos,hyp_ky_pos,hyp_kz_pos= les_prep1d(hyp_3d,nkx0_red,nky0_red,nkz0_red)

#    #Linear plot
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(kx_out_pos,hyp_kx_pos,'r-')
#    plt.xlabel(r'$k_x\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.title('Hyperdiffusion')
# #   plt.xlim((0,par['kxmax0']))
#    plt.show()
#
#
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(ky_out_pos,hyp_ky_pos,'r-')
#    plt.title('Hyperdiffusion')
#    plt.xlabel(r'$k_y\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
##    plt.xlim((0,par['kymax0']))
#    plt.show()
#
#
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.2)
#    plt.plot(kz_out_pos,hyp_kz_pos,'r-')
#    plt.title('Hyperdiffusion')
#    plt.xlabel(r'$k_z R$',size=13)
#    plt.figtext(0.02,0.02,time_string)
##    plt.xlim((0,par['kzmax0']))
#    plt.show()
#

    if data_out:
        file_name=par['diagdir'][1:-1]+'/hyp0_kx_red.dat'
        arr_out=np.empty((nkx0_red,2)) 
        arr_out[:,0]=kx_out_pos
        arr_out[:,1]=hyp_kx_pos
        np.savetxt(file_name,arr_out)
        file_name=par['diagdir'][1:-1]+'/hyp0_ky_red.dat'
        arr_out=np.empty((nky0_red/2,2)) 
        arr_out[:,0]=ky_out_pos
        arr_out[:,1]=hyp_ky_pos
        np.savetxt(file_name,arr_out)
        file_name=par['diagdir'][1:-1]+'/hyp0_kz_red.dat'
        arr_out=np.empty((nkz0_red/2,2)) 
        arr_out[:,0]=kz_out_pos
        arr_out[:,1]=hyp_kz_pos
        np.savetxt(file_name,arr_out)

#    #log-log plot
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(kx_out_pos,-1*hyp_kx_pos,basex=10,basey=10)
#    plt.title('Hyperdifusion')
#    plt.xlabel(r'$k_x\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(ky_out_pos,-1*hyp_ky_pos,basex=10,basey=10)
#    plt.title('Hyperdifusion')
#    plt.xlabel(r'$k_y\rho_i$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#    plt.figure(figsize=(4.5,3.0))
#    fig=plt.gcf()
#    fig.subplots_adjust(left=0.2)
#    fig.subplots_adjust(bottom=0.22)
#    plt.loglog(kz_out_pos,-1*hyp_kz_pos,basex=10,basey=10)
#    plt.title('Hyperdifusion')
#    plt.xlabel(r'$k_z R$',size=13)
#    plt.figtext(0.02,0.02,time_string)
#    plt.show()
#
def get_energy_tracesF(start_time=-1.0,end_time=-1.0,benchmark=False):
    """Plots energy time traces from g_out.dat"""
    time=get_time_from_gout()
    #print "Warning: Must benchmark with DNA output!!!!"
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    #coll_op=np.empty((par['nkx0'],par['nky0'],par['nkz0'],par['nv0']))
    #hypcoll_op=np.empty((par['nkx0'],par['nky0'],par['nkz0'],par['nv0']))
    #for i in range(par['nv0']):
    #    coll_op[:,:,:,i]=-np.pi**0.5*par['nu']*i
    #    i_nv=float(i)/float(par['nv0'])
    #    hypcoll_op[:,:,:,i]=-np.pi**0.5*i_nv**par['hypv_order']*par['hyp_v']
    #np.info(coll_op)

    J0a=get_gyroavg()
    rhs_t=np.empty(0)
    drive_t=np.empty(0)
    coll_t=np.empty(0)
    hcoll_t=np.empty(0)
    hyps_t=np.empty(0)
    energy_t=np.empty(0)
    enphi_t=np.empty(0)
    ntime=iend-istart+1
    for t in range(istart,iend+1):
        print 'time=',time[t],' of ',time[iend]
        gt0=read_time_step_g(t)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        #gt0=np.roll(gt0,par['nky0']/2-1,axis=1)
        #gt0=np.roll(gt0,par['nkz0']/2-1,axis=2)
        phi=get_phi(gt0)
        #print 'sum(phi)',np.sum(phi)
        E_op=np.pi**0.5*np.conj(gt0)
        for i in range(par['nkx0']):
            for j in range(par['nky0']):
                E_op[i,j,:,0]=E_op[i,j,:,0]+np.pi**0.25*J0a[i,j]*np.conj(phi[i,j,:])
        #Energy
        eterm=0.5*np.real(np.sum(np.pi**0.5*np.conj(gt0)*gt0,axis=3))
        energy_t=np.append(energy_t,ksum_3dF(eterm))
        #En phi
        for i in range(par['nkx0']):
            for j in range(par['nky0']):
                eterm[i,j,:]=np.real(0.5*np.pi**0.25*np.conj(phi[i,j,:])*gt0[i,j,:,0]*J0a[i,j])
        enphi_t=np.append(enphi_t,ksum_3dF(eterm))
        #Collisions
        rhs=get_rhs_lin(gt0,1)
        eterm=np.real(np.sum(E_op*rhs,axis=3))
        #eterm=np.real(np.sum(gt0*rhs,axis=3))
        coll_t=np.append(coll_t,ksum_3dF(eterm))
        #Drive
        rhs=get_rhs_lin(gt0,6)
        eterm=np.real(np.sum(E_op*rhs,axis=3))
        #eterm=np.real(E_op[:,:,:,2]*rhs[:,:,:,2])
        drive_t=np.append(drive_t,ksum_3dF(eterm))
        #hyper-collisions
        rhs=get_rhs_lin(gt0,2)
        eterm=np.real(np.sum(E_op*rhs,axis=3))
        #print 'sum(hcol)',np.sum(rhs)
        hcoll_t=np.append(hcoll_t,ksum_3dF(eterm))
        #hyp's
        rhs=get_rhs_lin(gt0,8)
        eterm=np.real(np.sum(E_op*rhs,axis=3))
        #print 'sum(hcol)',np.sum(rhs)
        hyps_t=np.append(hyps_t,ksum_3dF(eterm))
        #All
        rhs=get_rhs_lin(gt0,0)
        eterm=np.real(np.sum(E_op*rhs,axis=3))
        #print 'sum(rhs)',np.sum(rhs)
        rhs_t=np.append(rhs_t,ksum_3dF(eterm))

    plt.plot(time[istart:iend+1],rhs_t,'-o',label='rhs')
    plt.plot(time[istart:iend+1],drive_t,'-+',label='drive')
    plt.plot(time[istart:iend+1],coll_t,'-x',label='coll')
    plt.plot(time[istart:iend+1],hcoll_t,'-.',label='hypc')
    #plt.plot(time[istart:iend+1],coll_t)
    plt.xlabel('$t(R/v_T)$')
    plt.ylabel('Energy Balance')
    plt.legend()
    plt.show()

    if benchmark:
        en_in=np.genfromtxt('energy.dat')
        plt.plot(time[istart:iend+1],energy_t,'-x',label='entropy pdiag')
        plt.plot(en_in[:,0],en_in[:,1],'-+',label='entropy diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],enphi_t,'-x',label='enphi pdiag')
        plt.plot(en_in[:,0],en_in[:,2],'-+',label='enphi diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],rhs_t,'-x',label='rhs pdiag')
        plt.plot(en_in[:,0],en_in[:,3],'-+',label='rhs diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],drive_t,'-x',label='drive pdiag')
        plt.plot(en_in[:,0],en_in[:,4],'-+',label='drive diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],coll_t,'-x',label='coll pdiag')
        plt.plot(en_in[:,0],en_in[:,5],'-+',label='coll diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],hcoll_t,'-x',label='hcoll pdiag')
        plt.plot(en_in[:,0],en_in[:,6],'-+',label='hcoll diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],hyps_t,'-x',label='hyps pdiag')
        plt.plot(en_in[:,0],en_in[:,7],'-+',label='hyps diag')
        plt.legend()
        plt.show()

        enistart=np.argmin(abs(en_in[:,0]-start_time))
        eniend=np.argmin(abs(en_in[:,0]-end_time))
        entime=en_in[enistart:eniend+1,0]
        encoll=en_in[enistart:eniend+1,5]
        plt.plot(entime,coll_t/encoll,label='pcoll/coll')
        #plt.legend()
        #plt.show()
        enhcoll=en_in[enistart:eniend+1,6]
        plt.plot(entime,hcoll_t/enhcoll,label='phcoll/hcoll')
        plt.legend()
        plt.show()

    return rhs_t,drive_t,coll_t,hcoll_t

def get_energy_traces(start_time=-1.0,end_time=-1.0):
    """Plots energy time traces from energy.dat
    Returns an array in the same form as energy.dat"""

    en_in=np.genfromtxt('energy.dat')
    time=en_in[:,0]
    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop
    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    #coll_op=np.empty((par['nkx0'],par['nky0'],par['nkz0'],par['nv0']))
    #hypcoll_op=np.empty((par['nkx0'],par['nky0'],par['nkz0'],par['nv0']))
    #for i in range(par['nv0']):
    #    coll_op[:,:,:,i]=-np.pi**0.5*par['nu']*i
    #    i_nv=float(i)/float(par['nv0'])
    #    hypcoll_op[:,:,:,i]=-np.pi**0.5*i_nv**par['hypv_order']*par['hyp_v']
    #np.info(coll_op)

    plt.plot(time[istart:iend+1],en_in[istart:iend+1,1],label='Entropy')
    plt.plot(time[istart:iend+1],en_in[istart:iend+1,2],label='ES Energy')
    plt.plot(time[istart:iend+1],en_in[istart:iend+1,3],label='dE/dt Total')
    plt.plot(time[istart:iend+1],en_in[istart:iend+1,4],label='Drive')
    plt.plot(time[istart:iend+1],en_in[istart:iend+1,5],label='Collisions')
    plt.plot(time[istart:iend+1],en_in[istart:iend+1,6],label='Hyper-collisions')
    plt.plot(time[istart:iend+1],en_in[istart:iend+1,7],label='Hyper-diffusion')
    plt.plot(time[istart:iend+1],en_in[istart:iend+1,8],label='Nonlinearity')
    #plt.plot(time[istart:iend+1],coll_t)
    plt.xlabel('$t(R/v_T)$')
    plt.ylabel('Energy Balance')
    plt.legend()
    plt.show()

    return en_in

def get_energy_traces0(start_time=-1.0,end_time=-1.0,benchmark=False):
    """Plots energy time traces from energy3d.dat \n
    benchmark=True benchmarks agains energy.dat"""
    time=get_time_from_energy()
    #print "Warning: Must benchmark with DNA output!!!!"
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))

    #coll_op=np.empty((par['nkx0'],par['nky0'],par['nkz0'],par['nv0']))
    #hypcoll_op=np.empty((par['nkx0'],par['nky0'],par['nkz0'],par['nv0']))
    #for i in range(par['nv0']):
    #    coll_op[:,:,:,i]=-np.pi**0.5*par['nu']*i
    #    i_nv=float(i)/float(par['nv0'])
    #    hypcoll_op[:,:,:,i]=-np.pi**0.5*i_nv**par['hypv_order']*par['hyp_v']
    #np.info(coll_op)

    J0a=get_gyroavg()
    drive_t=np.empty(0)
    diss_t=np.empty(0)
    hyps_t=np.empty(0)
    energy_t=np.empty(0)
    ntime=iend-istart+1
    ntot=par['nkx0']*par['nky0']*par['nkz0']
    mem_tot=ntot*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('energy3d.dat','rb')
        f.seek(i*(8+4*mem_tot))
        ttemp=np.fromfile(f,dtype='float64',count=1)
        print ttemp
        time[i]=ttemp
        f.seek(8+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        energy_t=np.append(energy_t,ksum_3dF(en_in))

        f.seek(8+mem_tot+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        drive_t=np.append(drive_t,ksum_3dF(en_in))

        f.seek(8+2*mem_tot+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        diss_t=np.append(diss_t,ksum_3dF(en_in))
        #diss_sum=diss_sum+en_in

        f.seek(8+3*mem_tot+i*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        hyps_t=np.append(hyps_t,ksum_3dF(en_in))
        f.close()

    plt.plot(time[istart:iend+1],energy_t,'-.',label='energy')
    plt.plot(time[istart:iend+1],drive_t,'-+',label='drive')
    plt.plot(time[istart:iend+1],diss_t,'-x',label='diss')
    plt.plot(time[istart:iend+1],hyps_t,'-.',label='hyps')
    #plt.plot(time[istart:iend+1],coll_t)
    plt.xlabel('$t(R/v_T)$')
    plt.ylabel('Energy Balance')
    plt.legend()
    plt.show()

    if benchmark:
        en_in=np.genfromtxt('energy.dat')
        plt.plot(time[istart:iend+1],energy_t,'-x',label='free energy pdiag')
        plt.plot(en_in[:,0],en_in[:,1]+en_in[:,2],'-+',label='free energy diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],drive_t+diss_t+hyps_t,'-x',label='rhs pdiag')
        plt.plot(en_in[:,0],en_in[:,3],'-+',label='rhs diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],drive_t,'-x',label='drive pdiag')
        plt.plot(en_in[:,0],en_in[:,4],'-+',label='drive diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],diss_t,'-x',label='coll+h pdiag')
        plt.plot(en_in[:,0],en_in[:,5]+en_in[:,6],'-+',label='coll+h diag')
        plt.legend()
        plt.show()
        plt.plot(time[istart:iend+1],hyps_t,'-x',label='hyps pdiag')
        plt.plot(en_in[:,0],en_in[:,7],'-+',label='hyps diag')
        plt.legend()
        plt.show()

    #return rhs_t,drive_t,coll_t,hcoll_t

def get_entropy(g_in,sum_n=True):
   """Gets sum 0.5*pi**0.5 conjg(g)*g for 4D distribution function"""
   if sum_n:
       entropy=np.sum(0.5*np.pi**0.5*np.conj(g_in)*g_in,axis=3)
   else:
       entropy=0.5*np.pi**0.5*np.conj(g_in)*g_in
   return entropy

def get_fe_es(g_in):
   """Gets electrostatic part of the free energy."""
   phi_bar=get_phi_bar(g_in)
   fe_es=0.5*np.pi**0.25*np.conj(phi_bar)*g_in[:,:,:,0]
   return fe_es

def get_entropy_hermite(g_in,kzind=-1,include_kz0=True):
   """Gets entropy summed over kx,ky,kz but resolved in Hermite n. 
   kzind=-1 (default): sum over kz 
   kzind=kzind (default): selects kz index."""
   if include_kz0:
       ikz_start=0
   else:
       ikz_start=1

   if kzind==-1:
       entropy=2.0*np.real(0.5*np.pi*np.sum(np.sum(np.sum(np.conj(g_in[1:,:,ikz_start:,:])*g_in[1:,:,ikz_start:,:],axis=0),axis=0),axis=0))
       entropy=entropy+np.real(0.5*np.pi*np.sum(np.sum(np.conj(g_in[0,:,ikz_start:,:])*g_in[0,:,ikz_start:,:],axis=0),axis=0))
   else:
       g_kz=g_in[:,:,kzind,:]
       entropy=2.0*np.real(0.5*np.pi*np.sum(np.sum(np.conj(g_kz[1:,:,:])*g_kz[1:,:,:],axis=0),axis=0))
       entropy=entropy+np.real(0.5*np.pi*np.sum(np.conj(g_kz[0,:,:])*g_kz[0,:,:],axis=0))
   return entropy

def get_entropy_hermiten(g_in,kzind=-1,include_kz0=True):
   """Gets entropy summed over kx,ky,kz but resolved in Hermite n. 
   kzind=-1 (default): sum over kz 
   kzind=kzind (default): selects kz index."""
   if include_kz0:
       ikz_start=0
   else:
       ikz_start=1

   if kzind==-1:
       entropy=2.0*np.real(0.5*np.pi*np.sum(np.sum(np.sum(np.conj(g_in[1:,:,ikz_start:])*g_in[1:,:,ikz_start:],axis=0),axis=0),axis=0))
       entropy=entropy+np.real(0.5*np.pi*np.sum(np.sum(np.conj(g_in[0,:,ikz_start:])*g_in[0,:,ikz_start:],axis=0),axis=0))
   else:
       g_kz=g_in[:,:,kzind]
       entropy=2.0*np.real(0.5*np.pi*np.sum(np.sum(np.conj(g_kz[1:,:])*g_kz[1:,:],axis=0),axis=0))
       entropy=entropy+np.real(0.5*np.pi*np.sum(np.conj(g_kz[0,:])*g_kz[0,:],axis=0))
   return entropy

def get_rhs_lin(g_in,which_term):
    """Gets the part of the linear RHS operator determined by \'which_term\'.
    This routine is 4D.
    This should be benchmarked against the code occasionally.
    which_term:
    0:all
    1:collisions
    2:hyper collisions
    3:phase mixing
    4:omn term
    5:kz phi term
    6:omt term
    7:omt FLR term
    8:hyp_xyz 
    9:hyp_conv """

    #0==all
    #1==Collisions
    #2==hyper collisions
    #3==phase mixing
    #4==omn term
    #5==kz phi term
    #6==omt term
    #7==omt FLR term
    #8==hyps
    kx,ky,kz,herm=get_grids()
    #print kx
    #print ky
    #print kz
    J0a=get_gyroavg()
    rhs=np.zeros((par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),dtype='complex128')
    #Collisions
    if which_term==1 or which_term==0:
        for l in range(par['nv0']):
            if l > 2 or par['em_conserve'] != 'T':
                rhs[:,:,:,l]=rhs[:,:,:,l]-par['nu']*l*g_in[:,:,:,l]
    #Hyper collisions
    if which_term==2 or which_term==0:
        for l in range(par['nv0']):
            rhs[:,:,:,l]=rhs[:,:,:,l]-par['hyp_v']\
                    *(float(l)/float(par['nv0']-1))**\
                    par['hypv_order']*g_in[:,:,:,l]
    if which_term==3 or which_term==0:
        for k in range(par['nkz0']):
            for l in range(par['nv0']):
                if l != 0:
                    rhs[:,:,k,l]=rhs[:,:,k,l]+\
                        -1.0J*kz[k]*(herm[l])**0.5*g_in[:,:,k,l-1]
                if l != par['nv0']-1:
                    rhs[:,:,k,l]=rhs[:,:,k,l]+\
                        -1.0J*kz[k]*(herm[l+1])**0.5*g_in[:,:,k,l+1]
                elif par['nuno_closure']:
                    #print "closure"
                    g_np1=-1.0J*kz[k]*(par['nv0'])**0.5*g_in[:,:,k,par['nv0']-1]\
                            /(par['nu']*par['nv0']+par['hyp_v']*(float(par['nv0'])/float(par['nv0']-1))**\
                    par['hypv_order'])
                    rhs[:,:,k,l]=rhs[:,:,k,l]+\
                        -1.0J*kz[k]*(float(l+1))**0.5*g_np1
    if which_term==4 or which_term==0:
        phi=get_phi(g_in)
        for i in range(par['nkx0']):
            for j in range(par['nky0']):
                rhs[i,j,:,0]=rhs[i,j,:,0]+\
                        -1.0J*par['omn']*ky[j]*np.pi**(-0.25)*J0a[i,j]*phi[i,j,:]
    if which_term==5 or which_term==0:
        phi=get_phi(g_in)
        for i in range(par['nkx0']):
            for j in range(par['nky0']):
                for k in range(par['nkz0']):
                    rhs[i,j,k,1]=rhs[i,j,k,1]+\
                        -1.0J*kz[k]*np.pi**(-0.25)*J0a[i,j]*phi[i,j,k]
    if which_term==6 or which_term==0:
        phi=get_phi(g_in)
        for i in range(par['nkx0']):
            for j in range(par['nky0']):
                rhs[i,j,:,2]=rhs[i,j,:,2]+\
                        -2.0**(-0.5)*1.0J*ky[j]*np.pi**(-0.25)*J0a[i,j]\
                        *phi[i,j,:]*par['omt']
    if which_term==7 or which_term==0:
        phi=get_phi(g_in)
        for i in range(par['nkx0']):
            for j in range(par['nky0']):
                rhs[i,j,:,0]=rhs[i,j,:,0]+\
                        1.0J*ky[j]*np.pi**(-0.25)*0.5*(kx[i]**2+ky[j]**2)\
                        *J0a[i,j]*phi[i,j,:]*par['omt']
    if which_term==9 or which_term==0:
        if par['hyp_conv'] != 0.0:
            rhs[:,:,0,:]=rhs[:,:,0,:]-par['hyp_nu']*par['hyp_conv']*g_in[:,:,0,:]
            for nkhc in range(1,par['num_k_hyp_conv']+1):
                rhs[:,:,nkhc,:]=rhs[:,:,nkhc,:]-par['hyp_nu']*par['hyp_conv']*g_in[:,:,nkhc,:]
                rhs[:,:,-nkhc,:]=rhs[:,:,-nkhc,:]-par['hyp_nu']*par['hyp_conv']*g_in[:,:,-nkhc,:]
                if par['hyp_conv_ky']:
                    rhs[:,nkhc,:,:]=rhs[:,nkhc,:,:]-par['hyp_nu']*par['hyp_conv']*g_in[:,nkhc,:,:]
                    rhs[:,-nkhc,:,:]=rhs[:,-nkhc,:,:]-par['hyp_nu']*par['hyp_conv']*g_in[:,-nkhc,:,:]
    if which_term==8 or which_term==0:
        if par['hyp_x'] != 0.0 or par['hyp_y'] != 0.0:
            for i in range(par['nkx0']):
                for j in range(par['nky0']):
                    rhs[i,j,:,:]=rhs[i,j,:,:]-par['hyp_x']*g_in[i,j,:,:]*\
                        (kx[i]/par['kxmax0'])**par['hypx_order']\
                        -par['hyp_y']*g_in[i,j,:,:]*\
                        (ky[j]/par['kymax0'])**par['hypy_order']
        if par['hyp_z'] != 0.0 :
            for k in range(par['nkz0']):
                rhs[:,:,k,:]=rhs[:,:,k,:]-par['hyp_z']*g_in[:,:,k,:]*\
                        (kz[k]/par['kzmax0'])**par['hypz_order']

    bench=False
    if bench and which_term==0:
        print "Testing single k version (this might take a long time)"
        num_errors=0
        for i in range(par['nkx0']):
            for j in range(par['nky0']):
                for k in range(par['nkz0']):
                    rhsk=get_rhs_lin_single_k(g_in[i,j,k,:],0,\
                                            kx[i],ky[j],kz[k])
                    test=np.abs(np.sum(rhsk-rhs[i,j,k,:]))
                    if test/np.abs(np.sum(rhs[i,j,k,:])) > 1.0e-13:
                        print "Error, i,j,k",i,j,k,test/np.abs(np.sum(rhs[i,j,k,:]))
                        num_errors=num_errors+1
        print "Number of errors:",num_errors

    return rhs

def get_rhs_lin_single_k(g_in,which_term,kx_in,ky_in,kz_in,test_extra=False,zamp=1.0):
    """For the wavenumber defined by kx_in,ky_in,kz_in, this routine 
    gets the part of the linear RHS operator determined by \'which_term\'.
    This should be benchmarked against the code occasionally.
    which_term:
    0:all
    1:collisions
    2:hyper collisions
    3:phase mixing
    4:omn term
    5:kz phi term
    6:omt term
    7:omt FLR term
    7:hyps """
    #0==all
    #1==Collisions
    #2==hyper collisions
    #3==phase mixing
    #4==omn term
    #5==kz phi term
    #6==omt term
    #7==omt FLR term
    #8==hyps
    rhs=np.zeros((par['nv0']),dtype='complex')
    herm=np.arange(par['nv0'],dtype='float')

    #Collisions
    if which_term==1 or which_term==0:
        for l in range(par['nv0']):
            if l > 2 or par['em_conserve'] != 'T':
                rhs[l]=-par['nu']*l*g_in[l]
    #Hyper collisions
    if which_term==2 or which_term==0:
        for l in range(par['nv0']):
            rhs[l]=rhs[l]-par['hyp_v']\
                    *(float(l)/float(par['nv0']-1))**\
                    par['hypv_order']*g_in[l]
    if which_term==3 or which_term==0:
        for l in range(par['nv0']):
            if l != 0:
                rhs[l]=rhs[l]+\
                    -1.0J*kz_in*(herm[l])**0.5*g_in[l-1]
            if l != par['nv0']-1:
                rhs[l]=rhs[l]+\
                    -1.0J*kz_in*(herm[l+1])**0.5*g_in[l+1]
            elif par['nuno_closure'] and (par['nu']!=0.0 or par['hyp_v']!=0.0):
                #print "closure"
                g_np1=-1.0J*kz_in*(par['nv0'])**0.5*g_in[par['nv0']-1]\
                        /(par['nu']*par['nv0']+par['hyp_v']*(float(par['nv0'])/float(par['nv0']-1))**\
                        par['hypv_order'])
                rhs[l]=rhs[l]+\
                        -1.0J*kz_in*(par['nv0'])**0.5*g_np1
    if test_extra:   #Explore what creates zig-zag--roughly model nonlinearity
        for l in range(par['nv0']-1):
            #if l%2 == 0 :
            #    ll=float(l)
            #else:
            #    ll=float(l+1)
            #rhs[l]=rhs[l]+\
            #    -zamp*kz_in*(ll)**0.5*g_in[l]
            if l%2 == 0 :
                ll=float(l)
                rhs[l]=rhs[l]+\
                    -zamp*kz_in*(ll)**0.5*g_in[l]

    if which_term==4 or which_term==0:
        phi=get_phi_single_k(g_in,kx_in,ky_in)
        rhs[0]=rhs[0]+\
                -1.0J*ky_in*par['omn']*np.pi**(-0.25)*gyro_avg_k(kx_in,ky_in)*phi
    if which_term==5 or which_term==0:
        phi=get_phi_single_k(g_in,kx_in,ky_in)
        rhs[1]=rhs[1]+\
            -1.0J*kz_in*np.pi**(-0.25)*gyro_avg_k(kx_in,ky_in)*phi
    if which_term==6 or which_term==0:
        phi=get_phi_single_k(g_in,kx_in,ky_in)
        rhs[2]=rhs[2]+\
                -2.0**(-0.5)*1.0J*ky_in*np.pi**(-0.25)*gyro_avg_k(kx_in,ky_in)\
                *phi*par['omt']
    if which_term==7 or which_term==0:
        phi=get_phi_single_k(g_in,kx_in,ky_in)
        rhs[0]=rhs[0]+\
                1.0J*ky_in*np.pi**(-0.25)*0.5*(kx_in**2+ky_in**2)\
                *gyro_avg_k(kx_in,ky_in)*phi*par['omt']
    if which_term==8 or which_term==0:
        if par['hyp_conv'] != 0.0: 
                     #np.abs(kz_in-0.0) < 1.0e-8 or \
                     #np.abs(np.abs(kz_in)-par['kzmin']) < 1.0e-8 or \
                     #np.abs(np.abs(ky_in)-par['kymin']) < 1.0e-8):
            #print '!Warning :  check hyp_conv!!'
            if np.abs(kz_in)-par['kzmin']*par['num_k_hyp_conv'] < 1.0e-6:
                print 'hyp_conv', kx_in,ky_in,kz_in
                rhs[:]=rhs[:]-par['hyp_nu']*par['hyp_conv']*g_in[:]
            if np.abs(ky_in)-par['kymin']*par['num_k_hyp_conv'] < 1.0e-6 \
                   and np.abs(ky_in-0.0) > 1.0e-6 and par['hyp_conv_ky']:
                print 'hyp_conv', kx_in,ky_in,kz_in
                rhs[:]=rhs[:]-par['hyp_nu']*par['hyp_conv']*g_in[:]
        if par['hyp_x'] != 0.0 or par['hyp_y'] != 0.0:
            rhs[:]=rhs[:]-par['hyp_x']*g_in[:]*\
                        (kx_in/par['kxmax0'])**par['hypx_order']\
                        -par['hyp_y']*g_in[:]*\
                        (ky_in/par['kymax0'])**par['hypy_order']
        if par['hyp_z'] != 0.0:
            rhs[:]=rhs[:]-par['hyp_z']*g_in[:]*\
                        (kz_in/par['kzmax0'])**par['hypz_order']
    return rhs

def get_phi(g_in):
    """Returns the potential given the distribution function."""
    J0a=get_gyroavg()
    Gam0=get_gamma0()
    #print 'sum(J0a)',np.sum(J0a)
    #print 'sum(Gam0)',np.sum(Gam0)
    if par['etg_factor'] != 0.0:
        print "Diags only work for etg_factor=0.0 at this time!"
        stop
    phi=np.empty((par['nkx0'],par['nky0'],par['nkz0']),dtype='complex128')
    for i in range(par['nkx0']):
        for j in range(par['nky0']):
            #print "i,j in get_phi",i,j
            #print "J0",J0a[i,j]
            #print "g_in",g_in[i,j,:,0]
            #print "Gam0",Gam0[i,j]
            phi[i,j,:]=np.pi**(0.25)*J0a[i,j]*g_in[i,j,:,0]\
                    /(par['Ti0Te']+1.0-Gam0[i,j])
    #print 'sum(phi) in get_phi',np.sum(phi)
    return phi

def get_phi_bar(g_in):
    """Returns the gyro-averaged potential given the distribution function."""
    J0a=get_gyroavg()
    phi=get_phi(g_in)
    phi_bar=np.empty((par['nkx0'],par['nky0'],par['nkz0']),dtype='complex128')
    for i in range(par['nkx0']):
        for j in range(par['nky0']):
            phi_bar[i,j,:]=J0a[i,j]*phi[i,j,:]
    return phi_bar

def get_phi_single_k(g_in,kx_in,ky_in):
    """Returns the potential for kx_in, ky_in given the distribution function."""
    #phi=np.empty((par['nkx0'],par['nky0'],par['nkz0']),dtype='complex128')
    #for i in range(par['nkx0']):
    #    for j in range(par['nky0']):
    if par['etg_factor'] != 0.0:
        print "Diags only work for etg_factor=0.0 at this time!"
        stop
    phi=np.pi**(0.25)*gyro_avg_k(kx_in,ky_in)*g_in[0]\
            /(par['Ti0Te']+1.0-gamma0_k(kx_in,ky_in))
    return phi

def get_phi_prefactor(kx_in,ky_in):
    """Returns the prefactor necessary to change g_0 into phi."""
    if(ky_in==0.0):
        print "!!!!Not ready for ky=0!!!"
        stop
    else:
        prefactor=np.pi**(0.25)*np.e**(-0.5*(kx_in**2+ky_in**2))\
            /(par['Ti0Te']+1.0-gamma0_k(kx_in,ky_in))

    return prefactor

def get_gyroavg():
    """Returns the matrix for J0 as a function of kx,ky"""
    kx,ky,kz,herm=get_grids()
    J0a=np.empty((par['nkx0'],par['nky0']))
    for i in range(par['nkx0']):
        for j in range(par['nky0']):
            J0a[i,j]=np.e**(-0.5*(kx[i]**2+ky[j]**2))
    return J0a

def gyro_avg_k(kx_in,ky_in):
    """Returns J0 for a given kx_in, ky_in."""
    #kx,ky,kz,herm=get_grids()
    #J0a=np.empty((par['nkx0'],par['nky0']))
    J0=np.e**(-0.5*(kx_in**2+ky_in**2))
    return J0

def get_gamma0():
    """Returns Gamma0 as a function of kx,ky"""
    kx,ky,kz,herm=get_grids()
    Gam0=np.empty((par['nkx0'],par['nky0']))
    for i in range(par['nkx0']):
        for j in range(par['nky0']):
            Gam0[i,j]=sc.special.ive(0,kx[i]**2+ky[j]**2)
    return Gam0

def gamma0_k(kx_in,ky_in):
    """Returns Gamma0 for kx_in, ky_in"""
    Gam0=sc.special.ive(0,kx_in**2+ky_in**2)
    return Gam0

def get_kindex(k,kmin,nk0):
    """Given a k, this returns the appropriate index for accessing 
    this k in the FFT k-format (see get_grids())"""
    ik=np.rint(k/kmin)
    if k < 0.0:
        ik=nk0+ik
    return ik

def get_kxindex(k,kmin,nk0):
    """Given a kx, this returns the appropriate index for accessing 
    this kx in the FFT k-format (see get_grids())"""
    ik=np.rint(k/kmin)
    if k < 0.0:
        ik=nk0+ik
       
    return ik

def get_grids_shifted():
    """This returns grids in the -k==>k form. """
    kxgrid=np.arange((par['nkx0']))
    kxgrid=kxgrid*par['kxmin']

    kygrid=np.arange(par['nky0'])
    kygrid=kygrid*par['kymin']
    kygrid=kygrid-(par['ky_nyq']-par['kymin'])
    if np.min(abs(kygrid)) != 0.0:
        kygrid=kygrid-np.min(abs(kygrid))   

    kzgrid=np.arange(par['nkz0'])
    kzgrid=kzgrid*par['kzmin']
    kzgrid=kzgrid-(par['kz_nyq']-par['kzmin'])
    if np.min(abs(kzgrid)) != 0.0:
        kzgrid=kzgrid-np.min(abs(kzgrid))   

    herm_grid=np.arange(par['nv0'])
    herm_grid=1.0*herm_grid

    return kxgrid,kygrid,kzgrid,herm_grid

def get_shells(min_width=1.0):
    """Can't remember . . . \n
    but I think this returns the boundaries of the k-shells"""
    wid=min_width
    kx,ky,kz,herm=get_grids()  
    kmax=min(par['kxmax0'],par['kymax0'])
    dk0=kx[1]-kx[0]
    dky0=ky[1]-ky[0]
    kgrid=np.empty((par['nkx0']*par['nky0'],2))
    count=0
    for i in range(par['nkx0']):
        for j in range(par['nky0']):
            kgrid[count,0]=kx[i]
            kgrid[count,1]=ky[j]
            count=count+1

    #print kgrid
    #plt.scatter(kgrid[:,0],kgrid[:,1],s=1)
    #plt.show()
    if dk0 != dky0:
        print "Error! delta kx must equal delta ky"
        stop
    area=2.0*np.pi*(kmax-wid*dk0/2.0)*wid*dk0
    kbounds=np.empty(0)
    #kbounds=np.append(kbounds,kmax+0.5*wid*dk0)
    #print "kxmax+0.5*wid*dk0",kmax+0.5*wid*dk0
    #print "kbounds",kbounds
    angle=np.arange(200)/(199.0)*2.0*np.pi
    circlex=np.cos(angle)
    circley=np.sin(angle)

    next_bound=(area/np.pi)**0.5
    print "next boundary",next_bound
    kbounds=np.append(kbounds,next_bound)
    keep_going=1
    bcount=0
    while keep_going:
        next_bound=(kbounds[bcount]**2+area/np.pi)**0.5
        print "next boundary",next_bound
        if next_bound < kmax+dk0/2.0:
            kbounds=np.append(kbounds,next_bound)
            bcount=bcount+1
            plt.plot(next_bound*circlex,next_bound*circley)
        else:
            keep_going=0

    plt.scatter(kgrid[:,0],kgrid[:,1],s=1)
    prange=[-5,5,-5,5]
    plt.axis(prange)
    plt.show()

    num_shells=bcount+1
    print "Area per shell:"
    print kbounds[0],np.pi*kbounds[0]**2
    for i in range(1,num_shells):
        print kbounds[i],np.pi*(kbounds[i]**2-kbounds[i-1]**2)

    return kbounds

def get_time_from_eshells():
   """Get's time from eshells.dat"""
   f=open('shell_info.dat')
   s_info=f.read()
   s_lines=s_info.split('\n')
   num_shells=len(s_lines)-1#int(float(s_lines[-2].split()[0]))
   print num_shells
   f.close()

   #time=np.empty(0)
   #continue_read=1
   #i=0
   #for i in range(8*num_shells*par['nkz0']*par['nv0']):
   #    f.seek(i*8)
   #    i=i+1
   #    input=np.fromfile(f,dtype='float64',count=1)
   #    if input==0 or input:
   #        time = np.append(time,input)
   #        #print input
   #np.savetxt('timetest.dat',time)
   #stop

   f = open('eshells.dat','rb')
   ntot=num_shells*par['nkz0']*par['nv0']
   mem_tot=7*ntot*8+num_shells*par['nkz0']*8
   time=np.empty(0)
   continue_read=1
   i=0
   while (continue_read): 
     f.seek(i*(mem_tot+8))
     i=i+1
     input=np.fromfile(f,dtype='float64',count=1)
     if input==0 or input:
         time = np.append(time,input)
     else:
         continue_read=0

   f.close()
   return time

def eshells_test(start_time=-1,end_time=-1,nlt_comparison=False):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1#int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    energy=np.genfromtxt('energy.dat')
    eshells=np.zeros((ntime,8))
    ebal=np.zeros((ntime,4,num_shells))
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Flux (no v dependenc)

    es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    print "np.shape(es_in)"
    print np.shape(es_in)
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        eshells[i-istart,0]=np.sum(es_in)
        f.seek(8+i*(8+mem_tot)+ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        eshells[i-istart,1]=np.sum(es_in)
        ebal[i-istart,0,:]=distribute_ebal_shells(es_in)
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        eshells[i-istart,2]=np.sum(es_in)
        ebal[i-istart,1,:]=distribute_ebal_shells(es_in)
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        eshells[i-istart,3]=np.sum(es_in)
        ebal[i-istart,1,:]=ebal[i-istart,1,:]+distribute_ebal_shells(es_in)
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        eshells[i-istart,4]=np.sum(es_in)
        ebal[i-istart,3,:]=distribute_ebal_shells(es_in)
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        eshells[i-istart,5]=np.sum(es_in)
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        eshells[i-istart,6]=np.sum(es_in)
        f.seek(8+i*(8+mem_tot)+7*ntot*8)
        fs_in=np.zeros((par['nkz0'],num_shells))
        fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
        fs_in=np.reshape(fs_in,(par['nkz0'],num_shells))
        print "np.shape(fs_in)",np.shape(fs_in)
        eshells[i-istart,7]=np.sum(fs_in)
        for j in range(num_shells):
            ebal[i-istart,2,j]=np.sum(fs_in[:,j])
        f.close()
 
    plt.plot(time[istart:iend+1],eshells[:,0],'-x',label='Energy (from shells)')
    plt.plot(energy[:,0],energy[:,1]+energy[:,2],label='Energy')
    plt.legend()
    plt.show()
    plt.plot(time[istart:iend+1],eshells[:,1],'-x',label='dE/dt (from shells)')
    plt.plot(energy[:,0],energy[:,3],label='dE/dt')
    plt.legend()
    plt.show()
    plt.plot(time[istart:iend+1],eshells[:,2],'-x',label='coll (from shells)')
    plt.plot(energy[:,0],energy[:,5],label='coll')
    plt.legend()
    plt.show()
    plt.plot(time[istart:iend+1],eshells[:,3],'-x',label='hyps (from shells)')
    plt.plot(energy[:,0],energy[:,6]+energy[:,7]+energy[:,9],label='hyps')
    plt.legend()
    plt.show()
    plt.plot(time[istart:iend+1],eshells[:,4],'-x',label='NL (from shells)')
    plt.plot(energy[:,0],energy[:,8],label='NL')
    plt.legend()
    plt.show()
    plt.plot(time[istart:iend+1],eshells[:,5],label='PH1 (from shells)')
    plt.plot(time[istart:iend+1],eshells[:,6],label='PH2 (from shells)')
    plt.legend()
    plt.show()
    plt.plot(time[istart:iend+1],eshells[:,7],'-x',label='Flux (from shells)')
    plt.plot(energy[:,0],energy[:,4],label='Flux')
    plt.legend()
    plt.show()

    plt.title("Global Energetics",size=18)
    plt.plot(energy[:,0],energy[:,3],label='dE/dt')
    plt.plot(energy[:,0],energy[:,4],label='Flux')
    plt.plot(energy[:,0],energy[:,5],label='Coll+Diss')
    plt.plot(energy[:,0],energy[:,8],label='NL')
    plt.plot(energy[:,0],energy[:,6]+energy[:,7]+energy[:,9],label='hyps')
    plt.legend(loc='upper left')
    plt.xlabel(r'$t(R/v_{ti})$',size=18)
    plt.show()

    ntime=time[istart:iend+1]
    for i in range(num_shells):
        plt.title('Shell '+str(i+1))
        plt.plot(ntime,ebal[:,0,i],label='dE/dt')
        plt.plot(ntime,ebal[:,1,i],label='Coll+Diss')
        plt.plot(ntime,ebal[:,2,i],label='Flux')
        plt.plot(ntime,ebal[:,3,i],label='NL')
        plt.plot(ntime,np.sum(ebal[:,1:,i],axis=1),'x',label='Sum')
        plt.legend(loc='upper left')
        plt.xlabel(r'$t(R/v_{ti})$',size=18)
        plt.show()

    if nlt_comparison:
        nlt,nlt_time,nlt_shells_t=nlt_test()
        for i in range(num_shells):
            plt.plot(nlt_time,nlt_shells_t[i,:],'-x')
            plt.plot(ntime,ebal[:,3,i],'-x')
            plt.show()

def eshells_test_n(start_time=-1,end_time=-1,herm_n=0,nlt_comparison=False):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)


    ebal=np.zeros((ntime,6,num_shells))
    #order in ebal
    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2

    es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    print "np.shape(es_in)"
    print np.shape(es_in)
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        #RHS
        f.seek(8+i*(8+mem_tot)+ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,0,:]=distribute_ebal_shells_n(es_in,herm_n)
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,1,:]=distribute_ebal_shells_n(es_in,herm_n)
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,1,:]=ebal[i-istart,1,:]+distribute_ebal_shells_n(es_in,herm_n)
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,3,:]=distribute_ebal_shells_n(es_in,herm_n)
        #PH1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,4,:]=distribute_ebal_shells_n(es_in,herm_n)
        #PH2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,5,:]=distribute_ebal_shells_n(es_in,herm_n)
        if herm_n==2:
            f.seek(8+i*(8+mem_tot)+7*ntot*8)
            #Drive
            fs_in=np.zeros((par['nkz0'],num_shells))
            fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
            fs_in=np.reshape(fs_in,(par['nkz0'],num_shells))
            #print "np.shape(fs_in)",np.shape(fs_in)
            for j in range(num_shells):
                ebal[i-istart,2,j]=np.sum(fs_in[:,j])
        else:
            ebal[i-istart,2,:]=0.0
        f.close()

    ntime=time[istart:iend+1]
    for i in range(num_shells):
        plt.plot(ntime,ebal[:,0,i],label='dE/dt')
        plt.plot(ntime,ebal[:,1,i],label='Coll+Diss')
        plt.plot(ntime,ebal[:,2,i],label='Flux')
        plt.plot(ntime,ebal[:,3,i],label='NL')
        plt.plot(ntime,ebal[:,4,i],label='PM n-1')
        plt.plot(ntime,ebal[:,5,i],label='PM n+1')
        plt.plot(ntime,np.sum(ebal[:,1:,i],axis=1),'+',color='black',label='Sum')
        plt.title("Energetics: Hermite n = "+str(herm_n)+", Shell= "+str(i))
        #plt.title("Hermite: "+str(herm_n))
        #if herm_n==0 or herm_n==1:
        #    plt.annotate("Note: discrepancies for n=0 and n=1 cancel.",[0,0])
        plt.legend(loc='upper left')
        plt.xlabel(r'$t(R/v_{ti})$',size=18)
        plt.show()
        if herm_n==0 or herm_n==1 or herm_n==2:
            plt.plot(ntime,np.sum(ebal[:,1:,i],axis=1)-ebal[:,0,i],':x',label='Sum-dEdt')
            plt.legend()
            plt.title("Discrepancy")
            plt.show()

    if nlt_comparison:
        if herm_n==0:
            fname='nlt_shells_n0.dat'
        elif herm_n==1:
            fname='nlt_shells_n1.dat'
        elif herm_n==2:
            fname='nlt_shells_n2.dat'
        elif herm_n==3:
            fname='nlt_shells_n3.dat'
        elif herm_n==par['nv0']/8:
            fname='nlt_shells_n1o8.dat'
        elif herm_n==2*par['nv0']/8:
            fname='nlt_shells_n2o8.dat'
        elif herm_n==3*par['nv0']/8:
            fname='nlt_shells_n3o8.dat'
        elif herm_n==4*par['nv0']/8:
            fname='nlt_shells_n4o8.dat'
        elif herm_n==5*par['nv0']/8:
            fname='nlt_shells_n5o8.dat'
        elif herm_n==6*par['nv0']/8:
            fname='nlt_shells_n6o8.dat'
        elif herm_n==7*par['nv0']/8:
            fname='nlt_shells_n7o8.dat'
        else:
            print "Invalid herm_n for nlt comparison."
            stop

        nlt,nlt_time,nlt_shells_t=nlt_test(file_name=fname)
        for i in range(num_shells):
            plt.plot(nlt_time,nlt_shells_t[i,:],'-x')
            plt.plot(ntime,ebal[:,3,i],'-x')
            plt.show()

def eshells_test_n_kz(start_time=-1,end_time=-1,herm_n=0,kz=0.1,nlt_comparison=False):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    kzind=np.argmin(abs(kzgrid-kz))

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)


    ebal=np.zeros((ntime,7,num_shells))
    #order in ebal
    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Energy

    shells_in=np.genfromtxt('shell_info.dat')

    gamma=np.zeros(num_shells)
    for i in range(num_shells):
        if i==0:
            kperp=shells_in[0]/2.0
        else:
            kperp=shells_in[i-1]
        mat=get_lin_matrix(0.0,kperp,kz)
        om,gm,ev=get_ev_spectrum(mat)
        gamma[i]=max(gm)
        print "gamma",gamma[i]

    es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    print "np.shape(es_in)"
    print np.shape(es_in)
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,6,:]=distribute_ebal_shells_n_kz(es_in,herm_n,kz)
        #RHS
        f.seek(8+i*(8+mem_tot)+ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,0,:]=distribute_ebal_shells_n_kz(es_in,herm_n,kz)
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,1,:]=distribute_ebal_shells_n_kz(es_in,herm_n,kz)
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,1,:]=ebal[i-istart,1,:]+distribute_ebal_shells_n_kz(es_in,herm_n,kz)
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,3,:]=distribute_ebal_shells_n_kz(es_in,herm_n,kz)
        #PH1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,4,:]=distribute_ebal_shells_n_kz(es_in,herm_n,kz)
        #PH2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[i-istart,5,:]=distribute_ebal_shells_n_kz(es_in,herm_n,kz)
        if herm_n==2:
            f.seek(8+i*(8+mem_tot)+7*ntot*8)
            #Drive
            fs_in=np.zeros((par['nkz0'],num_shells))
            fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
            fs_in=np.reshape(fs_in,(par['nkz0'],num_shells))
            #print "np.shape(fs_in)",np.shape(fs_in)
            for j in range(num_shells):
                ebal[i-istart,2,j]=(fs_in[kzind,j])
        else:
            ebal[i-istart,2,:]=0.0
        f.close()

    ntime=time[istart:iend+1]
    gam_line=np.zeros(len(ntime))
    for i in range(num_shells):
        #Test dEdt from file vs numerical calc
        #temp=(ebal[:,6,i]-np.roll(ebal[:,6,i],1))
        #temp[0]=0.0
        #temp_dt=ntime-np.roll(ntime,1)
        #temp_dt[0]=0.0
        #dt=np.sum(temp_dt)/float((len(ntime)-1))
        #temp=temp/dt
        #plt.plot(ntime,ebal[:,0,i],label='dE/dt')
        #plt.plot(ntime,temp,label='dE/dt_numerical')
        #plt.show()

        plt.plot(ntime,ebal[:,0,i],label='dE/dt')
        plt.plot(ntime,ebal[:,1,i],label='Coll+Diss')
        plt.plot(ntime,ebal[:,2,i],label='Flux')
        plt.plot(ntime,ebal[:,3,i],label='NL')
        if herm_n==0 and not 'version_flag' in par.keys():
            Jt0=-1.0*np.sum(ebal[:,1:6,i],axis=1)+ebal[:,0,i]
            ebal[:,5,i]=ebal[:,5,i]+Jt0  #This term (J_phi) isn't included in eshells.dat
            plt.annotate("Note: early version (<54a) fix.",[0,0])
        elif herm_n==1 and not 'version_flag' in par.keys():
            Jt0=-1.0*np.sum(ebal[:,1:6,i],axis=1)+ebal[:,0,i]
            ebal[:,4,i]=ebal[:,4,i]+Jt0  #This term (J_phi) isn't included in eshells.dat
            plt.annotate("Note: early version (<54a) fix.",[0,0])
        plt.plot(ntime,ebal[:,4,i],label='PM n-1')
        plt.plot(ntime,ebal[:,5,i],label='PM n+1')
        plt.plot(ntime,np.sum(ebal[:,1:6,i],axis=1),'+',color='black',label='Sum')
        plt.title("Energetics: Hermite n = "+str(herm_n)+ r"$ k_z R=$"+str(kz)+"Shell= "+str(i))
        #plt.title("Hermite: "+str(herm_n))
        #if herm_n==0 or herm_n==1:
        #    plt.annotate("Note: discrepancies for n=0 and n=1 cancel.",[0,0])
        plt.legend(loc='upper left')
        plt.xlabel(r'$t(R/v_{ti})$',size=18)
        plt.show()
        #if herm_n==0 or herm_n==1:
        #    plt.plot(ntime,np.sum(ebal[:,1:6,i],axis=1)-ebal[:,0,i],':x',label='Sum-dEdt')
        #    plt.legend()
        #    plt.title("Discrepancy")
        #    plt.show()

        gam_line[:]=gamma[i]
        plt.plot(ntime,ebal[:,0,i]/ebal[:,6,i]/2.0,label='GR:dE/dt')
        plt.plot(ntime,ebal[:,1,i]/ebal[:,6,i]/2.0,label='GR:Coll+Diss')
        plt.plot(ntime,ebal[:,2,i]/ebal[:,6,i]/2.0,label='GR:Flux')
        plt.plot(ntime,ebal[:,3,i]/ebal[:,6,i]/2.0,label='GR:NL')
        plt.plot(ntime,ebal[:,4,i]/ebal[:,6,i]/2.0,label='GR:PM n-1')
        plt.plot(ntime,ebal[:,5,i]/ebal[:,6,i]/2.0,label='GR:PM n+1')
        plt.plot(ntime,np.sum(ebal[:,1:6,i],axis=1)/ebal[:,6,i]/2.0,'+',color='black',label='Sum')
        plt.plot(ntime,gam_line,'--',label='GR:linear',color='black')
        plt.title("Growth Rates: Hermite n = "+str(herm_n)+ r"$ k_z R=$"+str(kz)+"Shell= "+str(i))
        #plt.title("Hermite: "+str(herm_n))
        #if herm_n==0 or herm_n==1:
        #    plt.annotate("Note: discrepancies for n=0 and n=1 cancel.",[0,0])
        plt.legend(loc='upper left')
        plt.xlabel(r'$t(R/v_{ti})$',size=18)
        plt.show()
        #if herm_n==0 or herm_n==1:
        #    plt.plot(ntime,(np.sum(ebal[:,1:6,i],axis=1)-ebal[:,0,i])/ebal[:,6,i],':x',label='Sum-dEdt')
        #    plt.legend()
        #    plt.title("Discrepancy")
        #    plt.show()



    if nlt_comparison:
        if herm_n==0:
            fname='nlt_shells_n0.dat'
        elif herm_n==1:
            fname='nlt_shells_n1.dat'
        elif herm_n==2:
            fname='nlt_shells_n2.dat'
        elif herm_n==3:
            fname='nlt_shells_n3.dat'
        elif herm_n==par['nv0']/4:
            fname='nlt_shells_n1o4.dat'
        elif herm_n==par['nv0']/2:
            fname='nlt_shells_n2o4.dat'
        elif herm_n==3*par['nv0']/4:
            fname='nlt_shells_n3o4.dat'
        else:
            print "Invalid herm_n for nlt comparison."
            stop

        nlt,nlt_time,nlt_shells_t=nlt_test(file_name=fname)
        for i in range(num_shells):
            plt.plot(nlt_time,nlt_shells_t[i,:],'-x')
            plt.plot(ntime,ebal[:,3,i],'-x')
            plt.show()

def eshells_cb_kzn(start_time=-1,end_time=-1,which_shell=-1,eshells_test=False,start_shell=0,end_shell=-1):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
    if end_shell==-1:
        end_shell=num_shells
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)

    ebal=np.zeros((par['nv0'],par['nkz0'],num_shells,7))
    #order in ebal
    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy

    es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    print "np.shape(es_in)"
    print np.shape(es_in)
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,6]+=es_in
        #RHS
        f.seek(8+i*(8+mem_tot)+ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,0]+=es_in
        #ebal[:,:,0]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,1]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,1]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,3]+=es_in
        #ebal[:,:,3]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,4]+=es_in
        #ebal[:,:,4]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,5]+=es_in
        #ebal[:,:,5]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #if herm_n==2:
        f.seek(8+i*(8+mem_tot)+7*ntot*8)
        #Drive
        fs_in=np.zeros((par['nkz0'],num_shells))
        fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
        fs_in=np.reshape(fs_in,(par['nkz0'],num_shells))
        for j in range(num_shells):
            ebal[2,:,:,2]+=fs_in[:,:]
        f.close()

    ebal=ebal/float(ntime)
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    ebal_kzsum=np.sum(ebal,axis=1)    
    ebal_kpkzsum=np.sum(np.sum(ebal,axis=1),axis=1)    

    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy

    n0p5=herm_grid**(-0.5)
    n1p0=herm_grid**(-1.0)
    n1p5=herm_grid**(-1.5)

    ebal_akz=np.zeros((par['nv0'],par['nkz0']/2-1,num_shells,7),dtype='float')
    ebal_akz[:,0,:,:]=ebal[:,0,:,:]
    for k in range(1,par['nkz0']/2-1):
        ebal_akz[:,k,:,:]=ebal[:,k,:,:]+ebal[:,par['nkz0']-k,:,:]

    ebal_akz_kpsum=np.sum(ebal_akz[:,:,start_shell:end_shell,:],axis=2)

    if eshells_test:
        #ebal_kpkzsum=np.sum(np.sum(ebal_akz,axis=1),axis=1)
        ebal_kpkzsum=np.sum(ebal_akz_kpsum,axis=1)
        plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,1]+ebal_kpkzsum[:,3]+ebal_kpkzsum[:,4]+ebal_kpkzsum[:,5]),label='sum')
        plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,0]),label='RHS')
        plt.legend()
        plt.show()
        for k in range(par['nkz0']/2-1):
            plt.title('ikz '+str(k))
            plt.loglog(herm_grid,np.abs(ebal_akz_kpsum[:,k,1]+ebal_akz_kpsum[:,k,3]+ebal_akz_kpsum[:,k,4]+ebal_akz_kpsum[:,k,5]),label='sum')
            plt.loglog(herm_grid,np.abs(ebal_akz_kpsum[:,k,0]),label='RHS')
            plt.legend()
            plt.show()
            

    for k in range(par['nkz0']/2-1):
        plt.title('Free Energy'+str(k))
        if k != 0:
           n1p5=2.0*np.abs(ebal_akz_kpsum[15,k,6])/n1p5[15]*n1p5
           n0p5=2.0*np.abs(ebal_akz_kpsum[15,k,6])/n0p5[15]*n0p5
           n1p0=2.0*np.abs(ebal_akz_kpsum[15,k,6])/n1p0[15]*n1p0
        plt.loglog(herm_grid,np.abs(ebal_akz_kpsum[:,k,6]),'-',label='FE',basex=10,basey=10)
        temp=np.abs(ebal_akz_kpsum[:,k,4]+ebal_akz_kpsum[:,k,5])
        temp=ebal_akz_kpsum[15,k,6]/temp[15]*temp
        plt.loglog(herm_grid,temp,'x',label='scaled PM1+PM2',basex=10,basey=10)
        plt.loglog(herm_grid,n0p5,'--',label='-0.5',basex=10,basey=10)
        plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
        plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
        plt.legend()
        plt.show()

        temp[:]=0.0
        for n in range(par['nv0']-1):
            temp[n]=np.abs(kzgrid[k])*(np.sqrt(n)*ebal_akz_kpsum[n,k,6]-np.sqrt(n+1)*ebal_akz_kpsum[n+1,k,6])
        plt.title('|kz|='+str(par['kzmin']*k))
        plt.loglog(herm_grid,np.abs(ebal_akz_kpsum[:,k,3]),'-',label='NL',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(ebal_akz_kpsum[:,k,1]),'-',label='Diss',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(temp),'-',label='EqVI',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(ebal_akz_kpsum[:,k,4]+ebal_akz_kpsum[:,k,5]),'x',label='PM1+PM2',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(ebal_akz_kpsum[:,k,4]+ebal_akz_kpsum[:,k,5]+ebal_akz_kpsum[:,k,1]+ebal_akz_kpsum[:,k,3]),':',label='Sum',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(ebal_akz_kpsum[:,k,0]),'-.',label='RHS',basex=10,basey=10)
        if k != 0:
           n1p5=2.0*np.abs(ebal_akz_kpsum[15,k,4]+ebal_akz_kpsum[15,k,5])/n1p5[15]*n1p5
           n0p5=2.0*np.abs(ebal_akz_kpsum[15,k,4]+ebal_akz_kpsum[15,k,5])/n0p5[15]*n0p5
           n1p0=2.0*np.abs(ebal_akz_kpsum[15,k,4]+ebal_akz_kpsum[15,k,5])/n1p0[15]*n1p0
        plt.loglog(herm_grid,n0p5,'--',label='-0.5',basex=10,basey=10)
        plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
        plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
        plt.legend()
        plt.show()
    
    if which_shell != -1:
        for k in range(par['nkz0']/2-1):
            plt.title('|kz|='+str(par['kzmin']*k))
            plt.loglog(herm_grid,np.abs(ebal_akz[:,k,which_shell,3]),'-',label='NL shell '+str(which_shell),basex=10,basey=10)
            plt.loglog(herm_grid,np.abs(ebal_akz[:,k,which_shell,1]),'-',label='Diss shell '+str(which_shell),basex=10,basey=10)
            plt.loglog(herm_grid,np.abs(ebal_akz[:,k,which_shell,4]+ebal_akz[:,k,which_shell,5]),'x',label='PM1+PM2 shell '+str(which_shell),basex=10,basey=10)
            if k != 0:
               n1p5=2.0*np.abs(ebal_akz[2,k,which_shell,4]+ebal_akz[2,k,which_shell,5])/n1p5[2]*n1p5
               n0p5=2.0*np.abs(ebal_akz[2,k,which_shell,4]+ebal_akz[2,k,which_shell,5])/n0p5[2]*n0p5
               n1p0=2.0*np.abs(ebal_akz[2,k,which_shell,4]+ebal_akz[2,k,which_shell,5])/n1p0[2]*n1p0
            print "n1p5",n1p5
            plt.loglog(herm_grid,n0p5,'--',label='-0.5',basex=10,basey=10)
            plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
            plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
            plt.legend()
            plt.show()

def eshells_energy(start_time=-1,end_time=-1,show_plots=False):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)

    ebal=np.zeros((par['nv0'],par['nkz0'],num_shells,8))
    #order in ebal
    #0:RHS
    #1:coll 
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy
    #7:hyps

    es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    print "np.shape(es_in)"
    print np.shape(es_in)
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,6]+=es_in
        #RHS
        f.seek(8+i*(8+mem_tot)+ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,0]+=es_in
        #ebal[:,:,0]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,1]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,7]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,3]+=es_in
        #ebal[:,:,3]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,4]+=es_in
        #ebal[:,:,4]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,5]+=es_in
        #ebal[:,:,5]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #if herm_n==2:
        f.seek(8+i*(8+mem_tot)+7*ntot*8)
        #Drive
        fs_in=np.zeros((par['nkz0'],num_shells))
        fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
        fs_in=np.reshape(fs_in,(par['nkz0'],num_shells))
        for j in range(num_shells):
            ebal[2,:,:,2]+=fs_in[:,:]
        f.close()

    ebal=ebal/float(ntime)
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    ebal_kzsum=np.sum(ebal,axis=1)    
    ebal_kpkzsum=np.sum(np.sum(ebal[:,:,:,:],axis=1),axis=1)    
    ebal_kpsum=np.sum(ebal[:,:,:,:],axis=2)

    print "Enter the shell where hyp_x,y becomes large (identify from plot):\n"
    plt.plot(np.sum(ebal_kzsum,axis=0)[:,7],label='hyps')
    plt.legend()
    plt.show()
    hyp_shell=int(float(raw_input('Enter shell:')))

    ebal_kpkzsum_subhyp=np.sum(np.sum(ebal[:,:,0:hyp_shell,:],axis=1),axis=1)    
    ebal_kpsum_subhyp=np.sum(ebal[:,:,0:hyp_shell,:],axis=2)

    #Jptest=np.empty((par['nv0'],par['nkz0'],num_shells))
    #for k in range(par['nkz0']):
    #    for n in range(par['nv0']):
    #        Jptest[n,k,:]=2.0*np.sqrt(n)*kzgrid[k]*ebal[n,k,:,6]
    #ebal_kzsum=np.sum(ebal,axis=1)    
    #Jptest_kzsum=np.sum(Jptest,axis=1)    
    #ebal_kpkzsum=np.sum(np.sum(ebal,axis=1),axis=1)    
    #Jptest_kpkzsum=np.sum(np.sum(Jptest,axis=1),axis=1)    

    #print np.info(ebal_kpkzsum)
    #print ebal_kpkzsum[:,6]
    #print herm_grid
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,6]),':',label='Total FE',basex=10,basey=10)
    np.savetxt(par['diagdir'][1:-1]+'/ese_hspect.dat',ebal_kpkzsum[:,6])
    n1p5=herm_grid**(-1.5)
    n1p5=2.0*np.abs(ebal_kpkzsum[20,6])/n1p5[20]*n1p5
    n1p0=herm_grid**(-1.0)
    n1p0=5.0*np.abs(ebal_kpkzsum[20,6])/n1p0[20]*n1p0
    plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
    plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
    plt.legend()
    plt.show()

    plt.loglog(herm_grid,np.abs(ebal_kpkzsum_subhyp[:,6]),':',label='FE (sub-hyp)',basex=10,basey=10)
    np.savetxt(par['diagdir'][1:-1]+'/ese_hspect_subhyp.dat',ebal_kpkzsum_subhyp[:,6])
    np.savetxt(par['diagdir'][1:-1]+'/ese_cspect_subhyp.dat',ebal_kpkzsum_subhyp[:,1])
    np.savetxt(par['diagdir'][1:-1]+'/ese_hypspect_subhyp.dat',ebal_kpkzsum_subhyp[:,7])
    np.savetxt(par['diagdir'][1:-1]+'/ese_nlspect_subhyp.dat',ebal_kpkzsum_subhyp[:,3])
    n1p5=herm_grid**(-1.5)
    n1p5=2.0*np.abs(ebal_kpkzsum[20,6])/n1p5[20]*n1p5
    n1p0=herm_grid**(-1.0)
    n1p0=5.0*np.abs(ebal_kpkzsum[20,6])/n1p0[20]*n1p0
    plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
    plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
    plt.legend()
    plt.show()

    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy
    for i in range(num_shells):
        plt.title('shell '+str(i))
        plt.loglog(herm_grid,np.abs(ebal_kzsum[:,i,6]),label='FE',basex=10,basey=10)
        np.savetxt(par['diagdir'][1:-1]+'/ese_hspect_shell'+str(i)+'.dat',ebal_kzsum[:,i,6])
        np.savetxt(par['diagdir'][1:-1]+'/ese_cspect_shell'+str(i)+'.dat',ebal_kzsum[:,i,1])
        np.savetxt(par['diagdir'][1:-1]+'/ese_hypspect_shell'+str(i)+'.dat',ebal_kzsum[:,i,7])
        np.savetxt(par['diagdir'][1:-1]+'/ese_nlspect_shell'+str(i)+'.dat',ebal_kzsum[:,i,3])
        n1p5=5.0*np.abs(ebal_kzsum[20,i,6])/n1p5[20]*n1p5
        n1p0=5.0*np.abs(ebal_kzsum[20,i,6])/n1p0[20]*n1p0
        plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
        plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
        plt.legend()
        plt.show()

    kzn=np.empty((par['nv0']))
    kzn_subhyp=np.empty((par['nv0']))
    #plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,6]),':',label='FE',basex=10,basey=10)
    for n in range(par['nv0']):
        kzn[n]=np.sum(np.abs(kzgrid)*ebal_kpsum[n,:,6])/np.sum(ebal_kpsum[n,:,6])
        kzn_subhyp[n]=np.sum(np.abs(kzgrid)*ebal_kpsum_subhyp[n,:,6])/np.sum(ebal_kpsum_subhyp[n,:,6])
        if show_plots:
            plt.plot(kzgrid,ebal_kpsum[n,:,6],'x-',label='Free Energy')
            ax=plt.axis()
            plt.vlines(kzn[n],ax[2],ax[3])
            plt.title("Hermite n = "+str(n)+", Shell= "+str(i))
            plt.legend(loc='upper right')
            plt.xlabel(r'$k_z L$',size=18)
            plt.show()
    plt.plot(np.arange(par['nv0']),kzn[:])
    np.savetxt(par['diagdir'][1:-1]+'/ese_kzn.dat',kzn)
    plt.legend() 
    plt.title(r'$<k_z>_n$')
    plt.show()
    plt.plot(np.arange(par['nv0']),kzn[:])
    np.savetxt(par['diagdir'][1:-1]+'/ese_kzn_subhyp.dat',kzn_subhyp)
    plt.legend() 
    plt.title(r'$<k_z>_n$')
    plt.show()

    plt.plot(np.arange(par['nv0']),kzn[:])
    plt.plot(np.arange(par['nv0']),kzn[:]/np.sqrt(herm_grid))
    plt.legend() 
    plt.title(r'$<k_z>_n/\sqrt{n}$')
    plt.show()

    kzn_shell=np.empty((par['nv0'],num_shells))
    for i in range(num_shells):
        np.savetxt(par['diagdir'][1:-1]+'/ese_kzspect_shell'+str(i)+'.dat',np.transpose(ebal[:,:,i,6]))
        for n in range(par['nv0']):
            kzn_shell[n,i]=np.sum(np.abs(kzgrid)*ebal[n,:,i,6])/np.sum(ebal[n,:,i,6])
            #print kz0
            if show_plots:
                plt.plot(kzgrid,ebal[n,:,i,6],'x-',label='Free Energy')
                ax=plt.axis()
                plt.vlines(kzn_shell[n,i],ax[2],ax[3])
                plt.title("Hermite n = "+str(n)+", Shell= "+str(i))
                plt.legend(loc='upper right')
                plt.xlabel(r'$k_z L$',size=18)
                plt.show()
        plt.plot(np.arange(par['nv0']),kzn_shell[:,i],label=str(i))
        np.savetxt(par['diagdir'][1:-1]+'/ese_kzn_shell'+str(i)+'.dat',kzn_shell[:,i])
    plt.legend() 
    plt.title(r'$<k_z>_n$')
    plt.show()

    for i in range(num_shells):
        plt.plot(np.arange(par['nv0']),kzn_shell[:,i]/np.sqrt(herm_grid),label=str(i))
    plt.legend() 
    plt.title(r'$<k_z>_n/\sqrt{n}$')
    plt.show()


def eshells_nlflux(start_time=-1,end_time=-1):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)

    ebal=np.zeros((par['nv0'],par['nkz0'],num_shells,8))
    #order in ebal
    #0:RHS
    #1:coll 
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy
    #7:hyps

    es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,6]+=es_in
        #RHS
        f.seek(8+i*(8+mem_tot)+ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,0]+=es_in
        #ebal[:,:,0]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,1]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,7]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,3]+=es_in
        #ebal[:,:,3]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,4]+=es_in
        #ebal[:,:,4]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,5]+=es_in
        #ebal[:,:,5]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #if herm_n==2:
        f.seek(8+i*(8+mem_tot)+7*ntot*8)
        #Drive
        fs_in=np.zeros((par['nkz0'],num_shells))
        fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
        fs_in=np.reshape(fs_in,(par['nkz0'],num_shells))
        for j in range(num_shells):
            ebal[2,:,:,2]+=fs_in[:,:]
        f.close()

    ebal=ebal/float(ntime)
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    ebal_kzsum=np.sum(ebal,axis=1)    
    ebal_kpkzsum=np.sum(np.sum(ebal[:,:,:,:],axis=1),axis=1)    
    ebal_kpsum=np.sum(ebal[:,:,:,:],axis=2)

    print "Shell 0: 0.0", s_lines[1]
    for i in range(1,num_shells-1):
        print "Shell ",i,":",s_lines[i],s_lines[i+1]
    print "Shell 0:", s_lines[-1], "max"

    print "Enter the desired flux shell large (identify from plot):\n"
    print "(Note: sum will include all shells larger than selection)\n"
    plt.plot(np.sum(ebal_kzsum,axis=0)[:,7],label='hyps')
    plt.legend()
    plt.show()
    flux_shell=int(float(raw_input('Enter shell:')))

    nlflux_shell0=np.sum(np.sum(ebal[:,:,(flux_shell+1):,3],axis=1),axis=1)    
    energy_shell0=np.sum(ebal[:,:,flux_shell,6],axis=1)
    #ebal_kpsum_subhyp=np.sum(ebal[:,:,0:hyp_shell,:],axis=2)
    cdiss_ratio=np.sum(ebal[:,:,(flux_shell+1):,1])/np.sum(ebal[:,:,:,1])
    print "collisional dissipation (shell+) over total",cdiss_ratio
    diss_highk=np.sum(ebal[:,:,(flux_shell+1):,1]+ebal[:,:,(flux_shell+1):,7])
    print "Dissipation (shell+):",diss_highk
    diss_total=np.sum(ebal[:,:,:,1]+ebal[:,:,:,7])
    print "Dissipation (shell-):",diss_total-diss_highk
    print "Dissipation (total):",diss_total
    plt.plot(nlflux_shell0,label='flux shell')
    plt.legend()
    plt.show()
    plt.plot(nlflux_shell0/energy_shell0,label='omega nl')
    plt.legend()
    plt.show()

    print "Total omega_nl:", np.sum(nlflux_shell0)/np.sum(energy_shell0)







#def eshells_cb_test(start_time=-1,end_time=-1,show_plots=False,test_eshells=False,end_shell=-1):
#    f=open('shell_info.dat')
#    s_info=f.read()
#    s_lines=s_info.split('\n')
#    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
#    print num_shells
#    f.close()
#    if end_shell==-1:
#        end_shell=num_shells
# 
#    f = open('eshells.dat','rb')
#    ntot=num_shells*par['nkz0']*par['nv0']
#    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
#    time=np.empty(0)
#    continue_read=1
#    i=0
#    while (continue_read): 
#      f.seek(i*(mem_tot+8))
#      i=i+1
#      input=np.fromfile(f,dtype='float64',count=1)
#      if input==0 or input:
#          time = np.append(time,input)
#      else:
#          continue_read=0
 






def eshells_cb_test(start_time=-1,end_time=-1,show_plots=False,test_eshells=False,end_shell=-1,plot_hs_shells=True):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
    if end_shell==-1:
        end_shell=num_shells
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)

    ebal=np.zeros((par['nv0'],par['nkz0'],num_shells,7))
    #order in ebal
    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy

    es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    print "np.shape(es_in)"
    print np.shape(es_in)
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,6]+=es_in
        #RHS
        f.seek(8+i*(8+mem_tot)+ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,0]+=es_in
        #ebal[:,:,0]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,1]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,1]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,3]+=es_in
        #ebal[:,:,3]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,4]+=es_in
        #ebal[:,:,4]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,5]+=es_in
        #ebal[:,:,5]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #if herm_n==2:
        f.seek(8+i*(8+mem_tot)+7*ntot*8)
        #Drive
        fs_in=np.zeros((par['nkz0'],num_shells))
        fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
        fs_in=np.reshape(fs_in,(par['nkz0'],num_shells))
        for j in range(num_shells):
            ebal[2,:,:,2]+=fs_in[:,:]
        f.close()

    ebal=ebal/float(ntime)
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    ebal_kzsum=np.sum(ebal,axis=1)    
    ebal_kpkzsum=np.sum(np.sum(ebal[:,:,0:end_shell],axis=1),axis=1)    

    #Jptest=np.empty((par['nv0'],par['nkz0'],num_shells))
    #for k in range(par['nkz0']):
    #    for n in range(par['nv0']):
    #        Jptest[n,k,:]=2.0*np.sqrt(n)*kzgrid[k]*ebal[n,k,:,6]
    #ebal_kzsum=np.sum(ebal,axis=1)    
    #Jptest_kzsum=np.sum(Jptest,axis=1)    
    #ebal_kpkzsum=np.sum(np.sum(ebal,axis=1),axis=1)    
    #Jptest_kpkzsum=np.sum(np.sum(Jptest,axis=1),axis=1)    

    if end_shell != -1:
        plt.title('end shell = '+str(end_shell)+' of '+str(num_shells))
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,0]),':',label='RHS',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,3]),'-',label='NL',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,1]),'-',label='Diss',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,4]+ebal_kpkzsum[:,5]),'x-',label='PM1+PM2',basex=10,basey=10)
    n1p5=herm_grid**(-1.5)
    n1p5=2.0*np.abs(ebal_kpkzsum[20,4]+ebal_kpkzsum[20,5])/n1p5[20]*n1p5
    plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,4]\
                      +ebal_kpkzsum[:,5]\
                      +ebal_kpkzsum[:,3]\
                      +ebal_kpkzsum[:,1]),'-.'\
                      ,label='Sum')
    plt.legend()
    plt.show()


    plt.figure(figsize=(4.5,3.5))
    fig=plt.gcf()
    fig.subplots_adjust(bottom=0.15)
    #plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,0]),':',label='RHS',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,6]),'-',color='blue',label=r'$\varepsilon_n$',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,3]),'-.',color='red',label=r'$N_{{\bf k},n}^{(f)}$',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,4]+ebal_kpkzsum[:,5]),'--',color='green',label=r'$J_{{\bf k},n+1/2}-J_{{\bf k},n-1/2}$',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,1]),':',color='purple',label='Diss',basex=10,basey=10)
    n1p5=herm_grid**(-1.85)
    n1p5=10.0*np.abs(ebal_kpkzsum[20,6])/n1p5[20]*n1p5
    plt.loglog(herm_grid[10:100],n1p5[10:100],'-',color='black',basex=10,basey=10)
    plt.annotate(r'$\propto n^{-1.85}$',(25,12))
    plt.legend(loc='lower left',prop={'size':10}) 
    plt.xlabel('Hermite n')
    #plt.ylabel(r'$k_z L_n$')
    plt.title(r'$k_\perp \rho_i < 1$')
    plt.show()


    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy
    if test_eshells:
        ###########Temporary##############
        ###########Temporary##############
        ###########Temporary##############
        ebal_kzsum[2,:,:]=0.0
        ###########Temporary##############
        ###########Temporary##############
        ###########Temporary##############
        plt.plot(herm_grid,ebal_kpkzsum[:,3],'-',label='NL')
        plt.legend()
        plt.show()


        #plt.plot(herm_grid,Jptest_kpkzsum[:],'-',label='Test')
        #plt.plot(herm_grid,ebal_kpkzsum[:,5],'-',label='Pm2')
        #plt.legend()
        #plt.show()

        plt.plot(herm_grid,ebal_kpkzsum[:,0],'x-',label='RHS')
        plt.plot(herm_grid,ebal_kpkzsum[:,3],'-',label='NL')
        plt.plot(herm_grid,ebal_kpkzsum[:,1],'-',label='Diss')
        plt.plot(herm_grid,ebal_kpkzsum[:,4]+ebal_kpkzsum[:,5],'-',label='PM1+PM2')
        plt.plot(herm_grid,ebal_kpkzsum[:,4]\
                          +ebal_kpkzsum[:,5]\
                          +ebal_kpkzsum[:,3]\
                          +ebal_kpkzsum[:,1],'+-'\
                          ,label='Sum')
        plt.legend()
        plt.show()

        for i in range(num_shells):
            plt.title('shell '+str(i))
            plt.plot(herm_grid,ebal_kzsum[:,i,0],'x',label='RHS')
            plt.plot(herm_grid,ebal_kzsum[:,i,1],label='Diss')
            plt.plot(herm_grid,ebal_kzsum[:,i,3],label='NL')
            plt.plot(herm_grid,ebal_kzsum[:,i,4]\
                              +ebal_kzsum[:,i,5],label='PM1+PM2')
            plt.plot(herm_grid,ebal_kzsum[:,i,1]\
                              +ebal_kzsum[:,i,2]\
                              +ebal_kzsum[:,i,3]\
                              +ebal_kzsum[:,i,4]\
                              +ebal_kzsum[:,i,5]\
                              ,'+' ,label='Sum')
            plt.legend()
            plt.show()
            plt.title('shell '+str(i))
            plt.plot(herm_grid,ebal_kzsum[:,i,3]/ebal_kzsum[:,i,6],label='NL/E')
            plt.legend()
            plt.show()

    n1p5=herm_grid**(-1.5)

    if plot_hs_shells:
      for i in range(num_shells):
        plt.title('shell '+str(i))
        plt.loglog(herm_grid,np.abs(ebal_kzsum[:,i,1]),label='Diss',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(ebal_kzsum[:,i,3]),label='NL',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(ebal_kzsum[:,i,4]\
                          +ebal_kzsum[:,i,5]),'x',label='PM1+PM2',basex=10,basey=10)
        n1p5=5.0*np.abs(ebal_kzsum[20,i,4]+ebal_kzsum[20,i,5])/n1p5[20]*n1p5
        plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
        plt.legend()
        plt.show()

    kzn=np.empty((par['nv0'],num_shells))
    kzn_tot=np.empty((par['nv0']))
    kzn_tot=np.empty(par['nv0'])
    kznJm=np.empty((par['nv0'],num_shells))
    kznJp=np.empty((par['nv0'],num_shells))
    ebal_lsum=np.sum(ebal,axis=2)
    for n in range(par['nv0']):
        kzn_tot[n]=np.sum(np.abs(kzgrid)*ebal_lsum[n,:,6])/np.sum(ebal_lsum[n,:,6])
    plt.plot(herm_grid,kzn_tot,label='kzn tot.')
    plt.legend() 
    plt.xlabel(r'$k_z L$',size=18)
    plt.title(r'$<k_z>_n$')
    plt.show()



    for i in range(num_shells):
        for n in range(par['nv0']):
            kzn[n,i]=np.sum(np.abs(kzgrid)*ebal[n,:,i,6])/np.sum(ebal[n,:,i,6])
            kznJm[n,i]=np.sum(np.abs(kzgrid)*ebal[n,:,i,4])/np.sum(ebal[n,:,i,4])
            kznJp[n,i]=np.sum(np.abs(kzgrid)*ebal[n,:,i,5])/np.sum(ebal[n,:,i,5])
            #print kz0
            if show_plots:
                plt.plot(kzgrid,ebal[n,:,i,6],'x-',label='Free Energy')
                ax=plt.axis()
                plt.vlines(kzn[n,i],ax[2],ax[3])
                plt.title("Hermite n = "+str(n)+", Shell= "+str(i))
                plt.legend(loc='upper right')
                plt.xlabel(r'$k_z L$',size=18)
                plt.show()
        plt.plot(np.arange(par['nv0']),kzn[:,i],label=str(i))
    plt.legend() 
    plt.title(r'$<k_z>_n$')
    plt.show()

    

    color_list=['b','g','r','c','m','y','k','purple','maroon','teal']
    one=np.zeros(par['nv0'])+1.0
    for i in range(num_shells):
        kzeff=np.sum(kzn[20:60,i]/np.sqrt(herm_grid[20:60]))/40.0
        plt.plot(np.arange(par['nv0']),kzeff*one,label=str(i),color=color_list[i])
        plt.plot(np.arange(par['nv0']),kzn[0,i]*one,'--',label=str(i),color=color_list[i])
    plt.legend() 
    plt.title(r'$<k_z>_n/\sqrt{n}$')
    plt.show()

    for i in range(num_shells):
        kzeff=np.sum(kzn[20:60,i]/np.sqrt(herm_grid[20:60]))/40.0
        plt.plot(np.arange(par['nv0']),kzn[0,i]/kzeff*one,label=str(i),color=color_list[i])
    plt.legend() 
    plt.title(r'$<k_z>^{\phi}/<k_z>^{eff}$')
    plt.axis([0.0,100.0,0.0,6.0])
    plt.show()



    color_array=['b','g','r','c','m','black','orange','purple']
    shell_labels=[r'$k_\perp \rho_i = [0.0,0.2)$', r'$k_\perp \rho_i = [0.2,0.4)$', r'$k_\perp \rho_i = [0.4,0.6)$', r'$k_\perp \rho_i = [0.6,0.8)$', r'$k_\perp \rho_i = [0.8,1.0)$', r'$k_\perp \rho_i = [1.2,1.4)$', r'$k_\perp \rho_i \geq 1.4$']

    plt.figure(figsize=(4.5,3.5))
    fig=plt.gcf()
    fig.subplots_adjust(left=0.2,bottom=0.2)
    for i in range(1,5):
        plt.plot(np.arange(par['nv0']),kzn[:,i]/np.sqrt(herm_grid),color=color_array[i-1],label=shell_labels[i])
        if i==1:
          plt.plot(0,kzn[0,i],'x',markersize=8,markeredgewidth=2.0,color=color_array[i-1],label=r'$k_z^{(\phi)}$')
        else:
          plt.plot(0,kzn[0,i],'x',markersize=8,markeredgewidth=2.0,color=color_array[i-1])
    plt.legend(prop={'size':8}) 
    plt.xlabel('Hermite n')
    plt.ylabel(r'$k_z L_n$')
    plt.title(r'$<k_z>_n/\sqrt{n}$')
    plt.show()

    for i in range(num_shells):
        plt.plot(np.arange(par['nv0']),kzn[:,i]/np.sqrt(herm_grid)/kzn[0,i],label=str(i))
    plt.title('TEST!')
    plt.legend() 
    plt.title(r'$<k_z>_n/\sqrt{n}$')
    plt.show()


    for i in range(num_shells):
        plt.plot(np.arange(par['nv0']),kzn[:,i]/herm_grid,label=str(i))
    plt.legend() 
    plt.title(r'$<k_z>_n/n$')
    plt.show()

    sigmam=np.zeros((par['nv0'],par['nkz0'],num_shells),dtype='float')
    delm=np.zeros((par['nv0'],par['nkz0'],num_shells),dtype='float')
    for n in range(1,par['nv0']):
        delm[n,:,:]=np.sqrt(ebal[n-1,:,:,6])-np.sqrt(ebal[n,:,:,6])
        for k in range(1,par['nkz0']):
            sigmam[n,k,:]=ebal[n,k,:,4]/(-np.abs(kzgrid[k])*np.sqrt(float(n))\
                          *np.sqrt(ebal[n,k,:,6])*(np.sqrt(ebal[n,k,:,6])+delm[n,k,:]))
        sigmam[n,par['nkz0']/2,:]=0.0
            #print "sigmam[n,k,:]",sigmam[n,k,:]

    sigmam_nkp=np.sum(sigmam,axis=1)
    for i in range(num_shells):
        plt.plot(sigmam_nkp[:,i],label=str(i))
    plt.title('sigmam')
    plt.legend()
    plt.show()

    kzsigma_n=np.zeros((par['nv0'],num_shells),dtype='float')
    for i in range(num_shells):
        for n in range(par['nv0']):
            kzsigma_n[n,i]=np.sum(np.abs(kzgrid)*ebal[n,:,i,6]*sigmam[n,:,i])/np.sum(ebal[n,:,i,6])
            #for k in range(par['nkz0']):
            #    print k,sigmam[n,k,i]
            #print "numer1",np.sum(sigmam[n,:,i]) 
            #print "numer2",np.sum(np.abs(kzgrid)*sigmam[n,:,i]) 
            #print "numer",np.sum(np.abs(kzgrid)*ebal[n,:,i,6]*sigmam[n,:,i]) 
            #print "denom",np.sum(ebal[n,:,i,6])
            #plt.plot(herm_grid,sigmam[:,k,i],label=str(k*par['kzmin']))
            #plt.title('sigmam')
            #plt.legend()
            #plt.show()
    for i in range(num_shells):
        plt.title('|kzsigma_n|')
        #print "kzsigma_n",kzsigma_n[:,i]
        plt.plot(herm_grid,np.abs(kzsigma_n[:,i]),label=str(i))
    plt.legend()
    plt.show()

    for i in range(num_shells):
        plt.title('|kzsigma_n|/sqrt(n)')
        #print "kzsigma_n",kzsigma_n[:,i]
        plt.plot(herm_grid,np.abs(kzsigma_n[:,i])/np.sqrt(herm_grid),label=str(i))
    plt.legend()
    plt.show()

    for i in range(num_shells):
        plt.title('|kzsigma_n|/sqrt(n)')
        #print "kzsigma_n",kzsigma_n[:,i]
        plt.plot(herm_grid,np.abs(kzsigma_n[:,i])/np.sqrt(herm_grid)/kzn[0,i],label=str(i))
    plt.title('TEST!')
    plt.legend()
    plt.show()



    for i in range(num_shells):
        plt.title('|kzsigma_n|/(n)')
        #print "kzsigma_n",kzsigma_n[:,i]
        plt.plot(herm_grid,np.abs(kzsigma_n[:,i])/(herm_grid),label=str(i))
    plt.legend()
    plt.show()

    EqVI_kpkzsum=np.zeros(par['nv0'],dtype='float')
    EqV_kpkzsum=np.zeros(par['nv0'],dtype='float')
    for n in range(par['nv0']-1):
        EqVI_kpkzsum[n]=np.sum(kzn[n,:]*ebal_kzsum[n,:,6]-kzn[n+1]*ebal_kzsum[n+1,:,6])
        EqV_kpkzsum[n]=np.sum(kzsigma_n[n,:]*ebal_kzsum[n,:,6]-kzsigma_n[n+1]*ebal_kzsum[n+1,:,6])

    if end_shell != -1:
        plt.title('end shell = '+str(end_shell)+' of '+str(num_shells))
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,0]),':',label='RHS',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,3]),'-',label='NL',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,1]),'-',label='Diss',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(EqVI_kpkzsum),'-',label='EqVI',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(EqV_kpkzsum),'-',label='EqV',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,4]+ebal_kpkzsum[:,5]),'x-',label='PM1+PM2',basex=10,basey=10)
    n1p5=herm_grid**(-1.5)
    n1p5=2.0*np.abs(ebal_kpkzsum[20,4]+ebal_kpkzsum[20,5])/n1p5[20]*n1p5
    plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
    plt.loglog(herm_grid,np.abs(ebal_kpkzsum[:,4]\
                      +ebal_kpkzsum[:,5]\
                      +ebal_kpkzsum[:,3]\
                      +ebal_kpkzsum[:,1]),'-.'\
                      ,label='Sum')
    plt.legend()
    plt.show()

    kbounds=np.genfromtxt(par['diagdir'][1:-1]+'/shell_info.dat')
    ks_grid=np.empty(len(kbounds)+1)
    ks_grid[0]=kbounds[0]/2.0
    for i in range(1,len(ks_grid)-1):
        ks_grid[i]=kbounds[i-1]+(kbounds[i]-kbounds[i-1])/2.0
    dkperp=kbounds[-1]-kbounds[-2]
    ks_grid[-1]=kbounds[-1]+dkperp
    print "kbounds",kbounds
    print "ks_grid",ks_grid
    plt.plot(ks_grid,kzn[0,:],label=r'$<k_z^{\phi}>$')
    plt.xlabel(r'$k_{\perp} \rho_i$')
    plt.legend()
    plt.show()




    EqVI=np.zeros((par['nv0'],num_shells),dtype='float')
    EqV=np.zeros((par['nv0'],num_shells),dtype='float')
    for i in range(num_shells):
        for n in range(par['nv0']-1):
            EqVI[n,i]=np.sum(np.abs(kzgrid[:])*(ebal[n,:,i,6]*np.sqrt(n)-\
                             ebal[n+1,:,i,6]*np.sqrt(n+1)))
            EqV[n,i]=-np.sum(np.abs(kzgrid[:])*sigmam[n,:,i]*(ebal[n,:,i,6]*np.sqrt(n)-\
                             ebal[n+1,:,i,6]*np.sqrt(n+1)))
        #print "EqV(i)",EqV[:,i]
        plt.loglog(herm_grid,np.sum(ebal[:,:,i,4]+ebal[:,:,i,5],axis=1),label='PM',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(EqV[:,i]),label='EqV',basex=10,basey=10)
        plt.loglog(herm_grid,np.abs(EqVI[:,i]),label='EqVI',basex=10,basey=10)
        plt.legend()
        plt.show()

   
    #for i in range(num_shells):
    #    plt.plot(np.arange(par['nv0']),kznJm[:,i],label=str(i))
    #plt.legend() 
    #plt.title(r'$<k_z>_n(J_{n-1/2})$')
    #plt.show()
    #for i in range(num_shells):
    #    plt.plot(np.arange(par['nv0']),kznJm[:,i]/np.sqrt(herm_grid),label=str(i))
    #plt.legend() 
    #plt.title(r'$<k_z>_n/\sqrt{n}(J_{n-1/2})$')
    #plt.show()

    #for i in range(num_shells):
    #    plt.plot(np.arange(par['nv0']),kznJm[:,i]/herm_grid,label=str(i))
    #plt.legend() 
    #plt.title(r'$<k_z>_n/n(J_{n-1/2})$')
    #plt.show()




    if test_eshells:
        kzphi=np.empty(par['nv0'],dtype='float')
        for i in range(num_shells):
            kzphi[:]=kzn[0,i]
            plt.loglog(herm_grid,np.abs((ebal_kzsum[:,i,4])+(ebal_kzsum[:,i,5])),label=str(i)+' '+r'$J_{n-1/2}+J_{n+1/2}$')
            plt.loglog(herm_grid,np.abs(ebal_kzsum[:,i,3]),label='NL')
            plt.loglog(herm_grid,np.abs(ebal_kzsum[:,i,1]),label='Diss')
            plt.loglog(herm_grid,np.abs(ebal_kzsum[:,i,0]),label='RHS')
            plt.loglog(herm_grid,herm_grid**(-0.5),'--',label='n^-0.5')
            plt.loglog(herm_grid,herm_grid**(-1.0),'--',label='n^-1.0')
            plt.loglog(herm_grid,herm_grid**(-1.5),'--',label='n^-1.5')
            plt.legend() 
            plt.show()
            plt.plot(herm_grid,kzphi,label='kzphi')
            plt.plot(herm_grid,np.abs(ebal_kzsum[:,i,3]/ebal_kzsum[:,i,6]),label='N/E')
            plt.legend()
            plt.show()
        #    PM1=2.0*kzn[:,i]*np.sqrt(herm_grid)*ebal_kzsum[:,i,6]
        #    PM2=-np.roll(2.0*kzn[:,i]*np.sqrt(herm_grid)*ebal_kzsum[:,i,6],-1)
        #    plt.plot(np.arange(par['nv0']),PM1+PM2,label=str(i)+' '+r'$PM1+PM2$')
        #    plt.plot(herm_grid,ebal_kzsum[:,i,3],label='NL')
        #    plt.plot(herm_grid,ebal_kzsum[:,i,1],label='Diss')
        #    plt.plot(np.arange(par['nv0']),PM1+PM2+ebal_kzsum[:,i,3]+ebal_kzsum[:,i,1],'x',label='SUM')
        #    plt.plot(np.arange(par['nv0']),ebal_kzsum[:,i,0],label='RHS')
        #    plt.legend() 
        #    plt.show()

        #for i in range(num_shells):
        #    plt.plot(np.arange(par['nv0']),2.0*kzn[:,i]*np.sqrt(herm_grid)*ebal_kzsum[:,i,6],label=str(i)+' '+r'$<k_z>_n \varepsilon_n$')
        #    plt.plot(np.arange(par['nv0']),np.abs(ebal_kzsum[:,i,4]),label=str(i)+' '+r'$J_{n-1/2}$')
        #    plt.plot(np.arange(par['nv0']),np.abs(ebal_kzsum[:,i,4])/(2.0*kzn[:,i]*np.sqrt(herm_grid)*ebal_kzsum[:,i,6]),label=str(i)+' '+r'$2/1')
        #    plt.legend() 
        #    plt.show()
        #    plt.plot(np.arange(par['nv0']),2.0*kzn[:,i]*np.sqrt(herm_grid+1)*np.roll(ebal_kzsum[:,i,6],-1),label=str(i)+' '+r'$<k_z>_n \sqrt{n} \varepsilon_n$')
        #    plt.plot(np.arange(par['nv0']),np.abs(ebal_kzsum[:,i,5]),label=str(i)+' '+r'$J_{n+1/2}$')
        #    plt.plot(np.arange(par['nv0']),np.abs(ebal_kzsum[:,i,5])/(2.0*kzn[:,i]*np.sqrt(herm_grid+1)*np.roll(ebal_kzsum[:,i,6],-1)),label=str(i)+' '+r'$2/1$')
        #    plt.legend() 
        #    plt.show()


    #for i in range(num_shells):
    #    plt.plot(np.arange(par['nv0']),kzn[:,i]/np.sqrt(herm_grid),label=str(i))
    #plt.legend() 
    #plt.title(r'$<k_z>_n/\sqrt{n}$')
    #plt.show()

    for i in range(num_shells):
        temp=np.sum(ebal[:,:,i,5],axis=1)
        plt.plot(np.arange(par['nv0']),temp,label=str(i))
    plt.legend() 
    plt.title(r'$J_{n+1/2}$')
    plt.show()

    Etot=np.empty(par['nv0'])
    Jptot=np.empty(par['nv0'])
    Disstot=np.empty(par['nv0'])
    #kzntot=np.empty(par['nv0'])
    for i in range(par['nv0']):
        Etot[i]=np.sum(ebal[i,:,:,6])
        Jptot[i]=np.sum(ebal[i,:,:,5])
        Disstot[i]=np.sum(ebal[i,:,:,1])
    #    kzntot[i]=np.sum(kzn[i,:])


    n0p5=500.0*herm_grid**(-0.5)
    n1p0=500.0*herm_grid**(-1.0)
    n1p5=500.0*herm_grid**(-1.5)

    print np.shape(Etot)
    print np.shape(Jptot)
    #print np.shape(kzntot)
    #plt.plot(herm_grid,Etot*kzntot,label='Etot*kzntot')
    #plt.plot(herm_grid,-Jptot,label='Jptot')
    #plt.legend()
    #plt.show()
    plt.loglog(herm_grid,Etot,label='Etot',basex=10,basey=10)
    plt.loglog(herm_grid,n1p5,label='-1.5',basex=10,basey=10)
    plt.legend()
    plt.show()
    plt.loglog(herm_grid,-Jptot,label='-Jptot',basex=10,basey=10)
    plt.loglog(herm_grid,n0p5,'--',label='-0.5',basex=10,basey=10)
    #plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
    #plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
    plt.legend()
    plt.show()
    plt.loglog(herm_grid,-Disstot,label='-Disstot',basex=10,basey=10)
    plt.loglog(herm_grid,n0p5,'--',label='-0.5',basex=10,basey=10)
    plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
    plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
    plt.legend()
    plt.show()
    plt.loglog(herm_grid,Jptot-np.roll(Jptot,1),'x-',label='ddn-Jptot',basex=10,basey=10)
    plt.loglog(herm_grid,-Disstot,label='-Disstot',basex=10,basey=10)
    plt.loglog(herm_grid,n0p5,'--',label='-0.5',basex=10,basey=10)
    plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
    plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
    plt.legend()
    plt.show()
    plt.loglog(herm_grid,Jptot-np.roll(Jptot,1),'x-',label='ddn-Jptot',basex=10,basey=10)
    plt.loglog(herm_grid,herm_grid*Etot,'x-',label='n*E',basex=10,basey=10)
    #plt.loglog(herm_grid,herm_grid*Etot-np.roll(herm_grid*Etot,1),'x-',label='ddn n*E',basex=10,basey=10)
    plt.loglog(herm_grid,-Disstot,label='-Disstot',basex=10,basey=10)
    plt.loglog(herm_grid,n0p5,'--',label='-0.5',basex=10,basey=10)
    plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
    plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
    plt.legend()
    plt.show()
    plt.loglog(herm_grid,-Disstot/Etot,label='-Disstot/Etot',basex=10,basey=10)
    plt.loglog(herm_grid,n0p5,'--',label='-0.5',basex=10,basey=10)
    plt.loglog(herm_grid,n1p0,'--',label='-1.0',basex=10,basey=10)
    plt.loglog(herm_grid,n1p5,'--',label='-1.5',basex=10,basey=10)
    plt.legend()
    plt.show()

    for i in range(num_shells):
        temp=np.sum(ebal[:,:,i,6],axis=1)
        plt.loglog(np.arange(par['nv0']),temp,label=str(i),basex=10,basey=10)
        plt.loglog(np.arange(par['nv0']),n0p5,'--',label='-0.5',basex=10,basey=10)
        plt.loglog(np.arange(par['nv0']),n1p0,':',label='-1.0',basex=10,basey=10)
        plt.loglog(np.arange(par['nv0']),n1p5,'-.',label='-1.5',basex=10,basey=10)
        plt.legend() 
        plt.title(r'$\varepsilon_n$')
        plt.show()

    #for i in range(num_shells):
    #    temp=np.sum(ebal[:,:,i,5],axis=1)
    #    temp+=np.sum(ebal[:,:,i,4],axis=1)
    #    temp2=np.sum(ebal[:,:,i,3],axis=1)
    #    temp3=np.sum(ebal[:,:,i,1],axis=1)
    #    #temp4=np.sum(ebal[:,:,i,2],axis=1)
    #    plt.plot(np.arange(par['nv0']),temp,label='J')
    #    plt.plot(np.arange(par['nv0']),temp2,label='N')
    #    plt.plot(np.arange(par['nv0']),temp3,label='C')
    #    #plt.plot(np.arange(par['nv0']),temp3+temp+temp2+temp4,'x',label='Sum')
    #    plt.title('shell '+str(i))
    #    plt.legend() 
    #    plt.show()

def eshells_cb_plot(start_time=-1,end_time=-1,show_plots=False,test_eshells=False,end_shell=-1,plot_hs_shells=True):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
    if end_shell==-1:
        end_shell=num_shells
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)

    ebal=np.zeros((par['nv0'],par['nkz0'],num_shells,7))
    #order in ebal
    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy

    es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    print "np.shape(es_in)"
    print np.shape(es_in)
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,6]+=es_in
        #RHS
        f.seek(8+i*(8+mem_tot)+ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,0]+=es_in
        #ebal[:,:,0]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,1]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,1]+=es_in
        #ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,3]+=es_in
        #ebal[:,:,3]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,4]+=es_in
        #ebal[:,:,4]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
        ebal[:,:,:,5]+=es_in
        #ebal[:,:,5]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #if herm_n==2:
        f.seek(8+i*(8+mem_tot)+7*ntot*8)
        #Drive
        fs_in=np.zeros((par['nkz0'],num_shells))
        fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
        fs_in=np.reshape(fs_in,(par['nkz0'],num_shells))
        for j in range(num_shells):
            ebal[2,:,:,2]+=fs_in[:,:]
        f.close()

    ebal=ebal/float(ntime)
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    ebal_kzsum=np.sum(ebal,axis=1)    
    ebal_kpkzsum=np.sum(np.sum(ebal[:,:,0:end_shell],axis=1),axis=1)    

    #Jptest=np.empty((par['nv0'],par['nkz0'],num_shells))
    #for k in range(par['nkz0']):
    #    for n in range(par['nv0']):
    #        Jptest[n,k,:]=2.0*np.sqrt(n)*kzgrid[k]*ebal[n,k,:,6]
    #ebal_kzsum=np.sum(ebal,axis=1)    
    #Jptest_kzsum=np.sum(Jptest,axis=1)    
    #ebal_kpkzsum=np.sum(np.sum(ebal,axis=1),axis=1)    
    #Jptest_kpkzsum=np.sum(np.sum(Jptest,axis=1),axis=1)    

    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy
    if test_eshells:
        ###########Temporary##############
        ###########Temporary##############
        ###########Temporary##############
        ebal_kzsum[2,:,:]=0.0
        ###########Temporary##############
        ###########Temporary##############
        ###########Temporary##############


        #plt.plot(herm_grid,Jptest_kpkzsum[:],'-',label='Test')
        #plt.plot(herm_grid,ebal_kpkzsum[:,5],'-',label='Pm2')
        #plt.legend()
        #plt.show()


    n1p5=herm_grid**(-1.5)

    kzn=np.empty((par['nv0'],num_shells))
    kzn_tot=np.empty((par['nv0']))
    kzn_tot=np.empty(par['nv0'])
    kznJm=np.empty((par['nv0'],num_shells))
    kznJp=np.empty((par['nv0'],num_shells))
    ebal_lsum=np.sum(ebal,axis=2)
    for n in range(par['nv0']):
        kzn_tot[n]=np.sum(np.abs(kzgrid)*ebal_lsum[n,:,6])/np.sum(ebal_lsum[n,:,6])
#    plt.plot(herm_grid,kzn_tot,label='kzn tot.')
#    plt.legend() 
#    plt.xlabel(r'$k_z L$',size=18)
#    plt.title(r'$<k_z>_n$')
#    plt.show()



    for i in range(num_shells):
        for n in range(par['nv0']):
            kzn[n,i]=np.sum(np.abs(kzgrid)*ebal[n,:,i,6])/np.sum(ebal[n,:,i,6])
            kznJm[n,i]=np.sum(np.abs(kzgrid)*ebal[n,:,i,4])/np.sum(ebal[n,:,i,4])
            kznJp[n,i]=np.sum(np.abs(kzgrid)*ebal[n,:,i,5])/np.sum(ebal[n,:,i,5])
            #print kz0
#            if show_plots:
#                plt.plot(kzgrid,ebal[n,:,i,6],'x-',label='Free Energy')
#                ax=plt.axis()
#                plt.vlines(kzn[n,i],ax[2],ax[3])
#                plt.title("Hermite n = "+str(n)+", Shell= "+str(i))
#                plt.legend(loc='upper right')
#                plt.xlabel(r'$k_z L$',size=18)
#                plt.show()
#        plt.plot(np.arange(par['nv0']),kzn[:,i],label=str(i))
#    plt.legend() 
#    plt.title(r'$<k_z>_n$')
#    plt.show()

    

    color_list=['b','g','r','c','m','y','k','purple','maroon','teal']
    one=np.zeros(par['nv0'])+1.0
#    for i in range(num_shells):
#        kzeff=np.sum(kzn[20:60,i]/np.sqrt(herm_grid[20:60]))/40.0
#        plt.plot(np.arange(par['nv0']),kzeff*one,label=str(i),color=color_list[i])
#        plt.plot(np.arange(par['nv0']),kzn[0,i]*one,'--',label=str(i),color=color_list[i])
#    plt.legend() 
#    plt.title(r'$<k_z>_n/\sqrt{n}$')
#    plt.show()

    for i in range(num_shells):
        kzeff=np.sum(kzn[20:60,i]/np.sqrt(herm_grid[20:60]))/40.0
#        plt.plot(np.arange(par['nv0']),kzn[0,i]/kzeff*one,label=str(i),color=color_list[i])
#    plt.legend() 
#    plt.title(r'$<k_z>^{\phi}/<k_z>^{eff}$')
#    plt.axis([0.0,100.0,0.0,6.0])
#    plt.show()



    shell_labels=[r'$k_\perp \rho_i = [0.0,0.2)$', r'$k_\perp \rho_i = [0.2,0.4)$', r'$k_\perp \rho_i = [0.4,0.6)$', r'$k_\perp \rho_i = [0.6,0.8)$', r'$k_\perp \rho_i = [0.8,1.0)$', r'$k_\perp \rho_i = [1.2,1.4)$', r'$k_\perp \rho_i \geq 1.4$']

    plt.figure(figsize=(4.5,3.5))
    fig=plt.gcf()
    fig.subplots_adjust(left=0.15,bottom=0.15)
    for i in range(1,5):
        plt.plot(np.arange(par['nv0']),kzn[:,i]/np.sqrt(herm_grid),color=color_list[i-1],label=shell_labels[i])
        plt.plot(0,kzn[0,i],'x',markersize=8,markeredgewidth=2.0,color=color_list[i-1])
    plt.plot(0,-1,'x',markersize=8,markeredgewidth=2.0,color='black',label=r'$k_z^{(\phi)}$')
    plt.legend(prop={'size':10}) 
    plt.xlabel('Hermite n')
    plt.axis([-10,120,0.0,0.7])
    #plt.ylabel(r'$k_z L_n$')
    plt.ylabel(r'$<k_z>_n/\sqrt{n}$')
    plt.show()



def eshells_cb_test_old(start_time=-1,end_time=-1,herm_n=0):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 
    print "num_shells",num_shells 

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)


    ebal=np.zeros((num_shells,par['nkz0'],7))
    #order in ebal
    #0:RHS
    #1:coll and hyps
    #2:Drive
    #3:NL
    #4:PH1
    #5:PH2
    #6:Free energy

    es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    print "np.shape(es_in)"
    print np.shape(es_in)
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[:,:,6]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #RHS
        f.seek(8+i*(8+mem_tot)+ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[:,:,0]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[:,:,1]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[:,:,3]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[:,:,4]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        #PH2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        ebal[:,:,5]+=distribute_ebal_shells_n(es_in,herm_n,keep_kz=True)
        if herm_n==2:
            f.seek(8+i*(8+mem_tot)+7*ntot*8)
            #Drive
            fs_in=np.zeros((par['nkz0'],num_shells))
            fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
            fs_in=np.reshape(fs_in,(par['nkz0'],num_shells))
            #print "np.shape(fs_in)",np.shape(fs_in)
            for j in range(num_shells):
                ebal[:,:,2]+=fs_in[:,:]
        else:
            ebal[:,:,2]=0.0
        f.close()

    ebal=ebal/float(ntime)
    ebal0=np.empty((num_shells,par['nkz0']/2-1,7))
    ebal0[:,0,:]=ebal[:,0,:]
    for i in range(par['nkz0']/2-2):
        ebal0[:,i+1,:]=ebal[:,i+1,:]+ebal[:,par['nkz0']-(i+1),:]
    kzgrid=np.arange(par['nkz0']/2-1)*par['kzmin']

    for i in range(num_shells):
        kz0=np.sum(kzgrid*ebal0[i,:,6])/np.sum(ebal0[i,:,6])
        print kz0
        plt.plot(kzgrid,ebal0[i,:,6],'x-',label='Free Energy')
        ax=plt.axis()
        plt.vlines(kz0,ax[2],ax[3])
        plt.title("Hermite n = "+str(herm_n)+", Shell= "+str(i))
        #plt.title("Hermite: "+str(herm_n))
        #if herm_n==0 or herm_n==1:
        #    plt.annotate("Note: discrepancies for n=0 and n=1 cancel.",[0,0])
        plt.legend(loc='upper right')
        plt.xlabel(r'$k_z L$',size=18)
        plt.show()

def energetics_kp_kz_n(start_time=-1,end_time=-1):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)

    CD=np.zeros((num_shells,par['nkz0'],par['nv0']))
    CD2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))
    NL=np.zeros((num_shells,par['nkz0'],par['nv0']))
    NL2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))
    PM=np.zeros((num_shells,par['nkz0'],par['nv0']))
    PM2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))
    DR=np.zeros((num_shells,par['nkz0']))
    DR2=np.zeros((num_shells,par['nkz0']/2))
    FE=np.zeros((num_shells,par['nkz0'],par['nv0']))
    FE2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))
    HD=np.zeros((num_shells,par['nkz0'],par['nv0']))
    HD2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))

    #From nlt_shells files:
    NL0=np.zeros((num_shells,par['nkz0']/2,11))

    #es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        FE=FE+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        CD=CD+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        HD=HD+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        NL=NL+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #PM1
        #f.seek(8+i*(8+mem_tot)+5*ntot*8)
        #es_in=np.fromfile(f,dtype='float64',count=ntot)
        #PM=PM+np.abs(np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F'))
        #PM2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        PM=PM+np.abs(np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F'))
        #Drive
        f.seek(8+i*(8+mem_tot)+7*ntot*8)
        fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
        DR=DR+np.reshape(fs_in,(num_shells,par['nkz0']),order='F')

    FE=FE/float(ntime)
    CD=CD/float(ntime)
    HD=HD/float(ntime)
    NL=NL/float(ntime)
    PM=PM/float(ntime)
    DR=DR/float(ntime)
    for i in range(par['nkz0']/2):
        FE2[:,i,:]=FE[:,i,:]+FE[:,-i,:]
        CD2[:,i,:]=CD[:,i,:]+CD[:,-i,:]
        HD2[:,i,:]=HD[:,i,:]+HD[:,-i,:]
        NL2[:,i,:]=NL[:,i,:]+NL[:,-i,:]
        PM2[:,i,:]=PM[:,i,:]+PM[:,-i,:]
        DR2[:,i]=DR[:,i]+DR[:,-i]

    shell_grid=np.arange(num_shells)+0.5
    kzgrid=np.arange(par['nkz0']/2)/float(par['nkz0']/2-1)*par['kzmax0']
    nlt_ns=np.empty(11)
    nlt_ns[0]=0
    nlt_ns[1]=1
    nlt_ns[2]=2
    nlt_ns[3]=3
    nlt_ns[4]=par['nv0']/8
    nlt_ns[5]=2*par['nv0']/8
    nlt_ns[6]=3*par['nv0']/8
    nlt_ns[7]=4*par['nv0']/8
    nlt_ns[8]=5*par['nv0']/8
    nlt_ns[9]=6*par['nv0']/8
    nlt_ns[10]=7*par['nv0']/8
    nlt_files=list()
    nlt_files.append('nlt_shells_n0.dat')
    nlt_files.append('nlt_shells_n1.dat')
    nlt_files.append('nlt_shells_n2.dat')
    nlt_files.append('nlt_shells_n3.dat')
    nlt_files.append('nlt_shells_n1o8.dat')
    nlt_files.append('nlt_shells_n2o8.dat')
    nlt_files.append('nlt_shells_n3o8.dat')
    nlt_files.append('nlt_shells_n4o8.dat')
    nlt_files.append('nlt_shells_n5o8.dat')
    nlt_files.append('nlt_shells_n6o8.dat')
    nlt_files.append('nlt_shells_n7o8.dat')
    levels=[0.0,10.0]
    
    for i in range(11):
        n=nlt_ns[i]
        file=nlt_files[i]
        CDtemp=(CD2[:,:,n]/FE2[:,:,n])
        #CDtemp[:,par['nkz0']/2]=0.0
        #CDtemp=np.roll(CDtemp,par['nkz0']/2-1,axis=1)[:,0:-1]
        egmax=np.max(abs(CDtemp))
        for j in range(num_shells):
            for k in range(par['nkz0']/2):
                nlt0=nlt_test_l_kz(start_time=start_time,end_time=end_time,file_name=file,kz_in=k,shell_in=j,show_plots=False)
                #print np.shape(nlt0)
                #print np.shape(NL[:,:,i])
                #print NL[j,k,n],NL[j,-k,n]
                #print "nlt test: ",np.sum(nlt0),NL[j,k,n]+NL[j,-k,n]
                print "nlt test: ",np.sum(nlt0),NL2[j,k,n]
                NL0[j,k,i]=np.abs(np.sum(nlt0[nlt0 < 0.0]))  #Energy flowing out
                #NL0[j,k,i]=NL0[j,k,i]+np.abs(np.sum(nlt0[nlt0 > 0.0]))
                #NL0[j,k,i]=NL0[j,k,i]/2.0
                #########
                #print NL0[j,k,i]
                #nlsum=0.0
                #for p in range(num_shells):
                #    for q in range(par['nkz0']/2):
                #        if nlt0[p,q] > 0:
                #            nlsum+=nlt0[p,q]
                #print nlsum
                ##########

        NLtemp=(NL0[:,:,i]/FE2[:,:,n])
        #NLtemp[:,par['nkz0']/2]=0.0
        #NLtemp=np.roll(NLtemp,par['nkz0']/2-1,axis=1)[:,0:-1]
        egmax=max(egmax,np.max(abs(NLtemp)))

        PMtemp=(PM2[:,:,n]/FE2[:,:,n])
        #PMtemp[:,par['nkz0']/2]=0.0
        #PMtemp=np.roll(PMtemp,par['nkz0']/2-1,axis=1)[:,0:-1]
        egmax=max(egmax,np.max(abs(PMtemp)))
        clevels=np.arange(101)/100.0*2.0*egmax-egmax
        NL_PM=(NLtemp-PMtemp)
        CD_PM=(np.abs(CDtemp)-np.abs(PMtemp))

        if n==2:
            DRtemp=(DR2[:,:]/FE2[:,:,n])
            #DRtemp[:,par['nkz0']/2]=0.0
            #DRtemp=np.roll(DRtemp,par['nkz0']/2-1,axis=1)[:,0:-1]
            CD_DR=(CDtemp-DRtemp)
            NL_DR=(NLtemp-DRtemp)
            PM_DR=(PMtemp-DRtemp)
            egmax=max(egmax,np.max(abs(DRtemp)))
            clevels=np.arange(101)/100.0*2.0*egmax-egmax
            plt.title("DR")
            #plt.pcolor(DRtemp,cmap=my_color_map(),vmin=-egmax,vmax=egmax)
            plt.contourf(kzgrid,shell_grid,DRtemp,50)
            plt.xlabel('$k_z R$',size=18)
            plt.ylabel('$l$',size=18)
            plt.colorbar()
            plt.contour(kzgrid,shell_grid,NL_DR,levels,linewidths=3,cmap=white_color_map())
            plt.contour(kzgrid,shell_grid,PM_DR,levels,linewidths=3,cmap=green_color_map())
            plt.show()
            plt.title("DR (n="+str(n)+")")


        #plt.title("CD (n="+str(n)+")")
        ##print "min FE",np.min(FE[:,:,n])
        ##print "argmin FE",np.argmin(FE[:,:,n])
        ##plt.contour(CDtemp,clevels,cmap=my_color_map())
        ##plt.pcolor(CDtemp,cmap=my_color_map(),vmin=-egmax,vmax=egmax)
        #plt.contour(CDtemp)
        #plt.xlabel('$k_z R$',size=18)
        #plt.ylabel('$l$',size=18)
        #plt.colorbar()
        #plt.show()

        plt.title("PM (n="+str(n)+")")
        #temp=PM[:,:,n]/FE[:,:,n]
        #temp[:,par['nkz0']/2]=0.0
        #plt.pcolor(PMtemp,cmap=my_color_map(),vmin=-egmax,vmax=egmax)
        #print clevels
        #print PMtemp
        print np.shape(kzgrid)
        print np.shape(shell_grid)
        print np.shape(PMtemp)
        #plt.contourf(shell_grid,kzgrid,PMtemp,50)
        plt.contourf(kzgrid,shell_grid,PMtemp,50)
        plt.xlabel('$k_z R$',size=18)
        plt.ylabel('$l$',size=18)
        plt.colorbar()
        plt.contour(kzgrid,shell_grid,NL_PM,levels,linewidths=3,cmap=black_color_map())
        plt.contour(kzgrid,shell_grid,CD_PM,levels,linewidths=3,cmap=white_color_map())
        plt.show()


        plt.title("NL (n="+str(n)+")")
        #plt.pcolor(NLtemp,cmap=my_color_map(),vmin=-egmax,vmax=egmax)
        plt.contourf(kzgrid,shell_grid,NLtemp,50)
        plt.xlabel('$k_z R$',size=18)
        plt.ylabel('$l$',size=18)
        plt.colorbar()
        plt.contour(kzgrid,shell_grid,NL_PM,levels,linewidths=3,cmap=black_color_map())
        plt.show()

def energetics_3d_flow(start_time=-1,end_time=-1,which_norm='unit'):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)

    CD=np.zeros((num_shells,par['nkz0'],par['nv0']))
    CD2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))
    NL=np.zeros((num_shells,par['nkz0'],par['nv0']))
    NL2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))
    PM=np.zeros((num_shells,par['nkz0'],par['nv0']))
    PM2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))
    DR=np.zeros((num_shells,par['nkz0']))
    DR2=np.zeros((num_shells,par['nkz0']/2))
    FE=np.zeros((num_shells,par['nkz0'],par['nv0']))
    FE2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))
    HD=np.zeros((num_shells,par['nkz0'],par['nv0']))
    HD2=np.zeros((num_shells,par['nkz0']/2,par['nv0']))

    #From nlt_shells files:
    NL0=np.zeros((num_shells,par['nkz0']/2,11))

    #es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        FE=FE+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        CD=CD+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        HD=HD+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        NL=NL+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #PM1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        PM=PM+(np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F'))
        #PM2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        PM=PM+(np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F'))
        #Drive
        f.seek(8+i*(8+mem_tot)+7*ntot*8)
        fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
        DR=DR+np.reshape(fs_in,(num_shells,par['nkz0']),order='F')
    
    FE=FE/float(ntime)
    CD=CD/float(ntime)
    HD=HD/float(ntime)
    NL=NL/float(ntime)
    PM=PM/float(ntime)
    DR=DR/float(ntime)
    for i in range(par['nkz0']/2):
        FE2[:,i,:]=FE[:,i,:]+FE[:,-i,:]
        CD2[:,i,:]=CD[:,i,:]+CD[:,-i,:]
        HD2[:,i,:]=HD[:,i,:]+HD[:,-i,:]
        NL2[:,i,:]=NL[:,i,:]+NL[:,-i,:]
        PM2[:,i,:]=PM[:,i,:]+PM[:,-i,:]
        DR2[:,i]=DR[:,i]+DR[:,-i]

    shell_grid=np.arange(num_shells)+0.5
    kzgrid=np.arange(par['nkz0']/2)/float(par['nkz0']/2-1)*par['kzmax0']
    nlt_ns=np.empty(11)
    nlt_ns[0]=0
    nlt_ns[1]=1
    nlt_ns[2]=2
    nlt_ns[3]=3
    nlt_ns[4]=par['nv0']/8
    nlt_ns[5]=2*par['nv0']/8
    nlt_ns[6]=3*par['nv0']/8
    nlt_ns[7]=4*par['nv0']/8
    nlt_ns[8]=5*par['nv0']/8
    nlt_ns[9]=6*par['nv0']/8
    nlt_ns[10]=7*par['nv0']/8
    nlt_files=list()
    nlt_files.append('nlt_shells_n0.dat')
    nlt_files.append('nlt_shells_n1.dat')
    nlt_files.append('nlt_shells_n2.dat')
    nlt_files.append('nlt_shells_n3.dat')
    nlt_files.append('nlt_shells_n1o8.dat')
    nlt_files.append('nlt_shells_n2o8.dat')
    nlt_files.append('nlt_shells_n3o8.dat')
    nlt_files.append('nlt_shells_n4o8.dat')
    nlt_files.append('nlt_shells_n5o8.dat')
    nlt_files.append('nlt_shells_n6o8.dat')
    nlt_files.append('nlt_shells_n7o8.dat')
    levels=[0.0,10.0]

    kp_loc=np.empty(0)
    kz_loc=np.empty(0)
    n_loc=np.empty(0)
    kp_comp=np.empty(0)
    kz_comp=np.empty(0)
    n_comp=np.empty(0)
    sources_sinks=np.empty(0)
    
    for i in range(11):
        print i," of ",10
        n=nlt_ns[i]
        file=nlt_files[i]
        for j in range(num_shells):
            for k in range(par['nkz0']/2):
                nl_out,nkp0,nkz0=nlt_out_l_kz(start_time=start_time,end_time=end_time,file_name=file,kz_in=k,shell_in=j,show_plots=False)
                kp_loc=np.append(kp_loc,0.2*j+0.1)
                kz_loc=np.append(kz_loc,par['kzmin']*k)
                n_loc=np.append(n_loc,0)
                #n_loc=np.append(n_loc,n)
                if which_norm=='energy':
                    so_si=CD2[j,k,n]/FE[j,k,n]
                    kp_comp=np.append(kp_comp,nkp0/FE[j,k,n])
                    kz_comp=np.append(kz_comp,nkz0/FE[j,k,n])
                    n_comp=np.append(n_comp,PM2[j,k,n]/FE[j,k,n])
                    if n==2:
                        so_si+=DR2[j,k]/FE[j,k,n]
                elif which_norm=='none':
                    so_si=CD2[j,k,n]
                    kp_comp=np.append(kp_comp,nkp0)
                    kz_comp=np.append(kz_comp,nkz0)
                    n_comp=np.append(n_comp,PM2[j,k,n])
                    if n==2:
                        so_si+=DR2[j,k]
                elif which_norm=='unit':
                    so_si=CD2[j,k,n]
                    magnitude=(nkp0**2+nkz0**2)**0.5
                    kp_comp=np.append(kp_comp,nkp0/magnitude)
                    kz_comp=np.append(kz_comp,nkz0/magnitude)
                    n_comp=np.append(n_comp,PM2[j,k,n]/FE[j,k,n])
                    if n==2:
                        so_si+=DR2[j,k]
                else:
                    stop

                sources_sinks=np.append(sources_sinks,so_si)

    vtk_writer=VTK_XML_Serial_Unstructured()
    ntot=num_shells*par['nkz0']/2
    zeros=np.zeros(ntot)
    for i in range(11):
        n=nlt_ns[i]
        file=nlt_files[i]
        vtk_writer.snapshot("fenergy3d_"+str(n)+".vtu",kp_loc[i*ntot:(i+1)*ntot],kz_loc[i*ntot:(i+1)*ntot],n_loc[i*ntot:(i+1)*ntot],\
                            x_jump=kp_comp[i*ntot:(i+1)*ntot],y_jump=kz_comp[i*ntot:(i+1)*ntot],z_jump=zeros,radii=sources_sinks[i*ntot:(i+1)*ntot],colors=n_comp[i*ntot:(i+1)*ntot])
        #vtk_writer.snapshot("fenergy3d_"+str(n)+".vtu",kp_loc[i*ntot:(i+1)*ntot],kz_loc[i*ntot:(i+1)*ntot],n_loc[i*ntot:(i+1)*ntot],\
        #                    x_jump=kp_comp[i*ntot:(i+1)*ntot],y_jump=kz_comp[i*ntot:(i+1)*ntot],z_jump=n_comp[i*ntot:(i+1)*ntot],colors=sources_sinks[i*ntot:(i+1)*ntot])
    vtk_writer.writePVD("fenergy3d_test.pvd")

def energetics_kp_kz_n_old(start_time=-1,end_time=-1):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print num_shells
    f.close()
 
    f = open('eshells.dat','rb')
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
      else:
          continue_read=0
 
    f.close()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 


    #Order from eshells file
    #0:Energy
    #1:RHS
    #2:coll
    #3:hypercoll & hyps
    #4:NL
    #5:PH 1
    #6:PH 2
    #7:Drive (no v dependence)

    CD=np.zeros((num_shells,par['nkz0'],par['nv0']))
    NL=np.zeros((num_shells,par['nkz0'],par['nv0']))
    PM=np.zeros((num_shells,par['nkz0'],par['nv0']))
    DR=np.zeros((num_shells,par['nkz0']))
    FE=np.zeros((num_shells,par['nkz0'],par['nv0']))
    HD=np.zeros((num_shells,par['nkz0'],par['nv0']))

    #es_in=np.zeros((par['nv0'],par['nkz0'],num_shells))
    
    ntot=num_shells*par['nkz0']*par['nv0']
    mem_tot=7*ntot*8+num_shells*par['nkz0']*8

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open('eshells.dat','rb')
        #Energy
        f.seek(8+i*(8+mem_tot))
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        FE=FE+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #coll
        f.seek(8+i*(8+mem_tot)+2*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        CD=CD+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #hyps
        f.seek(8+i*(8+mem_tot)+3*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        HD=HD+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #NL
        f.seek(8+i*(8+mem_tot)+4*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        NL=NL+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #PH1
        f.seek(8+i*(8+mem_tot)+5*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        PM=PM+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #PH2
        f.seek(8+i*(8+mem_tot)+6*ntot*8)
        es_in=np.fromfile(f,dtype='float64',count=ntot)
        PM=PM+np.reshape(es_in,(num_shells,par['nkz0'],par['nv0']),order='F')
        #Drive
        f.seek(8+i*(8+mem_tot)+7*ntot*8)
        fs_in=np.fromfile(f,dtype='float64',count=num_shells*par['nkz0'])
        DR=DR+np.reshape(fs_in,(num_shells,par['nkz0']),order='F')

    FE=FE/float(ntime)
    CD=CD/float(ntime)
    HD=HD/float(ntime)
    NL=NL/float(ntime)
    PM=PM/float(ntime)
    DR=DR/float(ntime)

    shell_grid=np.arange(num_shells+1)+0.5
    
    for n in range(par['nv0']):
        CDtemp=(CD[:,:,n]/FE[:,:,n])
        CDtemp[:,par['nkz0']/2]=0.0
        CDtemp=np.roll(CDtemp,par['nkz0']/2-1,axis=1)[:,0:-1]
        egmax=np.max(abs(CDtemp))
        NLtemp=(NL[:,:,n]/FE[:,:,n])
        NLtemp[:,par['nkz0']/2]=0.0
        NLtemp=np.roll(NLtemp,par['nkz0']/2-1,axis=1)[:,0:-1]
        egmax=max(egmax,np.max(abs(NLtemp)))
        PMtemp=(PM[:,:,n]/FE[:,:,n])
        PMtemp[:,par['nkz0']/2]=0.0
        PMtemp=np.roll(PMtemp,par['nkz0']/2-1,axis=1)[:,0:-1]
        egmax=max(egmax,np.max(abs(PMtemp)))
        CD_NL=(CDtemp-NLtemp)
        CD_PM=(CDtemp-PMtemp)
        NL_PM=(NLtemp-PMtemp)
        levels=[0.0,10.0]
        clevels=np.arange(101)/100.0*2.0*egmax-egmax


        if n==2:
            DRtemp=(DR[:,:]/FE[:,:,n])
            DRtemp[:,par['nkz0']/2]=0.0
            DRtemp=np.roll(DRtemp,par['nkz0']/2-1,axis=1)[:,0:-1]
            CD_DR=(CDtemp-DRtemp)
            NL_DR=(NLtemp-DRtemp)
            PM_DR=(PMtemp-DRtemp)
            egmax=max(egmax,np.max(abs(DRtemp)))
            clevels=np.arange(101)/100.0*2.0*egmax-egmax
            plt.title("DR")
            plt.pcolor(DRtemp,cmap=my_color_map(),vmin=-egmax,vmax=egmax)
            plt.xlabel('$k_z R$',size=18)
            plt.ylabel('$l$',size=18)
            plt.colorbar()
            plt.contour(NL_DR,levels,linewidths=3,cmap=black_color_map())
            plt.show()
            plt.title("DR (n="+str(n)+")")


        #plt.title("CD (n="+str(n)+")")
        ##print "min FE",np.min(FE[:,:,n])
        ##print "argmin FE",np.argmin(FE[:,:,n])
        ##plt.contour(CDtemp,clevels,cmap=my_color_map())
        #plt.pcolor(CDtemp,cmap=my_color_map(),vmin=-egmax,vmax=egmax)
        #plt.xlabel('$k_z R$',size=18)
        #plt.ylabel('$l$',size=18)
        #plt.colorbar()
        #plt.show()

        plt.title("PM (n="+str(n)+")")
        #temp=PM[:,:,n]/FE[:,:,n]
        #temp[:,par['nkz0']/2]=0.0
        #plt.pcolor(PMtemp,cmap=my_color_map(),vmin=-egmax,vmax=egmax)
        #print clevels
        print PMtemp
        plt.contourf(PMtemp,clevels,cmap=my_color_map())
        plt.xlabel('$k_z R$',size=18)
        plt.ylabel('$l$',size=18)
        plt.colorbar()
        plt.contour(NL_PM,levels,linewidths=3,cmap=black_color_map())
        plt.show()


        plt.title("NL (n="+str(n)+")")
        #temp=NL[:,:,n]/FE[:,:,n]
        #temp[:,par['nkz0']/2]=0.0
        plt.pcolor(NLtemp,cmap=my_color_map(),vmin=-egmax,vmax=egmax)
        plt.xlabel('$k_z R$',size=18)
        plt.ylabel('$l$',size=18)
        plt.colorbar()
        plt.show()


        #plt.plot(np.sum(CD[:,:,n],axis=1),'x-')
        #plt.show()
        #plt.title("Free Energy (n="+str(n)+")")
        #plt.pcolor(shell_grid,kzgrid,FE[:,:,n],cmap=my_color_map())
        #plt.plot(np.sum(FE[:,:,n],axis=1),'x-')
        #plt.show()


    #if nlt_comparison:
    #    if herm_n==0:
    #        fname='nlt_shells_n0.dat'
    #    elif herm_n==1:
    #        fname='nlt_shells_n1.dat'
    #    elif herm_n==2:
    #        fname='nlt_shells_n2.dat'
    #    elif herm_n==3:
    #        fname='nlt_shells_n3.dat'
    #    elif herm_n==par['nv0']/4:
    #        fname='nlt_shells_n1o4.dat'
    #    elif herm_n==par['nv0']/2:
    #        fname='nlt_shells_n2o4.dat'
    #    elif herm_n==3*par['nv0']/4:
    #        fname='nlt_shells_n3o4.dat'
    #    else:
    #        print "Invalid herm_n for nlt comparison."
    #        stop

    #    nlt,nlt_time,nlt_shells_t=nlt_test(file_name=fname)
    #    for i in range(num_shells):
    #        plt.plot(nlt_time,nlt_shells_t[i,:],'-x')
    #        plt.plot(ntime,ebal[:,3,i],'-x')
    #        plt.show()

def distribute_ebal_shells(es_in):
    num_shells=len(es_in)/(par['nv0']*par['nkz0'])
    #print num_shells
    es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
    ebal=np.empty(num_shells)
    for i in range(num_shells):
        ebal[i]=np.sum(es_in[:,:,i])
    return ebal

def distribute_ebal_shells_n(es_in,herm_n,keep_kz=False):
    num_shells=len(es_in)/(par['nv0']*par['nkz0'])
    #print num_shells
    es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))

    if keep_kz:
        ebal=np.empty((num_shells,par['nkz0']))
    else:
        ebal=np.empty(num_shells)

    for i in range(num_shells):
        if keep_kz:
            ebal[i,:]=es_in[herm_n,:,i]
        else:
            ebal[i]=np.sum(es_in[herm_n,:,i])
    return ebal

def distribute_ebal_shells_n_kz(es_in,herm_n,kz):
    num_shells=len(es_in)/(par['nv0']*par['nkz0'])
    #print num_shells
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    kzind=np.argmin(abs(kzgrid-kz))
    es_in=np.reshape(es_in,(par['nv0'],par['nkz0'],num_shells))
    ebal=np.empty(num_shells)
    for i in range(num_shells):
        ebal[i]=es_in[herm_n,kzind,i]
    return ebal

def nlt_test(start_time=-1,end_time=-1,file_name='nlt_shells.dat'):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    print "Number of kperp shells:", num_shells
    f.close()

    if not os.path.isfile(file_name):
        print "File does not exist!:", file_name
        stop
    
    f = open(file_name,'rb')
    ntot=num_shells**2*(par['nkz0']/2)**2
    mem_tot=ntot*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
          print input
      else:
          continue_read=0
 
    f.close()

    #print "time",time

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime 

    nlt_shells_t=np.zeros((num_shells,ntime))

    nlt_in=np.zeros(ntot)
    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        f = open(file_name,'rb')
        f.seek(8+i*(8+mem_tot))
        nlt_temp=np.fromfile(f,dtype='float64',count=ntot)
        nlt_in=nlt_in+nlt_temp
        nlt_temp=np.reshape(nlt_temp,(num_shells,num_shells,par['nkz0']/2,par['nkz0']/2),order='F')
        nlt_shells_t[:,i-istart]=np.sum(np.sum(np.sum(nlt_temp,axis=3),axis=2),axis=1)
        f.close()

    nlt_in=nlt_in/(float(ntime)) 
    nlt_in=np.reshape(nlt_in,(num_shells,num_shells,par['nkz0']/2,par['nkz0']/2),order='F')

    #nlt_kperp=np.sum(np.sum(nlt_in,axis=3),axis=2)
    #plt.contour(shell_grid,shell_grid,nlt_kperp,200)
    #plt.colorbar()
    #plt.show()
    nlt_kperp=np.sum(np.sum(nlt_in,axis=3),axis=2)
    #plt.pcolor(nlt_kperp,cmap='RdBu')
    #plt.pcolor(nlt_kperp,cmap='hot')
    shell_grid=np.arange(num_shells+1)+0.5
    plt.pcolor(shell_grid,shell_grid,nlt_kperp,cmap=my_color_map2())
    plt.xlabel('$l^\prime$',size=18)
    plt.ylabel('$l$',size=18)
    plt.colorbar()
    plt.show()

    nltkz=np.sum(np.sum(nlt_in,axis=0),axis=0)
    kz_grid=np.arange(par['nkz0']/2)/float(par['nkz0']/2-1)*par['kzmax0']-par['kzmin']/2.0
    plt.pcolor(kz_grid,kz_grid,nltkz,cmap=my_color_map2())
    plt.xlabel('$k_z^\prime$',size=18)
    plt.ylabel('$k_z$',size=18)
    plt.colorbar()
    plt.show()

    #for i in range(num_shells):
    #    plt.plot(time[istart:iend+1],nlt_shells_t[i,:])
    #    plt.show()

    return nlt_in, time[istart:iend+1],nlt_shells_t

def nlt_test_l_kz(start_time=-1,end_time=-1,file_name='nlt_shells.dat',kz_in=1,shell_in=0,show_plots=False):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    #print "Number of kperp shells:", num_shells
    f.close()

    kxg,kyg,kzg,hg=get_grids()

    if not os.path.isfile(file_name):
        print "File does not exist!:", file_name
        stop
    
    f = open(file_name,'rb')
    ntot=num_shells**2*(par['nkz0']/2)**2
    mem_tot=ntot*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
          #print input
      else:
          continue_read=0
 
    f.close()

    #print "time",time

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    #print "ntime",ntime 

    nlt_shells_t=np.zeros((num_shells,ntime))

    nlt_in=np.zeros(ntot)
    for i in range(istart,iend+1):
        #print 'time=',time[i],' of ',time[iend]
        f = open(file_name,'rb')
        f.seek(8+i*(8+mem_tot))
        nlt_temp=np.fromfile(f,dtype='float64',count=ntot)
        nlt_in=nlt_in+nlt_temp
        nlt_temp=np.reshape(nlt_temp,(num_shells,num_shells,par['nkz0']/2,par['nkz0']/2),order='F')
        nlt_shells_t[:,i-istart]=np.sum(np.sum(np.sum(nlt_temp,axis=3),axis=2),axis=1)
        f.close()

    nlt_in=nlt_in/(float(ntime)) 
    nlt_in=np.reshape(nlt_in,(num_shells,num_shells,par['nkz0']/2,par['nkz0']/2),order='F')

    shell_grid=np.arange(num_shells+1)+0.5
    kz_grid=np.arange(par['nkz0']/2)/float(par['nkz0']/2-1)*par['kzmax0']-par['kzmin']/2.0
    temp=nlt_in[shell_in,:,kz_in,:]
    nltmax=np.max(abs(temp))

    if show_plots:
        plt.pcolor(kz_grid,shell_grid,temp,cmap=my_color_map(),vmin=-nltmax,vmax=nltmax)
        plt.ylabel('$l^\prime$',size=18)
        plt.xlabel('$k_z^\prime $',size=18)
        plt.colorbar()
        plt.show()

    return temp

def nlt_out_l_kz(start_time=-1,end_time=-1,file_name='nlt_shells.dat',kz_in=1,shell_in=0,show_plots=False):
    f=open('shell_info.dat')
    s_info=f.read()
    s_lines=s_info.split('\n')
    num_shells=len(s_lines)-1 #int(float(s_lines[-2].split()[0]))
    #print "Number of kperp shells:", num_shells
    f.close()

    kxg,kyg,kzg,hg=get_grids()

    if not os.path.isfile(file_name):
        print "File does not exist!:", file_name
        stop
    
    f = open(file_name,'rb')
    ntot=num_shells**2*(par['nkz0']/2)**2
    mem_tot=ntot*8
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(i*(mem_tot+8))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
          #print input
      else:
          continue_read=0
 
    f.close()

    #print "time",time

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    #print "ntime",ntime 

    nlt_shells_t=np.zeros((num_shells,ntime))

    nlt_in=np.zeros(ntot)
    for i in range(istart,iend+1):
        #print 'time=',time[i],' of ',time[iend]
        f = open(file_name,'rb')
        f.seek(8+i*(8+mem_tot))
        nlt_temp=np.fromfile(f,dtype='float64',count=ntot)
        nlt_in=nlt_in+nlt_temp
        nlt_temp=np.reshape(nlt_temp,(num_shells,num_shells,par['nkz0']/2,par['nkz0']/2),order='F')
        nlt_shells_t[:,i-istart]=np.sum(np.sum(np.sum(nlt_temp,axis=3),axis=2),axis=1)
        f.close()

    nlt_in=nlt_in/(float(ntime)) 
    nlt_in=np.reshape(nlt_in,(num_shells,num_shells,par['nkz0']/2,par['nkz0']/2),order='F')

    shell_grid=np.arange(num_shells+1)+0.5
    kz_grid=np.arange(par['nkz0']/2)/float(par['nkz0']/2-1)*par['kzmax0']-par['kzmin']/2.0
    temp=nlt_in[shell_in,:,kz_in,:]
    nltmax=np.max(abs(temp))
    Tout=np.zeros(np.shape(temp))
    for i in range(np.shape(temp)[0]):
        for j in range(np.shape(temp)[1]):
            if temp[i,j]<0:
                Tout[i,j]=-temp[i,j]
    nout_kz=np.sum(Tout,axis=0)
    nout_shell=np.sum(Tout,axis=1)

    if show_plots:
        plt.pcolor(kz_grid,shell_grid,temp,cmap=my_color_map(),vmin=-nltmax,vmax=nltmax)
        plt.ylabel('$l^\prime$',size=18)
        plt.xlabel('$k_z^\prime $',size=18)
        plt.colorbar()
        plt.show()

        #kz0_trans=0.0
        #nltkz_sum=0.0
        #for i in range(len(nltkz)):
        #    if nltkz[i] < 0:
        #        kz0_trans+=(kzg[i]-kzg[kz_in])*nltkz[i]
        #        nltkz_sum+=nltkz[i]

        plt.plot(kz_grid,nout_kz)
        plt.xlabel('kz')
        plt.show()

        #shell_trans=0.0
        #nltshell_sum=0.0
        #for i in range(len(nltshell)):
        #    if nltshell[i] < 0:
        #        shell_trans+=(shell_grid[i]-shell_grid[shell_in])*nltshell[i]
        #        nltshell_sum+=nltshell[i]

        plt.plot(shell_grid[0:-1],nout_shell)
        plt.xlabel('kperp')
        plt.show()
        #if nltkz_sum==0.0:
        #    kz0_trans=0.0
        #else:
        #    kz0_trans=kz0_trans/nltkz_sum
        #if nltshell_sum==0.0:
        #    shell_trans=0.0
        #else:
        #    shell_trans=shell_trans/nltshell_sum

    nout_magnitude=np.sum(Tout)
    ninverse=0.0
    nforward=0.0
    for i in range(len(nout_kz)):
        if i < kz_in:
            ninverse-=nout_kz[i]
        if i > kz_in:
            nforward+=nout_kz[i]
    if show_plots:
        print "kz inverse:",ninverse
        print "kz forward:",nforward
    nout_kz0=ninverse+nforward

    ninverse=0.0
    nforward=0.0
    for i in range(len(nout_shell)):
        if i < shell_in:
            ninverse-=nout_shell[i]
        if i > shell_in:
            nforward+=nout_shell[i]
    if show_plots:
        print "kp inverse:",ninverse
        print "kp forward:",nforward
    nout_shell0=ninverse+nforward


    nout_kz0=nout_magnitude*nout_kz0/(nout_kz0**2+nout_shell0**2)**0.5
    nout_shell0=nout_magnitude*nout_shell0/(nout_kz0**2+nout_shell0**2)**0.5

    if show_plots:
        print "magnitude",nout_magnitude
        print "kz component",nout_kz0
        print "kperp component",nout_shell0

    return nout_magnitude,nout_shell0,nout_kz0

def my_color_map():
    """Returns a red-white-blue color map.  Great, e.g., for nonlinear transfer functions."""
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
           'green': ((0.0, 0.0, 0.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 0.0, 0.0)),
            'blue': ((0.0, 1.0, 1.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 0.0, 0.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    #plt.pcolor(rand(10,10),cmap=my_cmap)
    #plt.colorbar()
    #plt.show()
    return my_cmap

def my_color_map2():
    """Returns a red-white-blue color map.  Great, e.g., for nonlinear transfer functions."""
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.25, 0.5, 0.5),
                     (0.5, 1.0, 1.0),
                     (0.75, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
           'green': ((0.0, 0.0, 0.0),
                     (0.25, 0.5, 0.5),
                     (0.5, 1.0, 1.0),
                     (0.75, 0.5, 0.5),
                     (1.0, 0.0, 0.0)),
            'blue': ((0.0, 1.0, 1.0),
                     (0.25, 1.0, 1.0),
                     (0.5, 1.0, 1.0),
                     (0.75, 0.5, 0.5),
                     (1.0, 0.0, 0.0))}

    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    #plt.pcolor(rand(10,10),cmap=my_cmap)
    #plt.colorbar()
    #plt.show()
    return my_cmap

def les_comp(nkx0_red=32,nky0_red=64,nkz0_red=32):
    nb = int(raw_input('Introduce number of comparisons: '))
 
    diss0_kx = np.zeros((nkx0_red,2,nb))
    diss0_ky = np.zeros((nky0_red/2,2,nb))
    diss0_kz = np.zeros((nkz0_red/2,2,nb)) 
    drive0_kx = np.zeros((nkx0_red,2,nb))
    drive0_ky = np.zeros((nky0_red/2,2,nb))
    drive0_kz = np.zeros((nkz0_red/2,2,nb))
    entk0_kx = np.zeros((nkx0_red,2,nb))
    entk0_ky = np.zeros((nky0_red/2,2,nb))
    entk0_kz = np.zeros((nkz0_red/2,2,nb))
    hyp0_kx = np.zeros((nkx0_red,2,nb))
    hyp0_ky = np.zeros((nky0_red/2,2,nb))
    hyp0_kz = np.zeros((nkz0_red/2,2,nb))
    d_arr = ['r','b','g','y','m','c','k','r+','b+','g+','y+','m+','c+','k+'] 
    for i in range(nb):     
        diagdir = raw_input('Introduce the directories in order: ') 
        diss0_kx_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/diss0_kx_red.dat'
        diss0_kx[:,:,i] = np.loadtxt(diss0_kx_name)
        diss0_ky_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/diss0_ky_red.dat'
        diss0_ky[:,:,i] = np.loadtxt(diss0_ky_name)
        diss0_kz_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/diss0_kz_red.dat'
        diss0_kz[:,:,i] = np.loadtxt(diss0_kz_name)
        drive0_kx_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/drive0_kx_red.dat'
        drive0_kx[:,:,i] = np.loadtxt(drive0_kx_name)
        drive0_ky_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/drive0_ky_red.dat'
        drive0_ky[:,:,i] = np.loadtxt(drive0_ky_name)
        drive0_kz_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/drive0_kz_red.dat'
        drive0_kz[:,:,i] = np.loadtxt(drive0_kz_name)
        entk0_kx_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/entk0_kx_red.dat'
        entk0_kx[:,:,i] = np.loadtxt(entk0_kx_name)
        entk0_ky_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/entk0_ky_red.dat'
        entk0_ky[:,:,i] = np.loadtxt(entk0_ky_name)
        entk0_kz_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/entk0_kz_red.dat'
        entk0_kz[:,:,i] = np.loadtxt(entk0_kz_name)
        hyp0_kx_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/hyp0_kx_red.dat'
        hyp0_kx[:,:,i] = np.loadtxt(hyp0_kx_name)
        hyp0_ky_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/hyp0_ky_red.dat'
        hyp0_ky[:,:,i] = np.loadtxt(hyp0_ky_name)
        hyp0_kz_name = '/gpfs/ipp/abanonna/dna/' + diagdir + '/hyp0_kz_red.dat'
        hyp0_kz[:,:,i] = np.loadtxt(hyp0_kz_name)
  
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        #log-log plot
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(drive0_kx[:,0,i],drive0_kx[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Drive')
        plt.xlabel(r'$k_x\rho_i$',size=13)
    plt.show()

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(drive0_kx[:,0,i],drive0_kx[:,1,i],d_arr[i])
        plt.xlabel(r'$k_x\rho_i$',size=13)
        plt.title('Drive')
    plt.show()

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        #log-log plot
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(diss0_kx[:,0,i],-1*diss0_kx[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Dissipation')
        plt.xlabel(r'$k_x\rho_i$',size=13)
    plt.show()

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(diss0_kx[:,0,i],diss0_kx[:,1,i],d_arr[i])
        plt.xlabel(r'$k_x\rho_i$',size=13)
        plt.title('Dissipation')
    plt.show()

    
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        #log-log plot
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(entk0_kx[:,0,i],entk0_kx[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Free energy')
        plt.xlabel(r'$k_x\rho_i$',size=13)
    plt.show()
 
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(entk0_kx[:,0,i],entk0_kx[:,1,i],d_arr[i])
        plt.xlabel(r'$k_x\rho_i$',size=13)
        plt.title('Free energy')
    plt.show()


    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        #log-log plot
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(hyp0_kx[:,0,i],-1*hyp0_kx[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Hyperdiffusion')
        plt.xlabel(r'$k_x\rho_i$',size=13)
    plt.show()
 
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(hyp0_kx[:,0,i],hyp0_kx[:,1,i],d_arr[i])
        plt.xlabel(r'$k_x\rho_i$',size=13)
        plt.title('Hyperdiffusion')
    plt.show()

   
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb):  
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(drive0_ky[:,0,i],drive0_ky[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Drive')
        plt.xlabel(r'$k_y\rho_i$',size=13)
    plt.show()
 

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(drive0_ky[:,0,i],drive0_ky[:,1,i],d_arr[i])
        plt.xlabel(r'$k_y\rho_i$',size=13)
        plt.title('Drive')
    plt.show()


    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb):  
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(diss0_ky[:,0,i],-1*diss0_ky[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Dissipation')
        plt.xlabel(r'$k_y\rho_i$',size=13)
    plt.show()
 
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(diss0_ky[:,0,i],diss0_ky[:,1,i],d_arr[i])
        plt.xlabel(r'$k_y\rho_i$',size=13)
        plt.title('Dissipation')
    plt.show()


    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb):  
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(entk0_ky[:,0,i],entk0_ky[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Free energy')
        plt.xlabel(r'$k_y\rho_i$',size=13)
    plt.show()
 
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(entk0_ky[:,0,i],entk0_ky[:,1,i],d_arr[i])
        plt.xlabel(r'$k_y\rho_i$',size=13)
        plt.title('Free energy')
    plt.show()


    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb):  
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(hyp0_ky[:,0,i],-1*hyp0_ky[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Hyperdiffusion')
        plt.xlabel(r'$k_y\rho_i$',size=13)
    plt.show()
 

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(hyp0_ky[:,0,i],hyp0_ky[:,1,i],d_arr[i])
        plt.xlabel(r'$k_y\rho_i$',size=13)
        plt.title('Hyperdiffusion')
    plt.show()


    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb):     
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(drive0_kz[:,0,i],drive0_kz[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Drive')
        plt.xlabel(r'$k_z R$',size=13)
    plt.show()

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(drive0_kz[:,0,i],drive0_kz[:,1,i],d_arr[i])
        plt.xlabel(r'$k_z\rho_i$',size=13)
        plt.title('Drive')
    plt.show()

    fig=plt.gcf()
    for i in range(nb):     
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(diss0_kz[:,0,i],-1*diss0_kz[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Dissipation')
        plt.xlabel(r'$k_z R$',size=13)
    plt.show()

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(diss0_kz[:,0,i],diss0_kz[:,1,i],d_arr[i])
        plt.xlabel(r'$k_z\rho_i$',size=13)
        plt.title('Dissipation')
    plt.show()


    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb):     
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(entk0_kz[:,0,i],entk0_kz[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Free energy')
        plt.xlabel(r'$k_z R$',size=13)
    plt.show()

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(entk0_kz[:,0,i],entk0_kz[:,1,i],d_arr[i])
        plt.xlabel(r'$k_z\rho_i$',size=13)
        plt.title('Free energy')
    plt.show()

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb):     
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(hyp0_kz[:,0,i],-1*hyp0_kz[:,1,i],d_arr[i],basex=10,basey=10)
        plt.title('Hyperdifusion')
        plt.xlabel(r'$k_z R$',size=13)
    plt.show()

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    for i in range(nb): 
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.2)
        plt.plot(hyp0_kz[:,0,i],hyp0_kz[:,1,i],d_arr[i])
        plt.xlabel(r'$k_z\rho_i$',size=13)
        plt.title('Hyperdiffusion')
    plt.show()

def fit_function(x,y,p0,startx=-1,endx=-1,which_func='power'):
    if startx==-1:
        sx0=0
    else:
        sx0=argmin(abs(x-startx))
    if endx==-1:
        ex0=len(x)
    else:
        ex0=argmin(abs(x-endx))

    if ex0 < sx0:
        stop

    x0=x[sx0:ex0+1]
    y0=y[sx0:ex0+1]
    print "sx0,ex0",sx0,ex0
    print "p0",p0

    if which_func=='power':
        p,success=scipy.optimize.leastsq(power_law_res_func,p0,args=(x0,y0))
    elif which_func=='power_exp':
        p,success=scipy.optimize.leastsq(power_law_exp_res_func,p0[2:],args=(x0,y0,p0[0:2]))
    else:
        print "invalid selection for which_func!"
        stop

    if success < 1 or success > 4:
        print "Error: fit failed to converge!"
        stop
    else:
        print "Fit result:", p
        return p

def power_law_res_func(pars,x,y):
    res=pars[0]*x**(pars[1])-y
    return res

def power_law_exp_res_func(p_in,x,y,p_pow):
    res=p_pow[0]*x**(p_pow[1])*np.e**(p_in[0]*x**p_in[1])-y
    #res=pars[0]*x**(pars[1])*np.e**(pars[2]*x**1.5)-y
    return res

def test_landau_vs_nl(ikx0=0,start_time=-1.0,end_time=-1.0,num_ky=-1,num_kz=-1,start_ky=1,start_kz=0,plot_time=False):
    """Routine for testing Landau damping vs nonlinear, etc. a la Plunk"""

    start_kz=0
    num_kz=par['nkz0']/2-1
    out_dir=par['diagdir'][1:-1]+'/landau_test'
    os.system('mkdir '+out_dir)
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    allk=False
    if allk:
        print "allk option not ready at this time."
    else:
        kx=kxgrid[ikx0]
        if num_ky==-1:
            num_ky=par['nky0']
        if num_kz==-1:
            num_kz=par['nkz0']

        time=get_time_from_gout()
        print "getting time from fmom3d"
        ptime=get_time_from_fmom3d()
        print "Done."

        if start_time==-1.0:
            start_time=time[0]
        if end_time==-1.0:
            end_time=time[len(time)-1]
        if start_time > end_time:
            stop

        istart=np.argmin(abs(time-start_time))
        iend=np.argmin(abs(time-end_time))
        ntime=iend-istart+1
        
        ipstart=np.argmin(abs(ptime-start_time))
        ipend=np.argmin(abs(ptime-end_time))
        nptime=ipend-ipstart+1

        print "ntime",ntime
        print "nptime",nptime
        print num_ky*num_kz*par['nv0']*ntime*16/1000000.0, "MB"
        print num_ky*num_kz*nptime*16/1000000.0, "MB"

        time0=np.zeros(ntime)
        ptime0=np.zeros(nptime)
        gkykz=np.empty((num_ky,num_kz,par['nv0'],ntime),dtype='complex')
        phikykz=np.empty((num_ky,num_kz,nptime),dtype='complex')
        kz0_ky=np.empty((num_ky),dtype='float')
        gamma_L=np.empty((num_ky,num_kz))
        gamma_lin=np.empty((num_ky,num_kz))
        gamma_Llin=np.empty((num_ky,num_kz))
        omega_nl=np.empty((num_ky,num_kz))

        for t in range(ipstart,ipend+1):
            tind=t-ipstart
            print 'time=',ptime[t],' of ',ptime[ipend]
            phi=read_phi_fmom3d(t)
            phi=np.reshape(phi,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
            phikykz[:,:,tind]=phi[ikx0,start_ky:start_ky+num_ky,start_kz:start_kz+num_kz]
            ptime0[tind]=ptime[t]

        phisq=np.sum(np.conj(phikykz)*phikykz,axis=2)
        for j in range(start_ky,start_ky+num_ky):
            kz0_ky[j]=np.sum(kzgrid[0:num_kz]*phisq[j,:])/np.sum(phisq[j,:])

        for t in range(istart,iend+1):
            tind=t-istart
            print 'time=',time[t],' of ',time[iend]
            gt0=read_time_step_g(t)
            gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
            #eterm=get_henergy_single_k(gt0,gt0,kx,ky,kz,9)
            #nlterm[tind]=eterm[0]
            gkykz[:,:,:,tind]=gt0[ikx0,start_ky:start_ky+num_ky,start_kz:start_kz+num_kz,:]
            time0[tind]=time[t]

        #GsinkL =np.zeros((num_ky,num_kz))
        phi=np.zeros(nptime,dtype='complex128')
        phasemix=np.zeros(ntime)
        omnterm=np.zeros(ntime)
        omtflr=np.zeros(ntime)
        fenergy=np.zeros(ntime)
        Erhs=np.zeros(ntime)

        for j in range(start_ky,start_ky+num_ky):
            for k in range(start_kz,start_kz+num_kz):
                j0=j-start_ky
                k0=k-start_kz
                ky=kygrid[j]
                kz=kzgrid[k]
                print "kx,ky,kz",kx,ky,kz
                #svdo = get_svd_spectrum(0.0,ky,kz,start_time=start_time,end_time=end_time,\
                #                        calc_evs=False,read_g=False,g_in=gkykz[j,k,:,:])
                print kx,ky,kz
                for t in range(nptime):
                    phi[t]=phikykz[j0,k0,t]
                for t in range(ntime):
                    eterm=get_henergy_single_k(gkykz[j0,k0,:,t],gkykz[j0,k0,:,t],kx,ky,kz,0)
                    Erhs[t]=eterm[0]
                    eterm=get_henergy_single_k(gkykz[j0,k0,:,t],gkykz[j0,k0,:,t],kx,ky,kz,3)
                    phasemix[t]=eterm[0]
                    eterm=get_henergy_single_k(gkykz[j0,k0,:,t],gkykz[j0,k0,:,t],kx,ky,kz,4)
                    omnterm[t]=eterm[0]
                    eterm=get_henergy_single_k(gkykz[j0,k0,:,t],gkykz[j0,k0,:,t],kx,ky,kz,7)
                    omtflr[t]=eterm[0]
                    eterm=get_henergy_single_k(gkykz[j0,k0,:,t],gkykz[j0,k0,:,t],kx,ky,kz,-1)
                    fenergy[t]=eterm[0]
                #np.savetxt(par['diagdir'][1:-1]+'/phi_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz),phi)
                ###########################
                #Calculate correlation time
                ###########################
                ########Temp#########
                phicorr,tau,corr_time=my_corr_func_complex(phi,phi,ptime0,show_plot=False)
                #corr_time=1.0
                ########Temp#########
                #temp=(fenergy-np.roll(fenergy,1))
                #temp[0]=0.0
                #temp_dt=time0-np.roll(time0,1)
                #temp_dt[0]=0.0
                #dt=np.sum(temp_dt)/float((len(time0)-1))
                #temp=temp/dt
                #plt.plot(time0,Erhs,'+-',label='Erhs')
                #plt.plot(time0,temp,label='dEdt')
                #plt.plot(time0,fenergy,label='E')
                #plt.plot(time0,phasemix+omnterm+omtflr,'x-',label='Lsum')
                #plt.legend()
                #plt.title('Erhs vs dEdt')
                #plt.show()
                print "correlation time:", corr_time
                omega_nl[j0,k0]=1.0/corr_time
                ###########################
                #Calculate gamma_L
                ###########################
                #Note: 0.5 factor due to energy terms being quadratic
                temp=phasemix/fenergy
                gamma_L[j0,k0]=0.5*np.sum(temp)/len(temp)
                gammaLt=0.5*temp
                ###########################
                #Calculate gamma_lin
                ###########################
                mat=get_lin_matrix(kx,ky,kz)
                om,gm,ev=get_ev_spectrum(mat)
                gamma_lin[j0,k0]=gm[0]
                ###########################
                #Calculate gamma_Llin
                ###########################
                #Note: 0.5 factor due to energy terms being quadratic
                pm_temp=get_henergy_single_k(ev[:,0],ev[:,0],ikx0,ky,kz,3)
                fe_temp=get_henergy_single_k(ev[:,0],ev[:,0],ikx0,ky,kz,-1)
                gamma_Llin[j0,k0]=0.5*pm_temp[0]/fe_temp[0]
                if plot_time:
                    print "Plotting time dependence."
                    plt.title('ky='+str(ky)+' kz='+str(kz)+' kx='+str(kx))
                    gLlin=(np.zeros(len(time0))+1)*gamma_Llin[j0,k0]
                    glin=(np.zeros(len(time0))+1)*gm[0]
                    plt.plot(time0,gammaLt,'x-',label=r'$\gamma_{L,nl}$')
                    plt.plot(time0,gLlin,'+-',label=r'$\gamma_{L,lin}$')
                    plt.plot(time0,glin,label=r'$\gamma_{lin}$')
                    #plt.plot(time0,omnterm,label='OMN')
                    #plt.plot(time0,omtflr,label='FLR')
                    plt.legend()
                    plt.show()


        kzout=kzgrid[start_kz:start_kz+num_kz]
        for j in range(start_ky,start_ky+num_ky):
            j0=j-start_ky
            ky=kygrid[j]
            plt.title(r'$k_y \rho_i = $'+str(ky)[:4]+' '+ r'$k_x \rho_i = $'+str(kx)[:4],size=18)
            plt.plot(kzout,gamma_L[j0,:],'x-',label=r'$\gamma_{L,nl}$')
            ########Temp#########
            #plt.plot(kzout,omega_nl[j0,:],'+-',label=r'$\omega_{nl}$')
            ########Temp#########
            plt.plot(kzout,gamma_lin[j0,:],label=r'$\gamma_{lin}$')
            plt.plot(kzout,gamma_Llin[j0,:],label=r'$\gamma_{L,lin}$')
            plt.plot(kzout,kzout,label=r'$k_z v_{ti}$')
            ax=plt.axis()
            plt.vlines(kz0_ky,ax[2],ax[3])
            plt.xlabel(r'$k_z R$',size=18)
            plt.ylabel(r'$\omega (v_{ti}/R)$',size=18)
            plt.legend()
            plt.show()

                

                
                #np.savetxt(par['diagdir'][1:-1]+'/landau_test'+'/phasemix_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz),phasemix)
                #np.savetxt(par['diagdir'][1:-1]+'/landau_test'+'/omnterm_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz),omnterm)
                #np.savetxt(par['diagdir'][1:-1]+'/landau_test'+'/omtflr_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz),omtflr)
                #np.savetxt(par['diagdir'][1:-1]+'/landau_test'+'/fenergy_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz),fenergy)
                #np.savetxt(par['diagdir'][1:-1]+'/landau_test'+'/nlterm_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz),nlterm)

def test_mode_growth():
    
    time=get_time_from_fmom3d()
    phi2kznot0=np.zeros(len(time))
    phi2kx0=np.zeros(len(time))
    phi2ky0=np.zeros(len(time))
    phi2kz0=np.zeros(len(time))
    phi2_nozeros=np.zeros(len(time))
    phi2kx0only=np.zeros(len(time))
    phi2ky0only=np.zeros(len(time))
    phi2kz0only=np.zeros(len(time))
    phi2tot=np.zeros(len(time))
    phi2allzero=np.zeros(len(time))
    phi2kx0kz0=np.zeros(len(time))
    phi2ky0kz0=np.zeros(len(time))
    phi2kx0ky0=np.zeros(len(time))
    phi2kz0_2d=np.zeros((par['nkx0'],par['nky0'],len(time)))
    phi2kxsum=np.zeros((par['nky0'],par['nkz0'],len(time)))
    phi2zf=np.zeros((par['nkx0'],len(time)))
    kxgrid,kygrid,kzgrid,hg=get_grids()
    for i in range(len(time)):
        print i, ' of ',len(time)
        phi=read_phi_fmom3d(i)
        phi=np.reshape(phi,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        phisq=np.real(np.conj(phi)*phi)
        phi2zf[:,i]=phisq[:,0,0]
        phi2kxsum[:,:,i]=np.sum(phisq[:,:,:],axis=0)
        phi2kz0_2d[:,:,i]=phisq[:,:,0]
        #print np.shape(phi2kz0_2d)
        #if i > 2:
        #    plt.contour(phi2kz0_2d[:,:,i],50)
        #plt.contour(kygrid,kxgrid,phi2kz0_2d[:,:,i])
        #plt.contour(kygrid,kxgrid,phi2kz0[:,:,i])
        plt.show()
        phi2tot[i]=np.sum(phisq)
        phi2kx0[i]=np.sum(phisq[0,:,:])
        phi2ky0[i]=np.sum(phisq[:,0,:])
        phi2kz0[i]=np.sum(phisq[:,:,0])
        for ikx in range(par['nkx0']):
            for iky in range(par['nky0']):
                for ikz in range(par['nkz0']):
                    if ikz!=0:
                        phi2kznot0[i]+=phisq[ikx,iky,ikz]
                    if ikx !=0 and iky !=0 and ikz !=0:
                        phi2_nozeros[i]+=phisq[ikx,iky,ikz]
                    if ikx==0 and iky !=0 and ikz !=0:
                        phi2kx0only[i]+=phisq[ikx,iky,ikz]
                    if ikx!=0 and iky ==0 and ikz !=0:
                        phi2ky0only[i]+=phisq[ikx,iky,ikz]
                    if ikx!=0 and iky !=0 and ikz ==0:
                        phi2kz0only[i]+=phisq[ikx,iky,ikz]
                    if ikx==0 and iky !=0 and ikz ==0:
                        phi2kx0kz0[i]+=phisq[ikx,iky,ikz]
                    if iky ==0 and ikz ==0:
                        phi2ky0kz0[i]+=phisq[ikx,iky,ikz]
                    if ikx==0 and iky ==0 and ikz !=0:
                        phi2kx0ky0[i]+=phisq[ikx,iky,ikz]

    plt.contour(time,kxgrid,phi2zf,50)
    plt.xlabel('time')
    plt.ylabel('kx')
    plt.show()

    plt.plot(time,phi2kznot0,label='phi2kznot0')
    plt.plot(time,phi2tot,label='phi2tot')
    plt.plot(time,phi2kz0,label='phi2kz0')
    plt.plot(time,phi2kz0+phi2kznot0,'x',label='phi2kz0')
    plt.xlabel('t')
    plt.legend()
    plt.show()

    phi2kxsum_kztime=np.sum(phi2kxsum,axis=0)
    phi2kxsum_kytime=np.sum(phi2kxsum,axis=1)
    plt.contour(time,kygrid,phi2kxsum_kytime,50)
    plt.title('kxsum')
    plt.xlabel('time')
    plt.ylabel('ky')
    plt.show()
    plt.contour(time,kzgrid,phi2kxsum_kztime,50)
    plt.title('kxsum')
    plt.xlabel('time')
    plt.ylabel('kz')
    plt.show()

    phi2kz0_kytime=np.sum(phi2kz0_2d,axis=0)
    phi2kz0_kxtime=np.sum(phi2kz0_2d,axis=1)
    #plt.contour(time,kygrid,phi2kz0_kytime,50)
    plt.contour(time,kygrid,phi2kz0_kytime,50)
    plt.xlabel('time')
    plt.ylabel('ky')
    plt.show()

    plt.contour(time,kxgrid,phi2kz0_kxtime,50)
    plt.xlabel('time')
    plt.ylabel('kx')
    plt.show()

    plt.plot(time,phi2ky0kz0,label='phi2ky0kz0')
    plt.plot(time,phi2kx0kz0,label='phi2kx0kz0')
    plt.plot(time,phi2kz0,label='phi2kz0')
    plt.plot(time,phi2ky0only,label='phi2kz0only')
    plt.xlabel('t')
    plt.legend()
    plt.show()


    plt.plot(time,phi2kx0ky0,label='phi2kx0ky0')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2ky0kz0,label='phi2ky0kz0')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2kx0kz0,label='phi2kx0kz0')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2allzero,label='phi2allzero')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2tot,label='phisq_tot')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2kx0,label='phi2kx0')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2ky0,label='phi2ky0')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2kz0,label='phi2kz0')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2_nozeros,label='phi2_nozeros')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2kx0only,label='phi2kx0only')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2ky0only,label='phi2ky0only')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    plt.plot(time,phi2ky0only,label='phi2kz0only')
    plt.xlabel('t')
    plt.legend()
    plt.show()

    for i in range(par['nky0']):
        plt.title('ky='+str(i*par['kymin']))
        plt.contour(time,kxgrid,phi2kz0_2d[:,i,:])
        plt.xlabel('time')
        plt.ylabel('kx')
        plt.show()






def test_all_rhs_terms(kx,ky,kz,herm_n,start_time=-1.0,end_time=-1.0,calc_nl=True):
    """Routine for testing all the RHS terms"""

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    time=get_time_from_gout()
    ikx=get_kindex(kx,par['kxmin'],par['nkx0'])
    iky=get_kindex(ky,par['kymin'],par['nky0'])
    ikz=get_kindex(kz,par['kzmin'],par['nkz0'])

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    time2=get_time_from_energy()
    istart2=np.argmin(abs(time2-start_time))
    iend2=np.argmin(abs(time2-end_time))
    ntime2=iend2-istart2+1

    time0=np.zeros(ntime)
    phasemix=np.zeros(ntime)
    omnterm=np.zeros(ntime)
    omtflr=np.zeros(ntime)
    omtterm=np.zeros(ntime)
    fenergy=np.zeros(ntime)
    diss=np.zeros(ntime)
    Erhs=np.zeros(ntime)
    kzphi=np.zeros(ntime)
    nlterm=np.zeros(ntime)
    drive=np.zeros(ntime2)
    fenergy2=np.zeros(ntime2)
    diss2=np.zeros(ntime2)
    time02=np.zeros(ntime2)
    
    ntot=par['nkx0']*par['nky0']*par['nkz0']
    mem_tot=ntot*8

    for t in range(istart2,iend2+1):


        tind=t-istart2

        time02[tind]=time2[t]
        f = open('energy3d.dat','rb')
        f.seek(t*(8+4*mem_tot))
        ttemp=np.fromfile(f,dtype='float64',count=1)
        print ttemp
        f.seek(8+t*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        fenergy2[tind]=en_in[ikx,iky,ikz]
        f.seek(8+mem_tot+t*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        drive[tind]=en_in[ikx,iky,ikz]
        f.seek(8+2*mem_tot+t*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        diss2[tind]=diss2[tind]+en_in[ikx,iky,ikz]
        f.seek(8+3*mem_tot+t*(8+4*mem_tot))
        en_in=np.fromfile(f,dtype='float64',count=ntot)
        en_in=np.reshape(en_in,(par['nkx0'],par['nky0'],par['nkz0']),order='F')
        diss2[tind]=diss2[tind]+en_in[ikx,iky,ikz]
        f.close()


    for t in range(istart,iend+1):
        tind=t-istart
        print 'time=',time[t],' of ',time[iend]
        gt0=read_time_step_g(t)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        time0[tind]=time[t]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,0)
        if herm_n==-1:
            Erhs[tind]=np.sum(eterm)
        else:
            Erhs[tind]=eterm[herm_n]
        print 'Erhs',Erhs[tind]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,3)
        if herm_n==-1:
            phasemix[tind]=np.sum(eterm)
        else:
            phasemix[tind]=eterm[herm_n]
        print 'phasemix',phasemix[tind]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,4)
        if herm_n==-1:
            omnterm[tind]=np.sum(eterm)
        else:
            omnterm[tind]=eterm[herm_n]
        print 'omnterm',omnterm[tind]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,7)
        if herm_n==-1:
            omtflr[tind]=np.sum(eterm)
        else:
            omtflr[tind]=eterm[herm_n]
        print 'omtflr',omtflr[tind]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,-1)
        if herm_n==-1:
            fenergy[tind]=np.sum(eterm)
        else:
            fenergy[tind]=eterm[herm_n]
        print 'fenergy',fenergy[tind]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,6)
        if herm_n==-1:
            omtterm[tind]=np.sum(eterm)
        else:
            omtterm[tind]=eterm[herm_n]
        print 'omtterm',omtterm[tind]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,1)
        if herm_n==-1:
            diss[tind]=diss[tind]+np.sum(eterm)
        else:
            diss[tind]=diss[tind]+eterm[herm_n]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,2)
        if herm_n==-1:
            diss[tind]=diss[tind]+np.sum(eterm)
        else:
            diss[tind]=diss[tind]+eterm[herm_n]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,8)
        if herm_n==-1:
            diss[tind]=diss[tind]+np.sum(eterm)
        else:
            diss[tind]=diss[tind]+eterm[herm_n]
        print 'diss',diss[tind]
        eterm=get_henergy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],kx,ky,kz,5)
        if herm_n==-1:
            kzphi[tind]=np.sum(eterm)
        else:
            kzphi[tind]=eterm[herm_n]
        print 'kzphi',kzphi[tind]
        if calc_nl:
          eterm=get_henergy_single_k(gt0,gt0,kx,ky,kz,9)
          if herm_n==-1:
              nlterm[tind]=np.sum(eterm)
          else:
              nlterm[tind]=eterm[herm_n]
        else:
          eterm=0.0
          nlterm[tind]=0.0
        print 'nlterm',nlterm[tind]

    dEdt=fenergy-np.roll(fenergy,1)
    dEdt[0]=0
    dt=time0-np.roll(time0,1)
    dEdt=dEdt/dt
    if herm_n==-1:
        plt.title("Energetics, all n")
    else:
        plt.title("Energetics, n="+str(herm_n))
    plt.plot(time0,omtterm,label='Q')
    plt.plot(time0,diss,label='C')
    plt.plot(time0,phasemix+kzphi,label='PM')
    if calc_nl:
        plt.plot(time0,nlterm,label='NL')
    #plt.plot(time0,dEdt,'x-',label='dEdt')
    #plt.plot(time0,omtterm+diss+phasemix+nlterm+kzphi+omnterm+omtflr,'x-',color='black',label='Sum')
    plt.plot(time0,dEdt-Erhs,'o-',color='black',label='Inferred NL')
    #plt.plot(time0,Erhs,'x-',label='Rlin')
    #plt.plot(time0,fenergy,'-',label='FE')
    #plt.plot(time0,kzphi,'--',label='L')
    #plt.plot(time0,omnterm,'-.',label='OMN')
    #plt.plot(time0,omtflr,':',label='OMTFLR')
    #plt.plot(time0,omtterm+diss+phasemix+nlterm+kzphi+omnterm+omtflr,'x-',color='black',label='Sum')
    #plt.plot(time0,omtterm+diss+phasemix+kzphi+omnterm+omtflr,'+-',color='red',label='Sum lin')
    #plt.plot(time02,fenergy2,'<',label='FE2')
    #plt.plot(time02,drive,'>',label='Q2')
    #plt.plot(time02,diss2,'|',label='C2')
    plt.legend()
    plt.show()

def get_rhs_nl_single_k(g_in,kx,ky,kz,verbose=False):

    ###print "!!!!Warning!!!! Not yet benchmarked!!!"
    ###Successfully benchmarked on 5/3/13
    phib=get_phi_bar(g_in)
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    rhsnl=np.zeros(par['nv0'])
    if verbose:
        print "Calculating nonlinear term."
    ikxprange=np.arange(par['nkx0']*2-1)-(par['nkx0']-1)
    #print ikxprange
    for i in ikxprange:
        if verbose:
            print i
        for j in range(par['nky0']):
            for k in range(par['nkz0']):
                kxp=i*par['kxmin']
                kyp=kygrid[j]
                kzp=kzgrid[k]
                kxpp=kx-kxp
                kypp=ky-kyp
                kzpp=kz-kzp
                if np.abs(np.rint(kxp/par['kxmin'])) <= np.rint(par['kxmax0']/par['kxmin']) and \
                   np.abs(np.rint(kyp/par['kymin'])) <= np.rint(par['kymax0']/par['kymin']) and \
                   np.abs(np.rint(kzp/par['kzmin'])) <= np.rint(par['kzmax0']/par['kzmin']) and \
                   np.abs(np.rint(kxpp/par['kxmin'])) <= np.rint(par['kxmax0']/par['kxmin']) and \
                   np.abs(np.rint(kypp/par['kymin'])) <= np.rint(par['kymax0']/par['kymin']) and \
                   np.abs(np.rint(kzpp/par['kzmin'])) <= np.rint(par['kzmax0']/par['kzmin']):

                    ip,jp,kp,conjg_p=get_all_indices(kxp,kyp,kzp)
                    ipp,jpp,kpp,conjg_pp=get_all_indices(kxpp,kypp,kzpp)
                    if conjg_p:
                        phib0=np.conj(phib[ip,jp,kp])
                    else:
                        phib0=phib[ip,jp,kp]
                    if conjg_pp:
                        g0=np.conj(g_in[ipp,jpp,kpp,:])
                    else:
                        g0=g_in[ipp,jpp,kpp,:]
                    rhsnl=rhsnl+(kxp*ky-kx*kyp)*phib0*g0
    #print "done"
    return rhsnl

#def get_nlt_single_k(g_in,kx,ky,kz):
#    phib=get_phi_bar(g_in)
#    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
#    kxgrid_nlt,kygrid_nlt,kzgrid_nlt=get_grids_nlt()
#    #print "kxgrid_nlt",kxgrid_nlt
#    #print "kygrid_nlt",kygrid_nlt
#    #print "kzgrid_nlt",kzgrid_nlt
#    nkx1=(2*par['nkx0']-1)
#    nky1=(par['nky0']-1)
#    nkz1=(par['nkz0']-1)
#    NLT=np.zeros((nkx1,nky1,nkz1))
#    #kx=kxgrid[ikx]
#    #ky=kygrid[iky]
#    #kz=kzgrid[ikz]
#    phibk=phib[ikx,iky,ikz]
#    gk=g_in[ikx,iky,ikz,:]
#    #The following loops over k prime
#    for i in range(nkx1):
#        print i, " of ", nkx1-1
#        for j in range(nky1):
#            for k in range(nkz1):
#                kxp=kxgrid_nlt[i]
#                kyp=kygrid_nlt[j]
#                kzp=kzgrid_nlt[k]
#                kxpp=kx-kxp
#                kypp=ky-kyp
#                kzpp=kz-kzp
#                if np.abs(kxp) < par['kxmax0'] and \
#                   np.abs(kyp) < par['kymax0'] and \
#                   np.abs(kzp) < par['kzmax0'] and \
#                   np.abs(kxpp) < par['kxmax0'] and \
#                   np.abs(kypp) < par['kymax0'] and \
#                   np.abs(kzpp) < par['kzmax0']:
#
#                    ikxp,ikyp,ikzp,take_conjg_p=get_all_indices(kxp,kyp,kzp)
#                    #print "kxp,kyp,kzp",kxp,kyp,kzp
#                    #print "ikxp,ikyp,ikzp,take_conjg",ikxp,ikyp,ikzp,take_conjg_p
#                    ikxpp,ikypp,ikzpp,take_conjg_pp=get_all_indices(kxpp,kypp,kzpp)
#                    #print "kxpp,kypp,kzpp",kxpp,kypp,kzpp
#                    #print "ikxpp,ikypp,ikzpp,take_conjg",ikxpp,ikypp,ikzpp,take_conjg_pp
#                    ckkp=kxp*ky-kx*kyp
#                    if take_conjg_p:
#                        phib_p=np.conj(phib[ikxp,ikyp,ikzp])
#                        g_p=np.conj(g_in[ikxp,ikyp,ikzp,:])
#                    else:
#                        phib_p=(phib[ikxp,ikyp,ikzp])
#                        g_p=(g_in[ikxp,ikyp,ikzp,:])
#                    if take_conjg_pp:
#                        phib_pp=np.conj(phib[ikxpp,ikypp,ikzpp])
#                        g_pp=np.conj(g_in[ikxpp,ikypp,ikzpp,:])
#                    else:
#                        phib_pp=(phib[ikxpp,ikypp,ikzpp])
#                        g_pp=(g_in[ikxpp,ikypp,ikzpp,:])
#
#
#                    NLT[i,j,k]=np.real(np.pi**0.25*ckkp*np.conj(phibk)*\
#                        phib_p*g_pp[0])
#                    for n in range(par['nv0']):
#                        NLT[i,j,k]=NLT[i,j,k]-np.real(np.pi**0.5*ckkp*np.conj(gk[n])*\
#                            phib_pp*g_p[n])
#                    #NLT[i,j,k]=np.pi**0.25*ckkp*np.conj(phibk)*\
#                    #    phib[ikxp,ikyp,ikzp]*g_in[ikxpp,ikypp,ikzpp,0]
#                    #for n in range(par['nv0']):
#                    #    NLT[i,j,k]=NLT[i,j,k]-np.pi**0.5*ckkp*np.conj(gk[n])*\
#                    #        phib[ikxpp,ikypp,ikzpp]*g_in[ikxp,ikyp,ikzp,n]
#    return NLT
#

def get_lin_matrix(kx,ky,kz,verbose=False,test_extra=False,zamp=1.0):
    mat=np.zeros((par['nv0'],par['nv0']),dtype='complex64')
    g_1=np.zeros(par['nv0'],dtype='complex64')
    for i in range(par['nv0']):
        if verbose:
            print i
        g_1[:]=0.0        
        g_1[i]=1.0+0J
        mat[:,i]=get_rhs_lin_single_k(g_1,0,kx,ky,kz,test_extra=test_extra,zamp=zamp)
    return mat 

def get_ev_spectrum(mat,show_plots=False,solver=1,plot_evecs=False,num_evecs_plot=8):
    evecs=np.empty((par['nv0'],par['nv0']),dtype='complex64')
    evs=np.empty(par['nv0'],dtype='complex64')
    #np.info(mat)
    if solver==1:
        evs,evecs=lin.eig(mat)
    elif solver==2:
        evs=lin.eigvals(mat)
    else:
        stop

    #print evs
    gam=np.real(evs)
    ome=np.imag(evs)
    gam0=np.empty(par['nv0'])
    ome0=np.empty(par['nv0'])
    evecs0=np.empty((par['nv0'],par['nv0']),dtype='complex64')
    if show_plots:
      #plt.scatter(ome,gam)
      #plt.show()
      #matplotlib.rc('font',size=10)
      plt.figure(figsize=[4.0,3])
      #plt.title(r'$k_x \rho_i =$'+str(0)+r'$, k_y \rho_s =$'+str(ky)+r'$,k_z R = $'+str(kz))
      #####Plot eigenvalues
      plt.scatter(ome,gam,color='black',marker='x')
      #plt.scatter(svdo.omega,svdo.gamma,c='white',marker=marker_list,s=circle_size,cmap='jet')
      #####Plot pseudo-contours
      plt.xlabel(r'$\omega (v_{ti}/R)$',size=12)
      plt.ylabel(r'$\gamma (v_{ti}/R)$',size=12)
      plt.show()
      #plt.savefig(par['diagdir'][1:-1]+'/ev_spectrum.ps')
      #plt.close()
    for i in range(par['nv0']):
        mloc=np.argmax(gam)
        gam0[i]=gam[mloc]
        gam[mloc]=-1.0e15
        ome0[i]=ome[mloc]
        evecs0[:,i]=evecs[:,mloc]
    if plot_evecs:
        herm_grid=np.arange(par['nv0'])
        for i in range(num_evecs_plot):
            plt.plot(herm_grid,abs(evecs0[:,i]),label=str(float(i)))
        plt.title('Eigenvectors')
        plt.legend(loc='Upper Right')
        plt.xlabel('Hermite n')
        plt.show()
        for i in range(num_evecs_plot):
            plt.loglog(herm_grid,abs(evecs0[:,i]),label=str(float(i)),basex=10,basey=10)
        plt.title('Eigenvectors')
        plt.legend(loc='Upper Right')
        plt.xlabel('Hermite n')
        plt.show()
    return ome0,gam0,evecs0

def get_ev_spectra0(kx,ky,kz,test_extra=False,zamp=1.0,show_plots=False,plot_evecs=False,num_evecs_plot=8):

    mat=get_lin_matrix(kx,ky,kz,test_extra=test_extra,zamp=zamp)
    a,b,c=get_ev_spectrum(mat,show_plots=show_plots,plot_evecs=plot_evecs,num_evecs_plot=num_evecs_plot)
    return a,b,c


def my_corr_func_complex(v1,v2,time,show_plot=False,v1eqv2=True):
    #print "len(time)",len(time)
    #print "len(v1)",len(v1)
    #print "len(v2)",len(v2)
    timen,v1n=interpolate_grid_complex(time,v1) #,plot_data=True)
    timen,v2n=interpolate_grid_complex(time,v2)
    dt=timen[1]-timen[0]
    print "dt:", dt
    N=len(timen)
    cfunc=np.zeros(N,dtype='complex')
    for i in range(N):
        i0=i+1
        cfunc[-i0]=np.sum(np.conj(v1n[-i0:])*v2n[:i0])
    tau=np.arange(N)
    tau=tau*dt
    if v1eqv2:
        cfunc=np.real(cfunc)
    max_corr=max(np.abs(cfunc))
    corr_time=0.0
    i=0
    while corr_time==0.0:
        if (abs(cfunc[i])-max_corr/np.e) > 0.0 and \
           (abs(cfunc[i+1])-max_corr/np.e) <= 0.0:
            slope=(cfunc[i+1]-cfunc[i])/(tau[i+1]-tau[i])
            zero=cfunc[i]-slope*tau[i]
            corr_time=(max_corr/np.e-zero)/slope
        i+=1

    if show_plot:
        plt.plot(tau,cfunc,'x-')
        ax=plt.axis()
        plt.vlines(corr_time,ax[2],ax[3])
        plt.show()
    return cfunc,tau,corr_time

def interpolate_grid_complex(grid,data,plot_data=False):
    grid_new,rdata=interpolate_grid(grid,np.real(data),plot_data=False)
    grid_new,idata=interpolate_grid(grid,np.imag(data),plot_data=False)
    ndata=rdata+1.0J*idata
    if plot_data:
        plt.plot(grid,np.real(data),'-x',label='real old data')
        plt.plot(grid_new,np.real(ndata),'-+',label='real new data')
        plt.plot(grid,np.imag(data),'-x',label='imag old data')
        plt.plot(grid_new,np.imag(ndata),'-+',label='imag new data')
        plt.legend()
        plt.show()
    return grid_new,ndata

def interpolate_grid(grid,data,plot_data=False):
    minval=grid[0]
    maxval=grid[-1]
    dmin=min(np.abs(grid-np.roll(grid,1)))
    if dmin <= 0:
        print "Error in interpolate grid!!!"
        stop
    length=maxval-minval
    grid_new=np.arange(minval,maxval+dmin,dmin)
    ng=len(grid_new)
    tck=interpolate.splrep(grid,data)
    ndata=interpolate.splev(grid_new,tck,der=0)
    if plot_data:
        plt.plot(grid,data,'-x',label='old data')
        plt.plot(grid_new,ndata,'-+',label='new data')
        plt.show()

    return grid_new,ndata

def get_frequency_spectrum(time,signal,show_plots=False):
  nt=len(time)
  if nt%2 != 0:
      nt=nt-1
      time=time[1:]
      signal=signal[1:]
  Lt=time[-1]-time[0]
  #print "Length of time:",Lt
  #a=np.arange(n)/float(n-1)*Lt
  wmin=2.0*np.pi/Lt
  wmax=wmin*(nt/2-1)
  fa=np.arange(nt)/float(nt-1)*(2*wmax+wmin)-wmax

  #plt.plot(fa,'x')
  #plt.show()

  #print "wmin",wmin
  #print "wmax",wmax

  fs=np.fft.fft(signal)/float(nt)
  fs=np.roll(fs,nt/2-1)

  if show_plots:
      plt.plot(fa,abs(fs),'-x')
      plt.show()

  return fa,fs


def get_index_from_kx(kx):
    ik=np.rint(kx/par['kxmin'])
    return ik

def get_index_from_ky(ky):
    ik=np.rint(ky/par['kymin'])
    if ky < 0.0 and np.abs(ky) > 1.0e-15:
        ik=par['nky0']+ik
    return ik

def get_index_from_kz(kz):
    ik=np.rint(kz/par['kzmin'])
    if kz < 0.0 and np.abs(kz) > 1.0e-15:
        ik=par['nkz0']+ik
    return ik

def get_all_indices(kx,ky,kz):
    take_conjg=0
    ikx=get_index_from_kx(kx)
    if ikx<0:
        ikx=np.abs(ikx)
        take_conjg=1
    if take_conjg:
        iky=get_index_from_ky(-1.0*ky)
        ikz=get_index_from_kz(-1.0*kz)
    else:
        iky=get_index_from_ky(ky)
        ikz=get_index_from_kz(kz)
    return ikx,iky,ikz,take_conjg

#def get_checkpoint_from_gout():
#    time=get_time_from_gout()
#    numg=len(time)
#    gt0=read_time_step_g(numg)
#    gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
#    f = open('checkpoint','wb')
#    ntot=par['nkx0']*par['nky0']*par['nkz0']*par['nv0']
#    mem_tot=ntot*16
#    f.write(itime)
#    f.write(dt)
#    f.write(nkx0)
#    f.write(nky0)
#    f.write(nkz0)
#    f.write(nv0)
#    f.write(time)
#    f.close()




def black_color_map():
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (1.0, 0.0, 0.0)),
           'green': ((0.0, 0.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (1.0, 0.0, 0.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    #plt.pcolor(rand(10,10),cmap=my_cmap)
    #plt.colorbar()
    #plt.show()
    return my_cmap

def white_color_map():
    cdict = {'red': ((0.0, 1.0, 1.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
           'green': ((0.0, 1.0, 1.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
            'blue': ((0.0, 1.0, 1.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    #plt.pcolor(rand(10,10),cmap=my_cmap)
    #plt.colorbar()
    #plt.show()
    return my_cmap

def green_color_map():
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (1.0, 0.0, 0.0)),
           'green': ((0.0, 1.0, 1.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
            'blue': ((0.0, 0.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (1.0, 0.0, 0.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    #plt.pcolor(rand(10,10),cmap=my_cmap)
    #plt.colorbar()
    #plt.show()
    return my_cmap

def show_hyps():

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    #nu_kycomparison=np.empty((par['nky0']-2)/2+1)
    nu_kycomparison=np.empty(par['nky0'])
    #nu_kzcomparison=np.empty((par['nkz0']-2)/2+1)
    nu_kzcomparison=np.empty((par['nkz0']))
    nu_kxcomparison=np.empty(par['nkx0'])
    
    nu_kxcomparison[:]=par['nu']
    nu_kycomparison[:]=par['nu']
    nu_kzcomparison[:]=par['nu']
    
    coll=par['nu']*herm_grid
    #print coll
    hyp_coll=par['hyp_v']*(herm_grid/np.float(par['nv0']-1))**par['hypv_order']
    #temp=np.abs((hyp_coll-0.1*coll)/hyp_coll)
    temp=np.abs((hyp_coll-par['nu']))
    print "critical n:",np.argmin(temp)
    print (par['nu']/par['hyp_v'])**(1.0/float(par['hypv_order']))*par['nv0']


    
    plt.plot(herm_grid,coll,'x-',label='Collisions')
    plt.plot(herm_grid,hyp_coll,'+-',label='Hyper Collisions')
    plt.legend(loc='upper left')
    plt.xlabel('Hermite n')
    plt.title('Collisions and Hyper Collisions')
    plt.show()
    
    print par['hypx_order']
    
    #print "test"
    #print np.size(nu_kxcomparison)
    #print nu_kxcomparison
    #print np.size(kxgrid)
    #print "test"
    kxmax=par['kxmax0']
    kymax=par['kymax0']
    kzmax=par['kzmax0']
    hypx_out=par['hyp_x']*(kxgrid/kxmax)**par['hypx_order']
    plt.plot(kxgrid,hypx_out,label='Hyp_x')
    plt.plot(kxgrid,nu_kxcomparison,label='Minimum Coll')
    plt.plot(kxgrid,(par['nv0']-1)*nu_kxcomparison,label='Maximum Coll')
    plt.legend(loc='upper left')
    plt.xlabel(r'$k_x \rho_i$',size=18)
    plt.show()
    #hypx_out=par['hyp_x']*(kxgrid/kxmax)**par['hypx_order']
    #plt.semilogy(kxgrid,hypx_out,label='Hyp_x')
    #plt.legend(loc='upper left')
    #plt.xlabel(r'$k_x \rho_i$',size=18)
    #plt.show()
    
    
    hypy_out=par['hyp_y']*(kygrid/kymax)**par['hypy_order']
    plt.plot(kygrid,hypy_out,label='Hyp_y')
    plt.plot(kygrid,nu_kycomparison,label='Minimum Coll')
    plt.plot(kygrid,(par['nv0']-1)*nu_kycomparison,label='Maximum Coll')
    plt.legend(loc='upper left')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.show()
    #hypy_out=par['hyp_y']*(kygrid/kymax)**par['hypy_order']
    #plt.semilogy(kygrid,hypy_out,label='Hyp_y')
    #plt.legend(loc='upper left')
    #plt.xlabel(r'$k_y \rho_i$',size=18)
    #plt.show()
    
    hypz_out=par['hyp_z']*(kzgrid/kzmax)**par['hypz_order']
    plt.plot(kzgrid,hypz_out,label='Hyp_z')
    plt.plot(kzgrid,nu_kzcomparison,label='Minimum Coll')
    plt.plot(kzgrid,(par['nv0']-1)*nu_kzcomparison,label='Maximum Coll')
    plt.legend(loc='upper left')
    plt.xlabel(r'$k_z R$',size=18)
    plt.show()
    #

def paraview_test():
    vtk_writer=VTK_XML_Serial_Unstructured()
    x=np.arange(20)
    y=np.arange(20)
    z=np.arange(20)
    xvec=np.zeros(20)
    yvec=np.zeros(20)
    yvec[:]=1.0
    zvec=yvec*2.0
    vtk_writer.snapshot("vtk_test1.vtu",x,y,z,x_jump=xvec,y_jump=yvec,z_jump=zvec)
    vtk_writer.snapshot("vtk_test2.vtu",x,y,z,x_jump=xvec,y_jump=yvec,z_jump=zvec)
    vtk_writer.snapshot("vtk_test3.vtu",x,y,z,x_jump=xvec,y_jump=yvec,z_jump=zvec)
    vtk_writer.snapshot("vtk_test4.vtu",x,y,z,x_jump=xvec,y_jump=yvec,z_jump=zvec)
    vtk_writer.snapshot("vtk_test5.vtu",x,y,z,x_jump=xvec,y_jump=yvec,z_jump=zvec)
    vtk_writer.snapshot("vtk_test6.vtu",x,y,z,x_jump=xvec,y_jump=yvec,z_jump=zvec)
    vtk_writer.writePVD("pvd_test.pvd")

def plot_D_kperp():
    kx,ky,kz,herm=get_grids()
    J0p=np.arange(par['nky0']/2,dtype='float64')
    Gam0=np.arange(par['nky0']/2,dtype='float64')
    D0=np.arange(par['nky0']/2,dtype='float64')
    for j in range(par['nky0']/2):
        print "ky",ky[j],j
        Gam0[j]=sc.special.ive(0,ky[j]**2)
        print "Gam0",Gam0[j],j
        J0p[j]=np.e**(-ky[j]**2/2.0)
        print "J0",J0p[j],j
        D0[j]=1.0/(par['Ti0Te']+1-Gam0[j])
        print "D0",D0[j],j

    plt.plot(ky[0:par['nky0']/2],J0p,label='J0')
    plt.legend(loc='upper right')
    plt.show()
    plt.plot(ky[0:par['nky0']/2],Gam0,label='Gam0')
    plt.legend(loc='upper right')
    plt.show()
    plt.plot(ky[0:par['nky0']/2],D0,label='D0')
    plt.legend(loc='upper right')
    plt.show()
    plt.plot(ky[0:par['nky0']/2],J0p*D0,label='J0D0') 
    plt.legend(loc='upper right')
    plt.show()
    plt.plot(ky[0:par['nky0']/2],J0p,label='J0')
    plt.plot(ky[0:par['nky0']/2],Gam0,label='Gam0')
    plt.plot(ky[0:par['nky0']/2],D0,label='D0')
    plt.plot(ky[0:par['nky0']/2],J0p*D0,label='J0D0') 
    plt.legend(loc='upper right')
    plt.show()
                
def plot_Gam0():
    kx,ky,kz,herm=get_grids()
    Gam0=np.arange(par['nky0']/2,dtype='float64')
    term0=np.arange(par['nky0']/2,dtype='float64')
    term1=np.arange(par['nky0']/2,dtype='float64')
    term2=np.arange(par['nky0']/2,dtype='float64')
    term3=np.arange(par['nky0']/2,dtype='float64')
    for j in range(par['nky0']/2):
        print "ky",ky[j],j
        Gam0[j]=sc.special.ive(0,ky[j]**2)
        print "Gam0",Gam0[j],j
        term0[j]=1.0/sc.special.gamma(1.0)
        term1[j]=1.0/math.factorial(1.0)/sc.special.gamma(1.0+1.0)*(ky[j]**2/2.0)**2.0*1.0
        term2[j]=1.0/math.factorial(2.0)/sc.special.gamma(2.0+1.0)*(ky[j]**2/2.0)**2.0*2.0
        term3[j]=1.0/math.factorial(3.0)/sc.special.gamma(3.0+1.0)*(ky[j]**2/2.0)**2.0*3.0

    plt.plot(ky[0:par['nky0']/2],Gam0,label='Gam0')
    plt.plot(ky[0:par['nky0']/2],term0,label='1 term')
    plt.plot(ky[0:par['nky0']/2],term0+term1,label='2 terms')
    plt.plot(ky[0:par['nky0']/2],term0+term1+term2,label='3 terms')
    plt.plot(ky[0:par['nky0']/2],term0+term1+term2+term3,label='4 terms')
    plt.legend(loc='upper right')
    plt.show()
 
def plot_model_spect():
    kx,ky,kz,herm=get_grids()
    Gam0=np.arange(par['nky0']/2,dtype='float64')
    model=np.arange(par['nky0']/2,dtype='float64')
    spect7o3=np.arange(par['nky0']/2,dtype='float64')
    spect8o3=np.arange(par['nky0']/2,dtype='float64')
    spect4o3=np.arange(par['nky0']/2,dtype='float64')
    for j in range(par['nky0']/2):
        print "ky",ky[j],j
        Gam0[j]=sc.special.ive(0,ky[j]**2)
        print "Gam0",Gam0[j],j
        model[j]=(2-sc.special.ive(0,ky[j]**2))**(2.0/3.0)*ky[j]**(-4.0/3.0)*np.e**(2.0*ky[j]**2/3.0)
        spect7o3[j]=ky[j]**(-7.0/3.0)
        spect8o3[j]=ky[j]**(-8.0/3.0)
        spect4o3[j]=ky[j]**(-4.0/3.0)


    plt.loglog(ky[0:par['nky0']/2],model,label='model',basex=10,basey=10)
    plt.loglog(ky[0:par['nky0']/2],spect7o3,label='7o3',basex=10,basey=10)
    plt.loglog(ky[0:par['nky0']/2],spect8o3,label='8o3',basex=10,basey=10)
    plt.loglog(ky[0:par['nky0']/2],spect4o3,label='4o3',basex=10,basey=10)
    plt.legend(loc='upper right')
    plt.show()

def test_nuno_closure(g_in,kz):
    g_test=np.zeros(par['nv0'],dtype='complex')
    #kx,ky,kz,hg=get_grids()
    for i in range(par['nv0']-1):
        numerator=-1.0J*kz*(i+1)**0.5*g_in[i]
        denominator=par['nu']*(i+1)+par['hyp_v']*(float(i+1)/float(par['nv0']-1))**par['hypv_order']
        #print 'numerator',numerator
        #print 'denominator',denominator
        g_test[i+1]=numerator/denominator
                    #-1.0J*kz*(i+1)**0.5*g_in[i]\
                    #        /(par['nu']*(i+1)+par['hyp_v']*(float(i+1)/float(par['nv0']-1))**\
                    #par['hypv_order'])
        print np.abs(g_test[i+1])
    #plt.plot(np.real(g_in),label='real g_in')
    #plt.plot(np.real(g_test),label='real g_test')
    #plt.plot(np.imag(g_in),label='imag g_in')
    #plt.plot(np.imag(g_test),label='imag g_test')
    #plt.legend()
    #plt.show()
    hg=np.arange(par['nv0'])
    #plt.plot(hg,np.abs(g_test),'x',label='abs g_test')
    #plt.legend()
    #plt.show()
    plt.loglog(hg,np.abs(g_test),'x',label='abs g_test')
    plt.legend()
    plt.show()
    plt.loglog(np.abs(g_in),'+',label='abs g_in')
    plt.loglog(np.abs(g_test),'x',label='abs g_test')
    plt.legend()
    plt.show()
    #plt.loglog(np.abs(np.real(g_in)),'+',label='abs real g_in')
    #plt.loglog(np.abs(np.real(g_test)),'x',label='abs real g_test')
    #plt.legend()
    #plt.show()
 
def test_closure(kx,ky,kz):
    #kx,ky,kz,hg=get_grids()
    if_closure=par['nuno_closure']
    g_test=np.zeros(par['nv0']+1,dtype='complex')

    par['nuno_closure']=False
    mat=get_lin_matrix(kx,ky,kz)
    a1,b1,c1=get_ev_spectrum(mat)
    par['nuno_closure']=True
    mat=get_lin_matrix(kx,ky,kz)
    a2,b2,c2=get_ev_spectrum(mat)

    for i in range(par['nv0']):
        numerator=-1.0J*kz*(i+1)**0.5*c2[i,0]
        denominator=par['nu']*(i+1)+par['hyp_v']*(float(i+1)/float(par['nv0']-1))**par['hypv_order']
        g_test[i+1]=numerator/denominator
        print np.abs(g_test[i+1])
    hg=np.arange(par['nv0'])
    #plt.loglog(np.abs(g_test),'x',label='abs g_test')
    #plt.legend()
    #plt.show()
    plt.loglog(np.abs(c2[:,0]),'d-',label='abs ev w/c')
    plt.loglog(np.abs(c1[:,0]),'o-',markersize=4,label='abs ev wo/c')
    plt.loglog(np.abs(g_test),'x',label='abs g_test')
    plt.legend()
    plt.show()
    np.savetxt(par['diagdir'][1:-1]+'/eg_wclosure',np.abs(c2[:,0]))
    np.savetxt(par['diagdir'][1:-1]+'/eg_woclosure',np.abs(c1[:,0]))
    plt.scatter(a1,b1,marker='x',label='evs woc')
    plt.scatter(a2,b2,marker='+',label='evs wc')
    plt.legend()
    plt.show()


#############TEMP############
def test_gnlout():
    ikx = 0
    iky = 6
    ikz = 5
    in0 = 3
    short_factor = 10
    time = get_time_from_gout(gnl = True)
    print "len(time)",len(time)

    gt = np.empty(len(time)/short_factor,dtype='complex')
    for i in range(len(time)/short_factor):
        print i, " of ",len(time)/short_factor 
        #j = i+len(time)-len(time)/short_factor
        j = i
        gt0 = read_time_step_g(j)
        gt0 = np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        gt[i] = gt0[ikx,iky,ikz,in0]

    gnlt = np.empty(len(time)/short_factor,dtype='complex')
    for i in range(len(time)/short_factor):
        print i, " of ",len(time)/short_factor 
        #j = i+len(time)-len(time)/short_factor
        j = i
        gt0 = read_time_step_g(j,gnl=True)
        gt0 = np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        gnlt[i] = gt0[ikx,iky,ikz,in0]

    timek = get_time_from_gkfile(ikx,iky,read_nl=True)
    gnlk = np.empty(len(timek)/short_factor,dtype='complex')
    for i in range(len(timek)/short_factor):
        #j = i+len(timek)-len(timek)/short_factor
        j = i
        gtemp = read_time_step_gkfile(ikx,iky,j,read_nl=True)
        gtemp=np.reshape(gtemp,(par['nkz0'],par['nv0']),order='F')
        gnlk[i] = gtemp[ikz,in0]

    gk = np.empty(len(timek)/short_factor,dtype='complex')
    for i in range(len(timek)/short_factor):
        #j = i+len(timek)-len(timek)/short_factor
        j = i
        gtemp = read_time_step_gkfile(ikx,iky,j)
        gtemp = np.reshape(gtemp,(par['nkz0'],par['nv0']),order='F')
        gk[i] = gtemp[ikz,in0]

    print "len(gnlt)",len(gnlt)
    print "len(time[0:len(time)/short_factor])",len(timek[0:len(time)/short_factor])
    print "len(gnlk)",len(gnlk)
    print "len(timek[0:len(timek)/short_factor])",len(timek[0:len(timek)/short_factor])
    
    plt.plot(time[0:len(time)/short_factor],np.real(gnlt),'x-')
    plt.plot(timek[0:len(timek)/short_factor],np.real(gnlk))
    plt.show()

    plt.plot(time[0:len(time)/short_factor],np.real(gt),'x-')
    plt.plot(timek[0:len(timek)/short_factor],np.real(gk))
    plt.show()

    plt.plot(time[0:len(time)/short_factor],np.real(gt),'x-')
    plt.plot(time[0:len(time)/short_factor],np.real(gnlt),'x-')
    plt.show()


