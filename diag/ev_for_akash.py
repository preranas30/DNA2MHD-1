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


