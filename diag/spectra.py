import numpy as np
import matplotlib.pyplot as plt
from config import *
from dna_diags import *
import os
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

class svecs:
    '''class for svd vectors containing various properties'''
    def __init__(self,kx,ky,kz,vecs_in=0.0+0.0J,svals=0.0,tvecs=0.0,omega_max=2.0,gamma_max=2.0,\
                 nomega=200,ngamma=200,svd_omega=0.0,svd_gamma=0.0,calc_omegas=False,\
                 ev_error=0.0,read_from_file=False,diagdir='none'):
        #initialize:   stuff=mvecs(vec_array,kx,ky,kz,tvecs, . . . )
        if read_from_file:
            if diagdir=='none':
                diagdir=par['diagdir'][1:-1]
            file_in=diagdir+'/svecs_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
            print file_in
            f=open(file_in,'r')
            input=f.readline()
            print input
            in0=input.split(':')
            self.vec_size=int(float(in0[1])) 
            input=f.readline()
            print input
            in0=input.split(':')
            self.n_vecs=int(float(in0[1])) 
            input=f.readline()
            print input
            in0=input.split(':')
            self.kx=float(in0[1]) 
            input=f.readline()
            print input
            in0=input.split(':')
            self.ky=float(in0[1]) 
            input=f.readline()
            print input
            in0=input.split(':')
            self.kz=float(in0[1]) 
            input=f.readline()
            print input
            for i in range(12):
                input=f.readline()
                print input
            big_array=np.genfromtxt(f)
            self.svals=big_array[0:self.n_vecs]
            self.ev_error=big_array[self.n_vecs:2*self.n_vecs]
            self.E=big_array[2*self.n_vecs:3*self.n_vecs]
            self.Q=big_array[3*self.n_vecs:4*self.n_vecs]
            self.C=big_array[4*self.n_vecs:5*self.n_vecs]
            self.Ees=big_array[5*self.n_vecs:6*self.n_vecs]
            self.Een=big_array[6*self.n_vecs:7*self.n_vecs]
            self.Erhs=big_array[7*self.n_vecs:8*self.n_vecs]
            self.n_avg=big_array[8*self.n_vecs:9*self.n_vecs]
            self.omega=big_array[9*self.n_vecs:10*self.n_vecs]
            self.gamma=big_array[10*self.n_vecs:11*self.n_vecs]
            self.vectors=np.empty((self.vec_size,self.n_vecs),dtype='complex')
            basenum=11*self.n_vecs
            for i in range(self.n_vecs):
                self.vectors[:,i]=big_array[basenum+i*2*self.vec_size:basenum+(i*2+1)*self.vec_size]+(1.0J)*big_array[basenum+(i*2+1)*self.vec_size:basenum+(i*2+2)*self.vec_size]
        else:
            self.vectors=vecs_in
            self.vec_size=shape(vecs_in)[0]
            self.n_vecs=len(svals)
            self.kx=kx
            self.ky=ky
            self.kz=kz
            self.svals=svals
            self.ev_error=ev_error

            herm_grid=np.arange(par['nv0'])

            self.E=np.empty(0)
            self.Q=np.empty(0)
            self.C=np.empty(0)
            self.Ees=np.empty(0)
            self.Een=np.empty(0)
            self.Erhs=np.empty(0)
            self.n_avg=np.empty(0)
            for i in range(self.n_vecs):
                self.E=np.append(self.E,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,-1))
                self.Q=np.append(self.Q,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,6)/self.E[i])
                self.C=np.append(self.C,(get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,1)+\
                                        +get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,2)+\
                                        +get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,8))\
                                        /self.E[i])
                self.Ees=np.append(self.Ees,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,-2))
                self.Een=np.append(self.Een,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,-3))
                self.Erhs=np.append(self.Erhs,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,0)/self.E[i])
                self.n_avg=np.append(self.n_avg,np.sum(abs(vecs_in[:,i])*herm_grid)/np.sum(abs(vecs_in[:,i])))

            self.tsvd=tvecs
            if calc_omegas:
                self.omega=np.zeros(self.n_vecs)
                self.gamma=np.zeros(self.n_vecs)
                for i in range(self.n_vecs):
                    pev,a,b,c=get_pseudo_eigenvalue(vecs_in[:,i],kx,ky,kz,ngamma=ngamma,nomega=nomega,\
                                            omega_max=omega_max,gamma_max=gamma_max,plot_data=False)
                    self.omega[i]=imag(pev)
                    self.gamma[i]=real(pev)
                    self.ev_error[i]=c
            else:
                self.omega=svd_omega
                self.gamma=svd_gamma

    def write_svecs(self,diagdir='none',output_time=False):
        if diagdir=='none':
            diagdir=par['diagdir'][1:-1]
        file_out=diagdir+'/svecs_kx'+str(self.kx)+'ky'+str(self.ky)+'kz'+str(self.kz)
        print file_out
        f=open(file_out,'w')
        f.write('vec_size:')
        f.write(str(self.vec_size)+'\n')
        f.write('n_vecs:')
        f.write(str(self.n_vecs)+'\n')
        f.write('kx:')
        f.write(str(self.kx)+'\n')
        f.write('ky:')
        f.write(str(self.ky)+'\n')
        f.write('kz:')
        f.write(str(self.kz)+'\n')
        f.write('The following arrays:\n')
        f.write('svals:\n')
        f.write('ev_error:\n')
        f.write('E:\n')
        f.write('Q:\n')
        f.write('C:\n')
        f.write('Ees:\n')
        f.write('Een:\n')
        f.write('Erhs:\n')
        f.write('n_avg:\n')
        f.write('omega:\n')
        f.write('gamma:\n')
        f.write('real then imag vectors:\n')
        np.savetxt(f,self.svals)
        np.savetxt(f,self.ev_error)
        np.savetxt(f,self.E)
        np.savetxt(f,self.Q)
        np.savetxt(f,self.C)
        np.savetxt(f,self.Ees)
        np.savetxt(f,self.Een)
        np.savetxt(f,self.Erhs)
        np.savetxt(f,self.n_avg)
        np.savetxt(f,self.omega)
        np.savetxt(f,self.gamma)
        for i in range(self.n_vecs):
            np.savetxt(f,np.real(self.vectors[:,i]))
            np.savetxt(f,np.imag(self.vectors[:,i]))
        if output_time:
            stop
            f.write('real vectors:\n')
            np.savetxt(f,np.real(self.vectors))
            f.write('imag vectors:\n')
            np.savetxt(f,np.imag(self.vectors))
        f.close()

    def compare_svecs(self,svec2):
        print "vec_size diff:"
        print self.vec_size-svec2.vec_size
        print "n_vecs diff:"
        print self.n_vecs-svec2.n_vecs
        print "kx diff:"
        print self.kx-svec2.kx
        print "ky diff:"
        print self.ky-svec2.ky
        print "kz diff:"
        print self.kz-svec2.kz
        print "svals diff:"
        print np.sum(self.svals-svec2.svals)
        print "ev_error diff:"
        print np.sum(self.ev_error-svec2.ev_error)
        print "E diff:"
        print np.sum(self.E-svec2.E)
        print "Q diff:"
        print np.sum(self.Q-svec2.Q)
        print "C diff:"
        print np.sum(self.C-svec2.C)
        print "Ees diff:"
        print np.sum(self.Ees-svec2.Ees)
        print "Een diff:"
        print np.sum(self.Een-svec2.Een)
        print "Erhs diff:"
        print np.sum(self.Erhs-svec2.Erhs)
        print "n_avg diff:"
        print np.sum(self.n_avg-svec2.n_avg)
        print "omega diff:"
        print np.sum(self.omega-svec2.omega)
        print "gamma diff:"
        print np.sum(self.gamma-svec2.gamma)
        for i in range(self.n_vecs):
            print "vectors diff:"
            print np.sum(np.abs(self.vectors[:,i]-svec2.vectors[:,i]))

class evecs:
    'class for linear eigenvectors containing various properties'
    def __init__(self,vecs_in,kx,ky,kz,ev_omega=0.0,ev_gamma=0.0):
        #initialize:   stuff=mvecs(vec_array,kx,ky,kz,tvecs, . . . )
        #type
        #vectors
        #omega
        #gamma
        #tsvd (for SVD)
        self.vectors=vecs_in
        self.vec_size=shape(vecs_in)[0]
        self.n_vecs=shape(vecs_in)[1]
        self.kx=kx
        self.ky=ky
        self.kz=kz

        herm_grid=np.arange(par['nv0'])

        self.E=np.empty(0)
        self.Q=np.empty(0)
        self.C=np.empty(0)
        self.Ees=np.empty(0)
        self.Een=np.empty(0)
        self.Erhs=np.empty(0)
        self.n_avg=np.empty(0)
        for i in range(self.n_vecs):
            self.E=np.append(self.E,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,-1))
            self.Q=np.append(self.Q,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,6)/self.E[i])
            self.C=np.append(self.C,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,1)/self.E[i])
            self.Ees=np.append(self.Ees,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,-2))
            self.Een=np.append(self.Een,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,-3))
            self.Erhr=np.append(self.Erhs,get_energy_single_k(vecs_in[:,i],vecs_in[:,i],kx,ky,kz,0)/self.E[i])
            self.n_avg=np.append(self.n_avg,np.sum(abs(vecs_in[:,i])*herm_grid)/np.sum(abs(vecs_in[:,i])))

        self.omega=ev_omega
        self.gamma=ev_gamma


def get_pseudo_spectra(kx,ky,kz,nomega=40,ngamma=50,omega_max=0.8,\
                       show_plots=False,gamma_max=1.0,nlevels=50,verbose=False,\
                      read_ev_file=False,ev_file='eigenvalues'):
    mat=get_lin_matrix(kx,ky,kz)
    herm_grid=np.arange(par['nv0'])
    nimag=nomega
    nreal=ngamma
    deltaR=2*gamma_max/float(nreal-1)
    deltaI=2*omega_max/float(nimag-1)
    Id=np.zeros((par['nv0'],par['nv0']),dtype='complex64')
    for i in range(par['nv0']):
        Id[i,i]=1.0
    pspect=np.zeros((nreal,nimag),dtype='float32')
    pvec_prop=np.zeros((nreal,nimag),dtype='float32')
    pvec_prop2=np.zeros((nreal,nimag),dtype='float32')
    raxis=np.arange(nreal)/float(nreal-1)*2*gamma_max-gamma_max
    iaxis=np.arange(nimag)/float(nimag-1)*2*omega_max-omega_max
    for i in range(nreal):
        for j in range(nimag):
            z=raxis[i]+1.0J*iaxis[j]
            rmat=z*Id-mat
            ru,rs,rv=np.linalg.svd(rmat)
            pspect[i,j]=np.min(rs)
            pv=np.conj(rv[par['nv0']-1,:])
            pvec_prop[i,j]=(get_energy_single_k(pv,pv,kx,ky,kz,6))
            pvec_prop[i,j]=pvec_prop[i,j]/abs(get_energy_single_k(pv,pv,kx,ky,kz,-1))
            pvec_prop2[i,j]=np.real(np.sum(herm_grid*np.conj(pv)*pv)/np.sum(np.conj(pv)*pv))
            #pvec_prop[i,j]=pvec_prop[i,j]/get_energy_single_k(pv,pv,kx,ky,kz,-1)
            #print "Check pseudo eigenvector."
            #print np.dot(np.conj(pv),pv)
            #temp=np.dot(rmat,np.conj(pv))
            #print np.real(np.sqrt(np.dot(np.conj(temp),temp)))
            if verbose:
                print z,pspect[i,j]
    if read_ev_file:
        evff=np.genfromtxt(ev_file)
        ome=evff[:,0]
        gam=evff[:,1]
    else:
        ome,gam,evec=get_ev_spectrum(mat)

    #eveco=evecs(evec,kx,ky,kz,ome,gam)
    #def __init__(self,vecs_in,kx,ky,kz,ev_omega=0.0,ev_gamma=0.0):

    #plt.contour(iaxis,raxis,pspect,nlevels)
    #plt.xlabel(r'$\omega (v_{ti}/R)$',size=18)
    #plt.ylabel(r'$\gamma (v_{ti}/R)$',size=18)
    #plt.colorbar()
    #plt.show()

    if show_plots:
        plt.figure(figsize=[4,3])
        plt.title('Linear Spectrum')
        plt.scatter(ome,gam,color='black',marker='x')
        fig=plt.gcf()
        fig.subplots_adjust(bottom=0.20)
        fig.subplots_adjust(left=0.18)
        plt.xlabel(r'$\omega (v_{ti}/R)$',size=12)
        plt.ylabel(r'$\gamma (v_{ti}/R)$',size=12)
        axis_limits=plt.axis()
        if abs(axis_limits[2]) < 1.0e-13:
            ylim1=-0.1
        else:
            ylim1=axis_limits[2]
        if abs(axis_limits[3]) < 1.0e-13:
            ylim2=0.1
        else:
            ylim2=axis_limits[3]
        axis_limits=[axis_limits[0],axis_limits[1],ylim1,ylim2]
        plt.axis(axis_limits)
        plt.show()


        plt.figure(figsize=[4,3])
        plt.title('Pseudospectrum')
        plt.scatter(ome,gam,color='black',marker='x')
        plt.contour(iaxis,raxis,pspect,nlevels)
        fig=plt.gcf()
        fig.subplots_adjust(bottom=0.20)
        fig.subplots_adjust(left=0.18)
        plt.xlabel(r'$\omega (v_{ti}/R)$',size=12)
        plt.ylabel(r'$\gamma (v_{ti}/R)$',size=12)
        plt.colorbar()
        plt.show()

        plt.figure(figsize=[4,3])
        plt.title('Pseudo Eigenvectors Q')
        plt.scatter(ome,gam,color='black',marker='x')
        levels=[-0.02,10.0]#np.arange(100)/99.0*1.0
        #plt.contour(iaxis,raxis,pspect,nlevels)
        plt.contour(iaxis,raxis,pvec_prop,levels,linewidths=3,cmap=black_color_map())
        #plt.contour(iaxis,raxis,pvec_prop,nlevels)
        fig=plt.gcf()
        fig.subplots_adjust(bottom=0.20)
        fig.subplots_adjust(left=0.18)
        plt.xlabel(r'$\omega (v_{ti}/R)$',size=12)
        plt.ylabel(r'$\gamma (v_{ti}/R)$',size=12)
        plt.colorbar()
        plt.show()

        plt.figure(figsize=[4,3])
        plt.title('Pseudo Eigenvectors Q')
        plt.scatter(ome,gam,color='black',marker='x')
        plt.contour(iaxis,raxis,pvec_prop,nlevels)
        fig=plt.gcf()
        fig.subplots_adjust(bottom=0.20)
        fig.subplots_adjust(left=0.18)
        plt.xlabel(r'$\omega (v_{ti}/R)$',size=12)
        plt.ylabel(r'$\gamma (v_{ti}/R)$',size=12)
        plt.colorbar()
        levels=[0.0,10.0]#np.arange(100)/99.0*1.0
        plt.contour(iaxis,raxis,pvec_prop,levels,linewidths=3,cmap=black_color_map())
        plt.show()

        plt.figure(figsize=[4,3])
        levels=np.arange(20)/19.0*10.0
        plt.title('Pseudo Eigenvectors <n>')
        plt.scatter(ome,gam,color='black',marker='x')
        plt.contour(iaxis,raxis,pvec_prop2,levels)
        fig=plt.gcf()
        fig.subplots_adjust(bottom=0.20)
        fig.subplots_adjust(left=0.18)
        plt.xlabel(r'$\omega (v_{ti}/R)$',size=12)
        plt.ylabel(r'$\gamma (v_{ti}/R)$',size=12)
        plt.colorbar()
        levels=[-0.02,10.0]#np.arange(100)/99.0*1.0
        plt.contour(iaxis,raxis,pvec_prop,levels,linewidths=3,cmap=black_color_map())
        plt.show()
        
    return raxis,iaxis,pspect,pvec_prop

def get_svd_spectrum(kx,ky,kz,start_time=-1,end_time=-1,show_plots=False,\
                     use_gk=False,no_pe=100,ng_pe=100,\
                    get_ps_ev=True,verbose=False,\
                    read_all_g=True,gamma_max=2.0,omega_max=2.0,\
                    calc_evs=True,read_g=True,g_in=0.0,\
                    free_energy_weight=True):
    """Routine for getting the SVD spectrum and calculating various related quantities. \n
    kx,ky,kz: the wavenumbers to analyze \n
    start_time:  start time for analysis \n
    end_time:  end time for analysis \n
    show_plots: flag for showing plots \
    use_gk: flag for using gk_*** files instead of g_out.dat \n
    no_pe: number points in the real (frequency) coordinate of the complex plane to analyze for the \n
    pseudo-eigenvalue calculations \n
    ng_pe: number points in the imaginary (growth rate) coordinate of the complex plane to analyze for the \n
    pseudo-eigenvalue calculations \n
    get_ps_ev: flag for calculating pseudo eigenvalues for SVD vectors (otherwise logarithmic derivative) \n
    verbose: flag for outputting more comments \n
    read_all_g: read data from g_out.dat \n
    gamma_max:  maximum gamma for calculation of pseudo-eigenvalues \n
    omega_max:  maximum omega for calculation of pseudo-eigenvalues \n
    calc_evs:  flag for calculating pseudo eigenvalues via method determined by get_ps_ev \n
    read_g: flag for reading data from file.  If false, g must be passed in the call to get_svd_spectrum \n
    g_in:  input distribution function if read_g=False. \n
    free_energy_weight:  weights (and unweights) g so that the scalar product is the free energy. \n
    Returns: an \'svecs\' object constructed from the SVD."""


    if kx < 0.0:
        print "Error! kx must be positive!"
        stop
    ikx=get_kindex(kx,par['kxmin'],par['nkx0'])
    iky=get_kindex(ky,par['kymin'],par['nky0'])
    ikz=get_kindex(kz,par['kzmin'],par['nkz0'])
    #print 'ikx',ikx
    #print 'iky',iky
    #print 'ikz',ikz

    if use_gk:
        time=get_time_from_gk(ikx,iky,ikz)
    else:
        time=get_time_from_gout()

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    print "Reading data for svd."
    #g_k: Total distribution function
    g_k=np.empty((par['nv0'],ntime),dtype='complex')

    if read_g:
        for t in range(istart,iend+1):
            tind=t-istart
            if verbose:
                print 'time=',time[t],' of ',time[iend]
            if use_gk:
                g_k[:,tind]=read_time_step_gkfile(ikx,iky,ikz,t)
            elif read_all_g:
                gt0=read_time_step_g(t)
                gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
                g_k[:,tind]=gt0[ikx,iky,ikz,:]
            else:
                g_k[:,tind]=read_time_step_gk(kx,ky,kz,t)
    else:
        g_k=g_in

    if(free_energy_weight):
        print "Weighting distribution function."
        for t in range(ntime):
            g_k[:,t]=g_free_energy_weight(kx,ky,g_k[:,t])

    print "Calculating svd."
    print "shape(g_k)",np.shape(g_k)
    ru,rs,rv=np.linalg.svd(g_k)
    if show_plots:
        plt.plot(rs)
        plt.show()

        plt.semilogy(rs)
        plt.show()

    nsv=min(par['nv0'],ntime)

    if(free_energy_weight):
        print "Weighting distribution function."
        for i in range(nsv):
            ru[:,i]=g_free_energy_unweight(kx,ky,ru[:,i])

    if show_plots:
        for i in range(min(10,nsv)):
            plt.plot(time[istart:iend+1],abs(rv[i,:]))
            plt.show()

    herm_grid=np.arange(par['nv0'])
    #for i in range(min(10,nsv)):
    if show_plots:
        for i in range(min(10,nsv)):
            #plt.plot(herm_grid,abs(ru[i,:]))
            plt.plot(herm_grid,abs(ru[:,i]))
            plt.show()

    #######################Test#####################
    #######################Test#####################
    #######################Test#####################
    #mat=get_lin_matrix(kx,ky,kz)
    #om,gm,evecs=get_ev_spectrum(mat,show_plots=False)
    #ru=evecs
    #######################Test#####################
    #######################Test#####################
    #######################Test#####################

    #Calculate effective growth rates
    gam=np.empty(0)
    ome=np.empty(0)
    ev_error=np.empty(0)
    print "Calculating svd info."
    if calc_evs:
        for i in range(nsv):
            if get_ps_ev:
                pev,rax,iax,diff_om=get_pseudo_eigenvalue(ru[:,i],kx,ky,kz,\
                                                   plot_data=False,ngamma=ng_pe,nomega=no_pe,\
                                                     gamma_max=gamma_max,omega_max=omega_max)
                ome=np.append(ome,np.imag(pev))
                gam=np.append(gam,np.real(pev))
                ev_error=np.append(ev_error,diff_om)
            else:
                #Everything
                gam=np.append(gam,get_energy_single_k(ru[:,i],ru[:,i],kx,ky,kz,0))
                gam[i]=gam[i]/get_energy_single_k(ru[:,i],ru[:,i],kx,ky,kz,-1)
                gam[i]=gam[i]/2.0 #since energy is quadratic
                #Energy sources/sinks
                #gam=np.append(gam,get_energy_single_k(ru[:,i],ru[:,i],kx,ky,kz,1))
                #gam[i]=gam[i]+get_energy_single_k(ru[:,i],ru[:,i],kx,ky,kz,2)
                #gam[i]=gam[i]+get_energy_single_k(ru[:,i],ru[:,i],kx,ky,kz,6)
                #gam[i]=gam[i]+get_energy_single_k(ru[:,i],ru[:,i],kx,ky,kz,8)
                #gam[i]=gam[i]/get_energy_single_k(ru[:,i],ru[:,i],kx,ky,kz,-1)
                ome=np.append(ome,get_frequency_logd(time[istart:iend+1]\
                                                     ,rv[i,:],plot_data=False))
        svec=svecs(kx,ky,kz,vecs_in=ru,svals=rs,tvecs=rv,svd_omega=ome,svd_gamma=gam,ev_error=ev_error)
    else:
        svec=svecs(kx,ky,kz,vecs_in=ru,svals=rs,tvecs=rv)


    if show_plots:
        plt.scatter(ome,gam)
        plt.show()

    return svec


    #####Frequency analysis tests####
    #####Frequency analysis tests####
    #####Frequency analysis tests####
    #tnew,rvnew=interpolate_grid_complex(time,rv[2,:],plot_data=True)

    #plt.semilogy(tnew,abs(rvnew),label='abs time trace')
    #plt.show()

    #om_grid,fspect=get_frequency_spectrum(tnew,rvnew)
    #lorentz_fit_par,cov=curve_fit(lorentzian,om_grid,abs(fspect)**2)
    #lor_fit=lorentzian(om_grid,lorentz_fit_par[0],lorentz_fit_par[1],\
    #                   lorentz_fit_par[2])
    #print lorentz_fit_par
    #plt.plot(om_grid,abs(fspect)**2)
    #plt.plot(om_grid,lor_fit)
    #plt.show()
    #omega_logd,gamma_logd=get_frequency_logd(tnew,rvnew)
    #plt.plot(tnew,omega_logd,label='omega')
    #plt.plot(tnew,gamma_logd,label='gamma')
    #plt.legend()
    #plt.show()

    #####Frequency analysis tests####
    #####Frequency analysis tests####
    #####Frequency analysis tests####

def g_free_energy_weight(kx,ky,g_in):
    phi_pre=get_phi_prefactor(kx,ky)
    #print 'phi_pre',phi_pre
    g_out=g_in*(0.5*np.pi**0.5)**0.5
    g_out[0]=g_out[0]*(1+np.pi**(-0.25)*phi_pre*np.e**(-0.5*(kx**2+ky**2)))**0.5
    ####Test####
    #print "Free energy:",get_energy_single_k(g_in,g_in,kx,ky,0.2,-1)
    #print "This calc:",np.sum(np.conj(g_out)*g_out)
    return g_out

def g_free_energy_unweight(kx,ky,g_in):
    phi_pre=get_phi_prefactor(kx,ky)
    #print 'phi_pre',phi_pre
    g_out=g_in/(0.5*np.pi**0.5)**0.5
    g_out[0]=g_out[0]/(1+np.pi**(-0.25)*phi_pre*np.e**(-0.5*(kx**2+ky**2)))**0.5
    ####Test####
    #print "Free energy:",get_energy_single_k(g_in,g_in,kx,ky,0.2,-1)
    #print "This calc:",np.sum(np.conj(g_out)*g_out)
    return g_out

def scalar_product_free_energy(g_1,g_2,kx,ky):
    """Returns the scalar product defined by the free energy"""
    vec1=g_free_energy_weight(kx,ky,g_1)
    vec2=g_free_energy_weight(kx,ky,g_2)
    return real(np.sum(np.conj(vec1)*vec2))

def plot_all_spectra(kx,ky,kz,nomega=40,ngamma=50,omega_max=0.8,gamma_max=1.5,\
                     nlevels=50,start_time=-1,end_time=-1,use_gk=False,\
                     get_ps_ev=True,\
                     no_pe=100,ng_pe=100,show_plots=False,verbose=False):
    print "kx,ky,kz",kx,ky,kz
    print "Getting eigenvalue spectrum."
    mat=get_lin_matrix(kx,ky,kz)
    omega,gamma,evec=get_ev_spectrum(mat)
    eveco=evecs(evec,kx,ky,kz,omega,gamma)
    print "Getting pseudo spectrum."
    raxis,iaxis,pspect,pvec_prop = get_pseudo_spectra(kx,ky,kz,nomega=nomega,ngamma=ngamma,omega_max=omega_max,gamma_max=gamma_max,nlevels=nlevels)
    #svdomega,svdgam,svec,tvec,svals=get_svd_spectrum(kx,ky,kz,start_time=start_time,end_time=end_time,\
    print "Getting svd spectrum."
    svdo=get_svd_spectrum(kx,ky,kz,start_time=start_time,end_time=end_time,\
                                          show_plots=False,use_gk=use_gk,\
                                          get_ps_ev=get_ps_ev,\
                                          ng_pe=ng_pe,no_pe=no_pe,read_all_g=True,\
                         gamma_max=gamma_max,omega_max=omega_max)

    mode_grid=np.arange(svdo.n_vecs)+1
    zero_array=np.zeros(svdo.n_vecs+1)
    zgrid=np.arange(svdo.n_vecs+1)

    multiple_symbols=True
    if multiple_symbols:
        marker_list=list('x' for i in range(len(svdo.omega)))
        marker_list[0]='o'
        marker_list[1]='^'
        marker_list[2]='s'
        marker_list[3]='d'
        marker_list[4]='p'
        marker_list[5]='h'
        marker_list[6]='<'
        marker_list[7]='>'
        marker_list[8]='v'
        marker_list[9]='+'

  
    if show_plots:
        gs=gridspec.GridSpec(4,1,height_ratios=[3,1,1,1])
        plt.subplots_adjust(hspace=0.4)
        plt.subplot(gs[0])
        plt.title(r'$k_x \rho_i =$'+str(0)+r'$, k_y \rho_s =$'+str(ky)+r'$,k_z R = $'+str(kz))
        plt.scatter(omega,gamma,color='black',marker='x')
        circle_size=svdo.ev_error*300.0
        #print circle_size
        plt.scatter(svdo.omega,svdo.gamma,c='white',s=circle_size,cmap='jet')
        svtemp=svdo.svals/svdo.svals[0]
        plt.scatter(svdo.omega,svdo.gamma,c=np.log10(svtemp),s=30,edgecolors='none',cmap='jet',\
                   vmin=(np.log10(svtemp[0])-3))
        plt.colorbar()
        plt.contour(iaxis,raxis,pspect,nlevels)
        levels=[-0.02,10.0]
        plt.xlabel(r'$\omega (v_{ti}/R)$',size=14)
        plt.ylabel(r'$\gamma (v_{ti}/R)$',size=14)
        plt.colorbar()
        #plt.contour(iaxis,raxis,pvec_prop,levels,cmap=black_color_map(),linewidths=3)
        plt.axis([-omega_max,omega_max,-1.0,abs(2.0*np.max(svdo.gamma))])
        plt.subplot(gs[1])
        plt.plot(mode_grid,svdo.Q,'x',markeredgewidth=2,label='Q/E')
        plt.plot(mode_grid,svdo.C,'+',markeredgewidth=2,label='C/E')
        plt.plot(zgrid,zero_array,color='black')
        plt.legend(loc='upper right',numpoints=1)
        axis_limits=plt.axis()
        axis_limits=[0.5,11.0,1.1*svdo.C[9],axis_limits[3]]
        plt.axis(axis_limits)
        plt.subplot(gs[2])
        temp=svdo.svals**2*(svdo.Q[0:svdo.n_vecs]+svdo.C[0:svdo.n_vecs])
        plt.axis([0.5,11,1.1*np.min(temp),-0.1*np.min(temp)])
        plt.plot(mode_grid,temp,'x',marker=marker_list,markeredgewidth=2)
        plt.plot(zgrid,zero_array,color='black')
        plt.ylabel(r'$\sigma_n^2(Q_n+C_n)$',size=14)
        #plt.plot(svdo.svals**2*(svdo.Q[0:svdo.n_vecs]),'+')
        #plt.plot(svdo.svals**2*(svdo.C[0:svdo.n_vecs]),'*')
        plt.subplot(gs[3])
        plt.plot(mode_grid,svdo.ev_error**2,'x',markeredgewidth=2)
        plt.axis([0.5,11,0.0,1.0])
        plt.ylabel(r'$\delta^2$',size=14)
        plt.xlabel('Mode number n')
        #plt.plot(svdo.C,'+')
        plt.show()
    else:
        matplotlib.rc('font',size=10)
        gs=gridspec.GridSpec(4,1,height_ratios=[5,1,1,1])
        plt.figure(figsize=[4.0,7])
        plt.subplots_adjust(hspace=0.2,left=0.17)
        plt.subplot(gs[0])
        plt.title(r'$k_x \rho_i =$'+str(0)+r'$, k_y \rho_s =$'+str(ky)+r'$,k_z R = $'+str(kz))
        #####Plot eigenvalues
        plt.scatter(omega,gamma,color='black',marker='x')
        #plt.scatter(svdo.omega,svdo.gamma,c='white',marker=marker_list,s=circle_size,cmap='jet')
        #####Plot pseudo-contours
        plt.contour(iaxis,raxis,pspect,nlevels)
        levels=[-0.02,10.0]
        plt.xlabel(r'$\omega (v_{ti}/R)$',size=12)
        plt.ylabel(r'$\gamma (v_{ti}/R)$',size=12)
        plt.colorbar(orientation='horizontal',pad=0.0,ticks=[0.20,0.4,0.6,0.8])
        #plt.contour(iaxis,raxis,pvec_prop,levels,cmap=black_color_map(),linewidths=3)

        #####Plot pseudo-eigenvalues
        circle_size=80
        if multiple_symbols:
            for i in range(len(svdo.omega)):
                plt.scatter(svdo.omega[i],svdo.gamma[i],c='white',marker=marker_list[i],s=circle_size,cmap='jet')
        else:
            plt.scatter(svdo.omega,svdo.gamma,c='white',marker=marker_list,s=circle_size,cmap='jet')

        svtemp=svdo.svals/svdo.svals[0]
        plt.scatter(svdo.omega,svdo.gamma,c=np.log10(svtemp),s=30,edgecolors='none',cmap='jet',\
                   vmin=(np.log10(svtemp[0])-3))
        plt.colorbar(orientation='horizontal',pad=0.2,ticks=[-2.5,-2.0,-1.5,-1.0,-0.5])
        ########
        plt.axis([-omega_max,omega_max,min(1.1*min(svdo.gamma[0:4]),-1.0),abs(2.0*np.max(svdo.gamma))])
        plt.annotate('A.',xy=(0.05,0.05),xycoords='axes fraction')
        plt.subplot(gs[1])
        if multiple_symbols:
            plt.plot(mode_grid[0],svdo.Q[0],color='blue',marker=marker_list[0],label='Q/E')
            plt.plot(mode_grid[0],svdo.C[0],color='green',marker=marker_list[0],label='C/E')
            for i in range(1,min(len(svdo.omega),10)):
                plt.plot(mode_grid[i],svdo.Q[i],color='blue',marker=marker_list[i])
                plt.plot(mode_grid[i],svdo.C[i],color='green',marker=marker_list[i])
        else:
            plt.plot(mode_grid,svdo.Q,'x',markeredgewidth=2,label='Q/E')
            plt.plot(mode_grid,svdo.C,'+',markeredgewidth=2,label='C/E')
        plt.plot(zgrid,zero_array,color='black')
        plt.legend(loc='upper right',numpoints=1,frameon=True)
        axis_limits=plt.axis()
        axis_limits=[0.5,11.0,1.1*svdo.C[9],axis_limits[3]]
        plt.axis(axis_limits)
        plt.annotate('B.',xy=(0.05,0.05),xycoords='axes fraction')
        plt.subplot(gs[2])
        temp=svdo.svals**2*(svdo.Q[0:svdo.n_vecs]+svdo.C[0:svdo.n_vecs])
        plt.axis([0.5,11,1.1*np.min(temp),-0.1*np.min(temp)])
        if multiple_symbols:
            plt.plot(mode_grid[0],temp[0],marker=marker_list[0],color='blue',label=r'$\sigma_n^2(Q_n+C_n)$')
            for i in range(1,min(len(svdo.omega),10)):
                plt.plot(mode_grid[i],temp[i],color='blue',marker=marker_list[i])
        else:
            plt.plot(mode_grid,temp,'x',markeredgewidth=2,label=r'$\sigma_n^2(Q_n+C_n)$')
        #plt.plot(mode_grid,temp,'x',markeredgewidth=2,label=r'$\sigma_n^2(Q_n+C_n)$')
        plt.plot(zgrid,zero_array,color='black')
        plt.legend(loc='lower right',numpoints=1,frameon=True)
        #plt.ylabel(r'$\sigma_n^2(Q_n+C_n)$',size=10)
        #plt.plot(svdo.svals**2*(svdo.Q[0:svdo.n_vecs]),'+')
        #plt.plot(svdo.svals**2*(svdo.C[0:svdo.n_vecs]),'*')
        plt.annotate('C.',xy=(0.05,0.05),xycoords='axes fraction')
        plt.subplot(gs[3])
        if multiple_symbols:
            plt.plot(mode_grid[0],svdo.ev_error[0]**2,marker=marker_list[0],color='blue')
            for i in range(1,min(len(svdo.omega),10)):
                plt.plot(mode_grid[i],svdo.ev_error[i]**2,marker=marker_list[i],color='blue')
        else:
            plt.plot(mode_grid,svdo.ev_error**2,'x',markeredgewidth=2)
        #plt.plot(mode_grid,svdo.ev_error**2,'x',markeredgewidth=2)
        plt.axis([0.5,11,0.0,1.0])
        plt.ylabel(r'$\delta^2$',size=12)
        plt.xlabel('Mode number n')
        plt.annotate('D.',xy=(0.05,0.75),xycoords='axes fraction')
        #plt.plot(svdo.C,'+')
        plt.savefig(par['diagdir'][1:-1]+'/spectra_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)+'.ps')
        plt.close()
 
    return svdo,eveco

def spectra_plot_scan():

    kx=0.0
    for i in range(16):
        for j in range(10):
            ky=(1+i)*par['kymin']
            kz=(1+j)*par['kzmin']
            omega_max=6.0*kz
            plot_all_spectra(kx,ky,kz,omega_max=omega_max,nomega=100,ngamma=100,start_time=150,\
                             no_pe=200,ng_pe=200,show_plots=False)

    #par_svec=np.empty(len(svec[:,0]))
    ##for i in range(len(svec[:,0])):
    #    #par_svec[i]=abs(get_energy_single_k(svec[:,i],svec[:,i],kx,ky,kz,-2)/\
    #    #            get_energy_single_k(svec[:,i],svec[:,i],kx,ky,kz,-3))
    #    #print par_svec[i]
    #for i in range(len(svec[:,0])):
    #    par_svec[i]=get_energy_single_k(svec[:,i],svec[:,i],kx,ky,kz,6)/\
    #                abs(get_energy_single_k(svec[:,i],svec[:,i],kx,ky,kz,1))

    #par_evec=np.empty(len(evec[:,0]))
    ##for i in range(len(evec[:,0])):
    ##    par_evec[i]=abs(get_energy_single_k(evec[:,i],evec[:,i],kx,ky,kz,-2)/\
    ##                get_energy_single_k(evec[:,i],evec[:,i],kx,ky,kz,-3))
    ##    print par_evec[i]
    #for i in range(len(evec[:,0])):
    #    par_evec[i]=get_energy_single_k(evec[:,i],evec[:,i],kx,ky,kz,6)/\
    #                abs(get_energy_single_k(evec[:,i],evec[:,i],kx,ky,kz,1))

    #plt.plot(par_svec,'x',label='SVEC')
    #plt.plot(par_evec,'+',label='EVEC')
    #plt.legend()
    #plt.show()
    #plt.semilogy(par_svec,'x',label='SVEC')
    #plt.semilogy(par_evec,'+',label='EVEC')
    #plt.legend()
    #plt.show()
    
def test_diag_gk(ikx,iky,ikz):
    #Note: test case for istep_gk=istep_gout
    #Successful for test case
    time=get_time_from_gout()
    timek=get_time_from_gk(ikx,iky,ikz)
    ntg=len(time[1:])
    ntgk=len(timek)
    ntime=max(ntg,ntgk)
    print "ntg,ntgk",ntg,ntgk
    if ntg != ntgk:
        stop
    for i in range(ntime):
        gk1=read_time_step_gkfile(ikx,iky,ikz,i)
        gt0=read_time_step_g(i+1)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        gk2=gt0[ikx,iky,ikz,:]
        diff=np.sum(abs(gk1-gk2))/np.sum(abs(gk2))
        print diff

def get_time_from_gk(ikx,iky,ikz):
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
    
    cikz=str(int(ikz))
    if len(cikz)==1:
        cikz='0'+cikz
    elif len(cikz) > 2:
        print "Error in get_time_from_gk"
        stop

    diagdir=par['diagdir'][1:-1]
    print diagdir
        
    file_name='g_kx'+cikx+'ky'+ciky+'kz'+cikz+'.dat'
    print "Reading file",diagdir+'/'+file_name
    file_exists=os.path.isfile(diagdir+'/'+file_name)
    if file_exists:
        pass
    else:
        print "File does not exist:",diagdir+'/'+file_name
        stop

    f=open(diagdir+'/'+file_name,'r')
    ntot=par['nv0']
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

def read_time_step_gkfile(ikx,iky,ikz,which_itime):
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
    
    cikz=str(int(ikz))
    if len(cikz)==1:
        cikz='0'+cikz
    elif len(cikz) > 2:
        print "Error in get_time_from_gk"
        stop

    diagdir=par['diagdir'][1:-1]
    #print diagdir
        
    file_name='g_kx'+cikx+'ky'+ciky+'kz'+cikz+'.dat'
    #print "Reading file",diagdir+'/'+file_name
    file_exists=os.path.isfile(diagdir+'/'+file_name)
    if file_exists:
        pass
    else:
        print "File does not exist:",diagdir+'/'+file_name
        stop


    f = open(diagdir+'/'+file_name,'rb')
    ntot=par['nv0']
    mem_tot=ntot*16
    gt0=np.empty((par['nv0']))
    f.seek(8+which_itime*(8+mem_tot))
    gt0=np.fromfile(f,dtype='complex128',count=ntot)
    #print sum(gt0)
    f.close()
    return gt0

def lorentzian(grid,h,gamma,x0):
    a=gamma*np.pi*h/2.0
    lor=0.5*gamma*a/np.pi/((grid-x0)**2+(0.5*gamma)**2)
    return lor

def get_frequency_logd(time,signal,plot_data=False):
    #Note: time and signal must be already interpolated onto even grid
    #print "np.info(signal)",np.info(signal)
    log_sig=np.log(signal)
    #print "np.info(log_sig)",np.info(log_sig)
    num_time=len(time)
    #plt.plot(time,'x')
    plt.show()
    dt=time[1]-time[0]
    dlog_sig=np.empty(num_time,dtype='complex64')
    dlog_sig[0]=0.0
    for i in range(num_time-1):
        dlog_sig[i+1]=(log_sig[i+1]-log_sig[i])/dt
        if abs(dlog_sig[i+1]-dlog_sig[i]) > 10.0:
            dlog_sig[i+1]=dlog_sig[i]
    omega=np.imag(dlog_sig)
    gamma=np.real(dlog_sig)
    omega_avg=np.sum(omega)/real(num_time)
    print "avg omega",omega_avg
    if plot_data:
        plt.plot(time,omega,label='omega')
        plt.plot(time,gamma,label='gamma')
        plt.legend()
        plt.show()

    return omega_avg

def get_pseudo_eigenvalue(g_in,kx,ky,kz,ngamma=400,nomega=400,omega_max=2.0,\
                      gamma_max=2.0,plot_data=False,verbose=False,energy_sp=True):
    nimag=nomega
    nreal=ngamma
    rhs=get_rhs_lin_single_k(g_in,0,kx,ky,kz)
    deltaI=2*omega_max/float(nimag-1)
    deltaR=2*gamma_max/float(nreal-1)
    diff_omega=np.zeros((nimag,nreal),dtype='float32')
    raxis=np.arange(nreal)/float(nreal-1)*2*gamma_max-gamma_max
    iaxis=np.arange(nimag)/float(nimag-1)*2*omega_max-omega_max
    minval=1.0e10
    for i in range(nimag):
        for j in range(nreal):
            z=raxis[j]+1.0J*iaxis[i]
            vec_temp=rhs-z*g_in
            if energy_sp:
                norm=np.sqrt(scalar_product_free_energy(rhs,rhs,kx,ky))+\
                     np.sqrt(scalar_product_free_energy(z*g_in,z*g_in,kx,ky))
                diff_omega[i,j]=np.abs(np.sqrt(scalar_product_free_energy(vec_temp,vec_temp,kx,ky))/norm)
            else:
                norm=np.sqrt(np.dot(np.conj(rhs),rhs))+np.sqrt(np.dot(np.conj(z*g_in),z*g_in))
                diff_omega[i,j]=np.abs(np.sqrt(np.dot(np.conj(vec_temp),vec_temp))/norm)
            if diff_omega[i,j] < minval:
                pev=z
                minval=diff_omega[i,j]
            #print z,diff_omega[i,j]

    ev_error=np.min(diff_omega)
    if verbose:
        print pev,ev_error
    if plot_data:
        plt.contour(raxis,iaxis,np.transpose(diff_omega),150)
        plt.show()
    return pev,raxis,iaxis,ev_error
    
def plot_gamma_kykz(kymax,kzmax):
    kx=0.0
    nky=int(kymax/par['kymin'])
    nkz=int(kzmax/par['kzmin'])
    kygrid=(np.arange(nky)+1.0)*par['kymin']
    kzgrid=(np.arange(nkz)+1.0)*par['kzmin']
    gmax=np.zeros((nky,nkz))
    zero_boundary=np.empty(0)
    maxloc_kz=np.empty(nky)
    for i in range(nky):
        for j in range(nkz):
            ky=(i+1)*par['kymin']
            kz=(j+1)*par['kzmin']
            mat=get_lin_matrix(kx,ky,kz)
            ome,gam,evecs=get_ev_spectrum(mat,show_plots=False)
            gmax[i,j]=np.max(gam)
            print ky,kz,gmax[i,j]
        maxloc_kz[i]=(np.argmax(gmax[i,:])+1)*par['kzmin']
    plt.scatter(maxloc_kz,kygrid)
    plt.contour(kzgrid,kygrid,gmax,50)
    plt.xlabel(r'$k_z R$',size=18)
    plt.ylabel(r'$k_y \rho_i$',size=18)
    plt.colorbar()
    plt.show()

def get_pseudo_eigenvector(kx,ky,kz,z_in):

    mat=get_lin_matrix(kx,ky,kz)
    Id=np.zeros((par['nv0'],par['nv0']),dtype='complex64')
    for i in range(par['nv0']):
        Id[i,i]=1.0
    #print z
    ru,rs,rv=np.linalg.svd((z_in*Id-mat))
    #print rs
    pspect=np.min(rs)
    print "pspect",pspect
    #print rs
    #print ru
    #print rv
    print shape(rv)
    plt.plot(abs(rv[par['nv0']-1,:]))
    plt.show()
    return pspect,rv[par['nv0']-1,:]

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

def svd_scan(allk=False,start_time=-1,end_time=-1,data_out=True,\
            num_kx=-1,num_ky=-1,num_kz=-1):

    if allk:
        print "allk option not ready at this time."
    else:

        kx=0.0
        if num_ky==-1:
            num_ky=par['nky0']
        if num_kz==-1:
            num_kz=par['nkz0']

        time=get_time_from_gout()

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
        print num_ky*num_kz*par['nv0']*ntime*16/1000000.0, "MB"

        gkykz=np.empty((num_ky,num_kz,par['nv0'],ntime),dtype='complex')
        for t in range(istart,iend+1):
            tind=t-istart
            print 'time=',time[t],' of ',time[iend]
            gt0=read_time_step_g(t)
            gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
            gkykz[:,:,:,tind]=gt0[0,0:num_ky,0:num_kz,:]

        Qsource=np.zeros((num_ky,num_kz)) #All source modes
        Csource=np.zeros((num_ky,num_kz))  
        Esource=np.zeros((num_ky,num_kz))  
        Gsource=np.zeros((num_ky,num_kz))  
        QsinkDW=np.zeros((num_ky,num_kz)) #All drift wave-like sink modes
        CsinkDW=np.zeros((num_ky,num_kz))
        EsinkDW=np.zeros((num_ky,num_kz))
        GsinkDW=np.zeros((num_ky,num_kz))
        QsinkL =np.zeros((num_ky,num_kz))  #All Landau-like sink modes
        CsinkL =np.zeros((num_ky,num_kz))
        EsinkL =np.zeros((num_ky,num_kz))
        GsinkL =np.zeros((num_ky,num_kz))

        for j in range(num_ky):
            for k in range(num_kz):
                ky=j*par['kymin']
                kz=k*par['kzmin']
                print "kx,ky,kz",0.0,ky,kz
                svdo = get_svd_spectrum(0.0,ky,kz,start_time=start_time,end_time=end_time,\
                                        calc_evs=False,read_g=False,g_in=gkykz[j,k,:,:])
                for i in range(svdo.n_vecs):
                    if svdo.Q[i]+svdo.C[i] > 0.0:
                        Qsource[j,k]+=svdo.svals[i]**2*svdo.Q[i]
                        Csource[j,k]+=svdo.svals[i]**2*svdo.C[i]
                        Esource[j,k]+=svdo.svals[i]**2*svdo.E[i]
                    else:
                        if abs(svdo.Q[i]/svdo.C[i]) > 0.05:
                            QsinkDW[j,k]+=svdo.svals[i]**2*svdo.Q[i]
                            CsinkDW[j,k]+=svdo.svals[i]**2*svdo.C[i]
                            EsinkDW[j,k]+=svdo.svals[i]**2*svdo.E[i]
                        else:
                            QsinkL[j,k]+=svdo.svals[i]**2*svdo.Q[i]
                            CsinkL[j,k]+=svdo.svals[i]**2*svdo.C[i]
                            EsinkL[j,k]+=svdo.svals[i]**2*svdo.E[i]
                if Esource[j,k] > 1.0e-12:
                    Gsource[j,k]=(Qsource[j,k]+Csource[j,k])/Esource[j,k]
                if EsinkDW[j,k] > 1.0e-12:
                    GsinkDW[j,k]=(QsinkDW[j,k]+CsinkDW[j,k])/EsinkDW[j,k]
                if EsinkL[j,k] > 1.0e-12:
                    GsinkL[j,k]=(QsinkL[j,k]+CsinkL[j,k])/EsinkL[j,k]


    if data_out:
        np.savetxt(par['diagdir'][1:-1]+'/Qsource.dat',Qsource)
        np.savetxt(par['diagdir'][1:-1]+'/Csource.dat',Csource)
        np.savetxt(par['diagdir'][1:-1]+'/Esource.dat',Esource)
        np.savetxt(par['diagdir'][1:-1]+'/QsinkDW.dat',QsinkDW)
        np.savetxt(par['diagdir'][1:-1]+'/CsinkDW.dat',CsinkDW)
        np.savetxt(par['diagdir'][1:-1]+'/EsinkDW.dat',EsinkDW)
        np.savetxt(par['diagdir'][1:-1]+'/QsinkL.dat',QsinkL)
        np.savetxt(par['diagdir'][1:-1]+'/CsinkL.dat',CsinkL)
        np.savetxt(par['diagdir'][1:-1]+'/EsinkL.dat',EsinkL)

    kxgrid,kygrid,kzgrid,hgrid=get_grids()
    print 'kxgrid',kxgrid
    print 'kygrid',kygrid
    print 'kzgrid',kzgrid
    for i in range(num_kz):
        plt.title(r'$k_z R=$'+str(i*par['kzmin']))
        plt.plot(kygrid[0:num_ky],Qsource[:,i],label='Q source')
        plt.plot(kygrid[0:num_ky],Csource[:,i],label='C source')
        plt.xlabel(r'$k_y \rho_i$',size=18)
        plt.legend()
        plt.show()
        plt.title(r'$k_z R=$'+str(i*par['kzmin']))
        plt.plot(kygrid[0:num_ky],QsinkDW[:,i],label='Q sink DW')
        plt.plot(kygrid[0:num_ky],CsinkDW[:,i],label='C sink DW')
        plt.plot(kygrid[0:num_ky],QsinkL[:,i],label='Q sink L')
        plt.plot(kygrid[0:num_ky],CsinkL[:,i],label='C sink L')
        plt.xlabel(r'$k_y \rho_i$',size=18)
        plt.legend()
        plt.show()

        plt.plot(kygrid[0:num_ky],QsinkDW[:,i]+QsinkL[:,i],label='Q sink')
        plt.plot(kygrid[0:num_ky],CsinkDW[:,i]+CsinkL[:,i],label='C sink')
        plt.xlabel(r'$k_y \rho_i$',size=18)
        plt.legend()
        plt.show()

        plt.title(r'$k_z R=$'+str(i*par['kzmin']))
        plt.plot(kygrid[0:num_ky],Gsource[:,i],label='gam source')
        plt.plot(kygrid[0:num_ky],GsinkDW[:,i],label='gam DW')
        plt.plot(kygrid[0:num_ky],GsinkL[:,i],label='gam L')
        plt.xlabel(r'$k_y \rho_i$',size=18)
        plt.legend()
        plt.show()

    return Qsource,Csource,Esource,Gsource,QsinkDW,CsinkDW,EsinkDW,GsinkDW,QsinkL,CsinkL,EsinkL,GsinkL

def svd_scan0(start_time=-1,end_time=-1,data_out=True,\
            verbose=True,xbins=1,plot_combo_spectra=False,\
             energy_balance=False,\
             xbin_start=0,ky_start_index=1,kz_start_index=0,\
             xbin_end=-1,ky_end_index=-1,kz_end_index=-1,\
             file_suffix='.dat'):

    time=get_time_from_gout()

    if ky_end_index==-1:
        ky_end_index=par['nky0']
    if kz_end_index==-1:
        kz_end_index=par['nkz0']
    if xbin_end==-1:
        xbin_end=xbins

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    if(verbose):
        print "Number of time points:",ntime

    if energy_balance:
        Qsource=np.zeros((par['nkx0'],par['nky0'],par['nkz0'])) #All source modes
        Csource=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))  
        Esource=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))  
        Gsource=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))  
        Qsink=np.zeros((par['nkx0'],par['nky0'],par['nkz0'])) #All drift wave-like sink modes
        Csink=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
        Esink=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
        Gsink=np.zeros((par['nkx0'],par['nky0'],par['nkz0']))
    else:
        Qsource=0
        Csource=0
        Esource=0
        Gsource=0
        Qsink=0
        Csink=0
        Esink=0
        Gsink=0

    print "ntime",ntime
    print par['nkx0']/xbins*par['nkz0']*par['nv0']*ntime*16/1000000.0, "MB"

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    nkx_bin=par['nkx0']/xbins

    gkxkz=np.empty((par['nkx0']/xbins,par['nkz0'],par['nv0'],ntime),dtype='complex')
    for i0 in range(xbin_start,xbin_end):
      for j in range(ky_start_index,ky_end_index):
      #for j in range(1,2):
        for t in range(istart,iend+1):
            tind=t-istart
            print 'time=',time[t],' of ',time[iend]
            gt0=read_time_step_g(t)
            gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
            gkxkz[:,:,:,tind]=gt0[i0*nkx_bin:(i0+1)*nkx_bin,j,:,:]

        if data_out and energy_balance:
            #Basically checkpoints
            print "Writing data."
            np.save(par['diagdir'][1:-1]+'/Qsource'+file_suffix,Qsource)
            np.save(par['diagdir'][1:-1]+'/Csource'+file_suffix,Csource)
            np.save(par['diagdir'][1:-1]+'/Esource'+file_suffix,Esource)
            np.save(par['diagdir'][1:-1]+'/Gsource'+file_suffix,Gsource)
            np.save(par['diagdir'][1:-1]+'/Qsink'+file_suffix,Qsink)
            np.save(par['diagdir'][1:-1]+'/Csink'+file_suffix,Csink)
            np.save(par['diagdir'][1:-1]+'/Esink'+file_suffix,Esink)
            np.save(par['diagdir'][1:-1]+'/Gsink'+file_suffix,Gsink)
            f=open(par['diagdir'][1:-1]+'/source_sink_info.dat','w')
            f.write("xbin,j: "+str(i0)+','+str(j))
            f.close()
            print "Done."

        ky=kygrid[j]
        for i in range(i0*nkx_bin,(i0+1)*nkx_bin):
        #for i in range(i0*nkx_bin,(i0)*nkx_bin+2):
            for k in range(kz_start_index,kz_end_index):
            #for k in range(4):
                kx=kxgrid[i]
                kz=kzgrid[k]
                print "kx,ky,kz",kx,ky,kz
                if energy_balance:
                    svdo = get_svd_spectrum(kx,ky,kz,start_time=start_time,end_time=end_time,\
                                        calc_evs=False,read_g=False,g_in=gkxkz[i-i0*nkx_bin,k,:,:])
                    for n in range(svdo.n_vecs):
                        if svdo.Q[n]+svdo.C[n] > 0.0:
                            Qsource[i,j,k]+=svdo.svals[n]**2*svdo.Q[n]
                            Csource[i,j,k]+=svdo.svals[n]**2*svdo.C[n]
                            Esource[i,j,k]+=svdo.svals[n]**2*svdo.E[n]
                        else:
                            Qsink[i,j,k]+=svdo.svals[n]**2*svdo.Q[n]
                            Csink[i,j,k]+=svdo.svals[n]**2*svdo.C[n]
                            Esink[i,j,k]+=svdo.svals[n]**2*svdo.E[n]
                    if Esource[i,j,k] > 1.0e-12:
                        Gsource[i,j,k]=(Qsource[i,j,k]+Csource[i,j,k])/Esource[i,j,k]
                    if Esink[i,j,k] > 1.0e-12:
                        Gsink[i,j,k]=(Qsink[i,j,k]+Csink[i,j,k])/Esink[i,j,k]
                if ((ky <= 1.4) and (ky > 0.0) and (kz <= 2.5) and (kz >=1.2) and (kx == 0) and plot_combo_spectra):
                        omega_max=min(1,6.0*kz)
                        plot_all_spectra(kx,ky,kz,omega_max=omega_max,nomega=100,ngamma=100,start_time=150,\
                                 no_pe=150,ng_pe=150,show_plots=False)


    if data_out and energy_balance:
        np.save(par['diagdir'][1:-1]+'/Qsource.dat',Qsource)
        np.save(par['diagdir'][1:-1]+'/Csource.dat',Csource)
        np.save(par['diagdir'][1:-1]+'/Esource.dat',Esource)
        np.save(par['diagdir'][1:-1]+'/Gsource.dat',Gsource)
        np.save(par['diagdir'][1:-1]+'/Qsink.dat',Qsink)
        np.save(par['diagdir'][1:-1]+'/Csink.dat',Csink)
        np.save(par['diagdir'][1:-1]+'/Esink.dat',Esink)
        np.save(par['diagdir'][1:-1]+'/Gsink.dat',Gsink)

    kxgrid,kygrid,kzgrid,hgrid=get_grids()
    print 'kxgrid',kxgrid
    print 'kygrid',kygrid
    print 'kzgrid',kzgrid

    #for i in range(num_kz):
    #    plt.title(r'$k_z R=$'+str(i*par['kzmin']))
    #    plt.plot(kygrid[0:num_ky],Qsource[:,i],label='Q source')
    #    plt.plot(kygrid[0:num_ky],Csource[:,i],label='C source')
    #    plt.xlabel(r'$k_y \rho_i$',size=18)
    #    plt.legend()
    #    plt.show()
    #    plt.title(r'$k_z R=$'+str(i*par['kzmin']))
    #    plt.plot(kygrid[0:num_ky],QsinkDW[:,i],label='Q sink DW')
    #    plt.plot(kygrid[0:num_ky],CsinkDW[:,i],label='C sink DW')
    #    plt.plot(kygrid[0:num_ky],QsinkL[:,i],label='Q sink L')
    #    plt.plot(kygrid[0:num_ky],CsinkL[:,i],label='C sink L')
    #    plt.xlabel(r'$k_y \rho_i$',size=18)
    #    plt.legend()
    #    plt.show()

    #    plt.plot(kygrid[0:num_ky],QsinkDW[:,i]+QsinkL[:,i],label='Q sink')
    #    plt.plot(kygrid[0:num_ky],CsinkDW[:,i]+CsinkL[:,i],label='C sink')
    #    plt.xlabel(r'$k_y \rho_i$',size=18)
    #    plt.legend()
    #    plt.show()

    #    plt.title(r'$k_z R=$'+str(i*par['kzmin']))
    #    plt.plot(kygrid[0:num_ky],Gsource[:,i],label='gam source')
    #    plt.plot(kygrid[0:num_ky],GsinkDW[:,i],label='gam DW')
    #    plt.plot(kygrid[0:num_ky],GsinkL[:,i],label='gam L')
    #    plt.xlabel(r'$k_y \rho_i$',size=18)
    #    plt.legend()
    #    plt.show()

    return Qsource,Csource,Esource,Gsource,Qsink,Csink,Esink,Gsink

def ev_scan(kx_in,kymax=2.0,kzmax=1.5,show_kz=False):
    kxgrid,kygrid,kzgrid,hermgrid=get_grids()
    j=0
    k=0
    nky=int(kymax/par['kymin'])
    nkz=int(kzmax/par['kzmin'])
    gam=np.zeros((nky,nkz))
    print kygrid
    for j in range(nky):
        for k in range(nkz):
            print kygrid[j],kzgrid[k]
            mat=get_lin_matrix(kx_in,kygrid[j],kzgrid[k])
            om,gm,ev=get_ev_spectrum(mat)
            gam[j,k]=gm[0]

    if show_kz:
        for k in range(nky):
            plt.plot(kzgrid[0:nkz],gam[k,:],'x-')
            plt.xlabel('kz')
            plt.title('ky='+str(kygrid[k]))
            plt.show()
    else:
        for k in range(nkz):
            plt.plot(kygrid[0:nky],gam[:,k],'x-')
            plt.xlabel('ky')
            plt.title('kz='+str(kzgrid[k]))
            plt.show()

def ev_ky_plot(kx_in,kymax=2.0,kzmax=1.5,ikz_start=0):
    kxgrid,kygrid,kzgrid,hermgrid=get_grids()
    j=0
    k=0
    nky=int(kymax/par['kymin'])
    nkz=int(kzmax/par['kzmin'])
    gam=np.zeros((nky,nkz))
    ome=np.zeros((nky,nkz))
    gam_max_ky=np.zeros(nky)
    ome_max_ky=np.zeros(nky)
    print kygrid
    for j in range(nky):
        for k in range(ikz_start,nkz):
            print kygrid[j],kzgrid[k]
            mat=get_lin_matrix(kx_in,kygrid[j],kzgrid[k])
            om,gm,ev=get_ev_spectrum(mat)
            gam[j,k]=gm[0]
            ome[j,k]=om[0]

    for j in range(nky):
        gam_max_ky[j]=max(gam[j,:])
        loc=np.argmax(gam[j,:])
        ome_max_ky[j]=ome[j,loc]

    plt.plot(kygrid[0:nky],gam_max_ky,'x-',label=r'$\gamma R/v_{ti} $')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.ylabel(r'$\gamma R/v_{ti}$',size=18)
    plt.show()
    plt.plot(kygrid[0:nky],ome_max_ky,'x-',label=r'$\omega R/v_{ti} $')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.ylabel(r'$\omega R/v_{ti}$',size=18)
    plt.show()

    plt.title(r'$\gamma R/v_{ti}$')
    plt.contour(kzgrid[0:nkz],kygrid[0:nky],gam,100)
    plt.xlabel(r'$k_z R$',size=18)
    plt.ylabel(r'$k_y \rho_i$',size=18)
    plt.show()

    plt.title(r'$\omega R/v_{ti}$')
    plt.contour(kzgrid[0:nkz],kygrid[0:nky],ome,100)
    plt.xlabel(r'$k_z R$',size=18)
    plt.ylabel(r'$k_y \rho_i$',size=18)
    plt.show()

def ev_ky_scan(ikx,ikz,kymax=1.0,evnum=0):
    kxgrid,kygrid,kzgrid,hermgrid=get_grids()
    kx=kxgrid[ikx]
    kz=kzgrid[ikz]
    nky=int(kymax/par['kymin'])
    gam=np.zeros((nky))
    ome=np.zeros((nky))
    for j in range(nky):
            print kygrid[j],kz
            mat=get_lin_matrix(kx,kygrid[j],kz)
            om,gm,ev=get_ev_spectrum(mat)#,show_plots=True)
            gam[j]=gm[evnum]
            ome[j]=om[evnum]

    plt.plot(kygrid[0:nky],gam,'x-',label=r'$\gamma L/v_{ti} $')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.ylabel(r'$\gamma L/v_{ti}$',size=18)
    plt.show()
    plt.plot(kygrid[0:nky],ome,'x-',label=r'$\omega L/v_{ti} $')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.ylabel(r'$\omega L/v_{ti}$',size=18)
    plt.show()

def ev_kx_scan(iky,ikz,kxmax=1.0,evnum=0):
    kxgrid,kygrid,kzgrid,hermgrid=get_grids()
    ky=kygrid[iky]
    kz=kzgrid[ikz]
    nkx=int(kxmax/par['kxmin'])
    gam=np.zeros((nkx))
    ome=np.zeros((nkx))
    for j in range(nkx):
            print kxgrid[j]
            mat=get_lin_matrix(kxgrid[j],ky,kz)
            om,gm,ev=get_ev_spectrum(mat)
            gam[j]=gm[evnum]
            ome[j]=om[evnum]

    plt.plot(kxgrid[0:nkx],gam,'x-',label=r'$\gamma L/v_{ti} $')
    plt.xlabel(r'$k_x \rho_i$',size=18)
    plt.ylabel(r'$\gamma L/v_{ti}$',size=18)
    plt.show()
    plt.plot(kxgrid[0:nkx],ome,'x-',label=r'$\omega L/v_{ti} $')
    plt.xlabel(r'$k_x \rho_i$',size=18)
    plt.ylabel(r'$\omega L/v_{ti}$',size=18)
    plt.show()

def ev_kz_scan(ikx,iky,kzmax=1.0,evnum=0):
    kxgrid,kygrid,kzgrid,hermgrid=get_grids()
    kx=kxgrid[ikx]
    ky=kygrid[iky]
    nkz=int(kzmax/par['kzmin'])
    gam=np.zeros((nkz))
    ome=np.zeros((nkz))
    for j in range(nkz):
            print kzgrid[j]
            mat=get_lin_matrix(kx,ky,kzgrid[j])
            om,gm,ev=get_ev_spectrum(mat)
            gam[j]=gm[evnum]
            ome[j]=om[evnum]

    plt.plot(kzgrid[0:nkz],gam,'x-',label=r'$\gamma L/v_{ti} $')
    plt.xlabel(r'$k_z L$',size=18)
    plt.ylabel(r'$\gamma L/v_{ti}$',size=18)
    plt.show()
    plt.plot(kygrid[0:nkz],ome,'x-',label=r'$\omega L/v_{ti} $')
    plt.xlabel(r'$k_z L$',size=18)
    plt.ylabel(r'$\omega L/v_{ti}$',size=18)
    plt.show()

def get_energy_spectra_from_svdscan(file_suffix='.dat'):
    """For plotting energy output from svd_scan0"""

    Qsource=np.load(par['diagdir'][1:-1]+'/Qsource'+file_suffix+'.npy')
    Csource=np.load(par['diagdir'][1:-1]+'/Csource'+file_suffix+'.npy')
    Esource=np.load(par['diagdir'][1:-1]+'/Esource'+file_suffix+'.npy')
    Gsource=np.load(par['diagdir'][1:-1]+'/Gsource'+file_suffix+'.npy')
    Qsink=np.load(par['diagdir'][1:-1]+'/Qsink'+file_suffix+'.npy')
    Csink=np.load(par['diagdir'][1:-1]+'/Csink'+file_suffix+'.npy')
    Esink=np.load(par['diagdir'][1:-1]+'/Esink'+file_suffix+'.npy')
    Gsink=np.load(par['diagdir'][1:-1]+'/Gsink'+file_suffix+'.npy')
    
    plot_energy_quantity(Qsource+Csource,"Q+C Source Modes.")
    plot_energy_quantity(Qsink+Csink,"Q+C Sink Modes.")
    plot_energy_quantity(Qsource+Csource,"Energy Source Modes.")
    plot_energy_quantity(Qsink+Csink,"Energy Sink Modes.")
    plot_energy_quantity(Qsource+Csource,"Gamma Source Modes.")
    plot_energy_quantity(Qsink+Csink,"Gamma Sink Modes.")

def plot_energy_quantity(en_in3d,energy_title,plotlog=False):
    """Plots contour and k spectra for a given energy quantity \n
    en_in3d: 3d input array \n
    energy_title: string defining energy quantity."""


    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    kxgrid_full=np.empty(par['nkx0']*2-1)
    for i in range(par['nkx0']):
        kxgrid_full[i+par['nkx0']-1]=kxgrid[i]
    for i in range(par['nkx0']-1):
        kxgrid_full[i]=-1.0*kxgrid[par['nkx0']-1-i]
    print kxgrid_full

    kx_out,ky_out,kz_out,herm_out=get_grids_shifted()
    e2d_plot=np.empty((par['nkx0']*2-1,par['nky0']))

    en_kza=np.sum(en_in3d,2)
    np.info(en_kza)
    e2d_plot=plot_prep2d(en_kza)

    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    fig.subplots_adjust(left=0.2)
    fig.subplots_adjust(bottom=0.2)
    plt.contourf(kxgrid_full,ky_out,np.transpose(e2d_plot),200)
    plt.title(energy_title)
    plt.xlabel(r'$k_x\rho_i$',size=13)
    plt.ylabel(r'$k_y\rho_i$',size=13)
    plt.colorbar()
    plt.show()

    en_kx=np.sum(en_kza,1)
    en_kx[1:]=2.0*en_kx[1:]
    en_ky=np.sum(en_kza[:,:],0)
    for i in range(1,par['nky0']/2):
        en_ky[i]=en_ky[i]+en_ky[par['nky0']-i]
        en_ky[par['nky0']-i]=en_ky[i]
    en_ky=np.roll(en_ky,par['nky0']/2-1,axis=0)
    en_kya=np.sum(en_in3d,1)
    en_kz=np.sum(en_kya[:,:],0)
    for i in range(1,par['nkz0']/2):
        en_kz[i]=en_kz[i]+en_kz[par['nkz0']-i]
        en_kz[par['nkz0']-i]=en_kz[i]
    en_kz=np.roll(en_kz,par['nkz0']/2-1,axis=0)

    #Linear plot
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    fig.subplots_adjust(left=0.2)
    fig.subplots_adjust(bottom=0.2)
    plt.plot(kx_out,en_kx)
    plt.xlabel(r'$k_x\rho_i$',size=13)
    plt.title(energy_title)
    plt.xlim((0,par['kxmax0']))
    plt.show()



    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    fig.subplots_adjust(left=0.2)
    fig.subplots_adjust(bottom=0.2)
    plt.plot(ky_out,en_ky)
    plt.title(energy_title)
    plt.xlabel(r'$k_y\rho_i$',size=13)
    plt.xlim((0,par['kymax0']))
    plt.show()
    plt.figure(figsize=(4.5,3.0))
    fig=plt.gcf()
    fig.subplots_adjust(left=0.2)
    fig.subplots_adjust(bottom=0.2)
    plt.plot(kz_out,en_kz)
    plt.title(energy_title)
    plt.xlabel(r'$k_z R$',size=13)
    plt.xlim((0,par['kzmax0']))
    plt.show()

    #log-log plot
    if plotlog:
        plt.figure(figsize=(4.5,3.0))
        fig=plt.gcf()
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(kx_out,-1.0*en_kx,basex=10,basey=10)
        plt.title(energy_title)
        plt.xlabel(r'$k_x\rho_i$',size=13)
        plt.show()
        plt.figure(figsize=(4.5,3.0))
        fig=plt.gcf()
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(ky_out,-1.0*en_ky,basex=10,basey=10)
        plt.title(energy_title)
        plt.xlabel(r'$k_y\rho_i$',size=13)
        plt.show()
        plt.figure(figsize=(4.5,3.0))
        fig=plt.gcf()
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.22)
        plt.loglog(kz_out,-1.0*en_kz,basex=10,basey=10)
        plt.title(energy_title)
        plt.xlabel(r'$k_z R$',size=13)
        plt.show()

def plot_vec_QC_E(veco):
    zgrid=np.arange(64)
    zero_array=np.zeros(64)
    mode_grid=np.arange(veco.n_vecs)+1
    plt.plot(mode_grid,veco.Q,'x',markeredgewidth=2,label='Q/E')
    plt.plot(mode_grid,veco.C,'+',markeredgewidth=2,label='C/E')
    plt.plot(zgrid,zero_array,color='black')
    plt.legend(loc='upper right',numpoints=1)
    axis_limits=plt.axis()
    axis_limits=[0.5,11.0,1.1*veco.C[9],axis_limits[3]]
    plt.axis(axis_limits)
    plt.ylabel(r'$(Q_n,C_n)E_n$',size=14)
    plt.xlabel('Mode number n')
    plt.show()

    herm_grid=np.arange(par['nv0'])
    for i in range(4,9):
        plt.plot(herm_grid,abs(veco.vectors[:,i]),label=str(float(i)))
    plt.title('Eigenvectors')
    plt.legend(loc='Upper Right')
    plt.xlabel('Hermite n')
    plt.show()

