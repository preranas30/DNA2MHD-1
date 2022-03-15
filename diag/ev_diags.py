#!/usr/bin/env python
# File: ev_diags.py

import numpy as np
import matplotlib.pyplot as plt
from config import *
from dna_diags import *

def ev_projection(time,kx,ky,kz,start_time=-1.0,end_time=-1.0,num_keep=6):

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    new_dir=par['diagdir'][1:len(par['diagdir'])-1]+'/ev_files/'
    ev_file='ev_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    print ev_file
    evals=np.genfromtxt(new_dir+ev_file)
    print evals

    ###Read in left and right eigenvectors
    revec_file='revec_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    rf=open(new_dir+revec_file,'rb')
    #evf=open(new_dir+ev_file,'rb')
    revecs=np.empty((par['nv0'],par['nv0']),dtype='complex')
    ev_num=np.array(range(par['nv0']))
    for i in range(par['nv0']):
        rf.seek(i*4+i*par['nv0']*16)
        ev_num[i]=np.fromfile(rf,dtype='int32',count=1)
        rf.seek((i+1)*4+i*par['nv0']*16)
        revecs[:,i]=np.fromfile(rf,dtype='complex128',count=par['nv0'])
        revecs[:,i]=revecs[:,i]/\
                np.sqrt(np.sum(np.conj(revecs[:,i])*revecs[:,i]))

    ortho_tol=1.0e-10
    #num_keep=6

    ru,rs,rv=np.linalg.svd(revecs[:,0:num_keep])
    #print "svals of revec",rs
    revecs_new=np.empty((par['nv0'],par['nv0']),dtype='complex')
    revecs_new[:,0:num_keep]=revecs[:,0:num_keep]
    revecs_new[:,num_keep:]=ru[:,num_keep:]
    ru,rs,rv=np.linalg.svd(revecs_new)
    print "svals of revec_new",rs
    levecs_new=np.linalg.pinv(revecs_new,rcond=10.0**-8)
    levecs_new=np.conj(levecs_new)

    products=np.empty((par['nv0'],par['nv0']),dtype='float')
    prod_quality=np.empty(par['nv0'])
    cutoff_test=0
    cutoff=0
    for i in range(par['nv0']):
        ###Normalize left eigenvectors
        #prod=abs(sum(np.conj(levecs[:,i])*revecs[:,i]))
        prod=(np.sum(np.conj(levecs_new[i,:])*revecs_new[:,i]))
        levecs_new[i,:]=levecs_new[i,:]/prod
        for j in range(par['nv0']):
            ###Note: Must use one of the following
            products[i,j]=np.real(np.sum(np.conj(levecs_new[i,:])*revecs_new[:,j]))
            #products[i,j]=abs(np.sum(np.conj(levecs[:,i])*revecs[:,j]))
        
        prod_quality[i]=np.sum(products[i,:])-products[i,i]
        if prod_quality[i] < ortho_tol and cutoff_test==0:
            cutoff=i
        else:
            cutoff_test=1
        
    ngood=cutoff+1

    print "Number of good Left EV's",ngood
    nv0=par['nv0']
    #xaxis=np.arange((par['nv0']))+1
    #yaxis=np.arange(cutoff+4)+1
    xaxis=np.arange(par['nv0'])+1
    yaxis=np.arange(par['nv0'])+1
    levels=np.arange(0.0,1.01,0.01)
    plt.contour(xaxis,yaxis,products,levels)
    #plt.title('Number of good (err < 10^-3) eigenvectors: '+str(cutoff+1))
    plt.title('Right vs. Left Eigenvectors from pseudo inverse')
    plt.xlabel('Right Eigenvector')
    plt.ylabel('Left Eigenvector')
    plt.colorbar()
    plt.show()
    #Get indices corresponding to k's
    if kx < 0.0:
        print "Error! kx must be positive!"
        stop
    ikx=get_kindex(kx,par['kxmin'],par['nkx0'])
    iky=get_kindex(ky,par['kymin'],par['nky0'])
    ikz=get_kindex(kz,par['kzmin'],par['nkz0'])
    print 'ikx',ikx
    print 'iky',iky
    print 'ikz',ikz

    #g_k: Total distribution function
    g_k=np.empty((par['nv0'],ntime),dtype='complex')
    #g_kr: Total distribution function reproduced from decomposition
    g_kr=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_knk: Sum of all ev's used
    g_knk=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_landau: Sum of all Landau ev's used
    g_landau=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_res: Residual distribution function: g_k-g_knk
    g_res=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_1: dist. func. for ITG mode
    g_1=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_2: dist. func. for Mode 2 - i.e. drift wave
    g_2=np.zeros((par['nv0'],ntime),dtype='complex')
    ar_n=np.empty((nv0,ntime),dtype='complex')
    for t in range(istart,iend+1):
        tind=t-istart
        print 'time=',time[t],' of ',time[iend]
        gt0=read_time_step_g(t)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        g_k[:,tind]=gt0[ikx,iky,ikz,:]
        ###Find coefficients in REV basis
        for evn in range(nv0):
            ar_n[evn,tind]=np.sum(np.conj(levecs_new[evn,:])*g_k[:,tind])
            g_kr[:,tind]=g_kr[:,tind]+ar_n[evn,tind]*revecs_new[:,evn]

        for evn in range(num_keep):        
            #Get sum of the modes we are keeping
            g_knk[:,tind]=g_knk[:,tind]+ar_n[evn,tind]*revecs_new[:,evn]

        g_1[:,tind]=ar_n[0,tind]*revecs_new[:,0]
        g_2[:,tind]=ar_n[1,tind]*revecs_new[:,1]
        #Get sum of the Landau modes we are keeping
        g_landau[:,tind]=g_knk[:,tind]-g_1[:,tind]-g_2[:,tind]
        #Residual distribution function
        g_res[:,tind]=g_k[:,tind]-g_knk[:,tind]
        print np.abs(np.sum(np.conj(g_res[:,tind])*g_res[:,tind]))/np.abs(np.sum(np.conj(g_k[:,tind])*g_k[:,tind]))

    g_kr=g_1+g_2+g_landau+g_res
    diff=g_k-g_kr
    #plt.plot(time[istart:iend+1],np.sum(np.abs(diff),axis=0),color='black',label='Error')
    #plt.show()
    print np.sum(np.abs(diff),axis=0)
    diff=g_k-g_knk-g_res
    print np.sum(np.abs(diff),axis=0)

    #Get Entropy stuff
    ent_diag=np.empty((num_keep,ntime))
    ent_tot=np.empty(ntime)
    ent_12=np.empty(ntime)
    ent_1L=np.empty(ntime)
    ent_2L=np.empty(ntime)
    ent_Ldiag=np.empty(ntime)
    ent_res=np.empty(ntime)
    ent_keep_res=np.empty(ntime)
    for t in range(istart,iend+1):
        tind=t-istart
        for evn in range(num_keep):
            ent_diag[evn,tind]=0.5*np.pi**0.5*np.real(\
                               np.sum(np.conj(ar_n[evn,tind]*revecs_new[:,evn])*\
                                        ar_n[evn,tind]*revecs_new[:,evn]))
        ent_12[tind]=0.5*np.pi**0.5*np.real(\
                               np.sum(np.conj(g_1[:,tind])*g_2[:,tind])+\
                               np.sum(np.conj(g_2[:,tind])*g_1[:,tind]))
        ent_1L[tind]=0.5*np.pi**0.5*np.real(\
                               np.sum(np.conj(g_1[:,tind])*g_landau[:,tind])+\
                               np.sum(np.conj(g_landau[:,tind])*g_1[:,tind]))
        ent_2L[tind]=0.5*np.pi**0.5*np.real(\
                               np.sum(np.conj(g_2[:,tind])*g_landau[:,tind])+\
                               np.sum(np.conj(g_landau[:,tind])*g_2[:,tind]))
        ent_Ldiag[tind]=0.5*np.pi**0.5*np.real(np.sum(np.conj(g_landau[:,tind])*g_landau[:,tind]))-\
                                               np.real(np.sum(ent_diag[2:,tind]))
        ent_keep_res[tind]=0.5*np.pi**0.5*np.real(\
                               np.sum(np.conj(g_knk[:,tind])*g_res[:,tind])+\
                               np.sum(np.conj(g_res[:,tind])*g_knk[:,tind]))
        ent_res[tind]=0.5*np.pi**0.5*np.real(np.sum(np.conj(g_res[:,tind])*g_res[:,tind]))
        ent_tot[tind]=0.5*np.pi**0.5*np.real(np.sum(np.conj(g_k[:,tind])*g_k[:,tind]))

        #print np.abs(np.sum(np.conj(g_k[:,tind])*g_k[:,tind])),\
        #      np.abs(np.sum(np.conj(diffr)*diffr)),\
        #      np.abs(np.sum(np.conj(al_n[:,tind])*ar_n[:,tind]))

    grs=evals[0:num_keep,1] 
    frs=evals[0:num_keep,0] 
    ents=0.5*np.pi**0.5*np.sum(np.abs(ar_n[:num_keep,:]**2),axis=1)
    #plt.scatter(frs,grs,c=ents,s=200,cmap='hsv',linewidth=0.1)
    #plt.scatter(frs,grs,c=ents,s=200,cmap='jet',linewidth=0.1)
    #plt.scatter(frs,grs,c=ents,s=200,cmap='autumn',linewidth=0.1)
    #plt.scatter(frs,grs,c=ents,s=200,cmap='RdBu',linewidth=0.1)
    #plt.scatter(frs,grs,c=ents,s=200,cmap='RdYlGn',linewidth=0.1)
    plt.scatter(frs,grs,c=np.log(ents/ents[0]),s=200,cmap='gist_rainbow',linewidth=0.1)
    for i in range(num_keep):
        plt.annotate(str(ents[i]/ents[0])[:6],(frs[i]+0.03,grs[i]))
    #s=size of symbol
    #linewidth= width of symbol outline
    #help(plt.scatter) for more info
    plt.xlabel('$\omega (v_T/R)$',size=22)
    plt.ylabel('$\gamma (v_T/R)$',size=22)
    plt.title('Entropy (normalized to ITG)')
    cb=plt.colorbar()
    cb.ax.set_ylabel("$ln(E/E_{ITG})$")
    plt.show()

    #print np.shape(ent_diag)

    plt.title("Entropy Components")
    plt.plot(time[istart:iend+1],ent_tot,color='black',label='Total Entropy')
    for evn in range(num_keep):
        plt.plot(time[istart:iend+1],ent_diag[evn,:],label='Ent mode '+str(evn+1))
    #    plt.plot(time[istart:iend+1],ent_diag[1,:],label='Ent mode 2')
    #    plt.plot(time[istart:iend+1],ent_diag[2,:],label='Ent mode 3')
    #    plt.plot(time[istart:iend+1],ent_diag[3,:],label='Ent mode 4')
    #    plt.plot(time[istart:iend+1],ent_diag[4,:],label='Ent mode 5')
    #    plt.plot(time[istart:iend+1],ent_diag[5,:],label='Ent mode 6')
    plt.plot(time[istart:iend+1],abs(ent_12),marker='+',linestyle='-.',label='Abs Ent 12 Cross')
    plt.plot(time[istart:iend+1],abs(ent_1L),marker='+',linestyle='-.',label='Abs Ent 1L Cross')
    plt.plot(time[istart:iend+1],abs(ent_2L),marker='+',linestyle='-.',label='Abs Ent 2L Cross')
    plt.plot(time[istart:iend+1],abs(ent_Ldiag),marker='+',linestyle='-.',label='Abs Ent L Cross')
    plt.plot(time[istart:iend+1],abs(ent_keep_res),marker='+',linestyle='-.',label='Abs Ent KR Cross')
    plt.plot(time[istart:iend+1],abs(ent_res),label='Abs Ent Res.')
    plt.plot(time[istart:iend+1]\
             ,np.sum(ent_diag,axis=0)+ent_12+ent_1L+ent_2L+ent_res+ent_keep_res+ent_Ldiag \
             ,marker='+'\
             ,color='black',linestyle='-.',label='Sum of All')
    plt.plot(time[istart:iend+1]\
             ,np.sum(ent_diag,axis=0)+ent_12+ent_1L+ent_2L+ent_res+ent_keep_res+ent_Ldiag-ent_tot \
             ,marker='+'\
             ,color='red',linestyle='-.',label='Diff')
    plt.legend(loc='lower right')
    plt.show()

    plt.title("ITG, DW, and Res.")
    for evn in range(2):
        plt.plot(time[istart:iend+1],ent_diag[evn,:],label='Ent mode '+str(evn+1))
    plt.plot(time[istart:iend+1],abs(ent_res),label='Abs Ent Res.')
    plt.plot(time[istart:iend+1],ent_res+ent_diag[0,:]+ent_diag[0,:],label='ITG+DW+Res')
    plt.plot(time[istart:iend+1],ent_tot,label='Total')
    plt.legend()
    plt.show()

    plt.title("Entropy Components")
    plt.semilogy(time[istart:iend+1],ent_tot,color='black',label='Total Entropy')
    for evn in range(num_keep):
        plt.semilogy(time[istart:iend+1],ent_diag[evn,:],label='Ent mode '+str(evn+1))
    #plt.semilogy(time[istart:iend+1],ent_diag[1,:],label='Ent mode 2')
    #plt.semilogy(time[istart:iend+1],ent_diag[2,:],label='Ent mode 3')
    #plt.semilogy(time[istart:iend+1],ent_diag[3,:],label='Ent mode 4')
    #plt.semilogy(time[istart:iend+1],ent_diag[4,:],label='Ent mode 5')
    #plt.semilogy(time[istart:iend+1],ent_diag[5,:],label='Ent mode 6')
    plt.semilogy(time[istart:iend+1],abs(ent_12),marker='+',linestyle='-.',label='Abs Ent 12 Cross')
    plt.semilogy(time[istart:iend+1],abs(ent_1L),marker='+',linestyle='-.',label='Abs Ent 1L Cross')
    plt.semilogy(time[istart:iend+1],abs(ent_2L),marker='+',linestyle='-.',label='Abs Ent 2L Cross')
    plt.semilogy(time[istart:iend+1],abs(ent_Ldiag),marker='+',linestyle='-.',label='Abs Ent L Cross')
    plt.semilogy(time[istart:iend+1],abs(ent_keep_res),marker='+',linestyle='-.',label='Abs Ent KR Cross')
    plt.semilogy(time[istart:iend+1],abs(ent_res),label='Abs Ent Res.')
    plt.semilogy(time[istart:iend+1]\
             ,np.sum(ent_diag,axis=0)+ent_12+ent_1L+ent_2L+ent_res+ent_keep_res+ent_Ldiag \
             ,marker='+'\
             ,color='black',linestyle='-.',label='Sum of All')
    plt.semilogy(time[istart:iend+1]\
             ,np.sum(ent_diag,axis=0)+ent_12+ent_1L+ent_2L+ent_res+ent_keep_res+ent_Ldiag-ent_tot \
             ,marker='+'\
             ,color='red',linestyle='-.',label='Diff')
    plt.legend(loc='lower right')
    plt.show()


    #plt.plot(time[istart:iend+1],(ent_12),marker='+',linestyle='-.',label='Ent 12 Cross')
    for evn in range(2,num_keep):
        plt.plot(time[istart:iend+1],ent_diag[evn,:],label='Ent mode '+str(evn+1))
    #plt.plot(time[istart:iend+1],ent_diag[3,:],label='Ent mode 4')
    #plt.plot(time[istart:iend+1],ent_diag[4,:],label='Ent mode 5')
    #plt.plot(time[istart:iend+1],ent_diag[5,:],label='Ent mode 6')
    plt.plot(time[istart:iend+1],(ent_1L),label='Ent 1L Cross')
    plt.plot(time[istart:iend+1],(ent_2L),label='Ent 2L Cross')
    plt.plot(time[istart:iend+1],(ent_Ldiag),label='Ent L Cross')
    plt.legend(loc='lower right')
    plt.show()
    
    plt.plot(time[istart:iend+1],ent_Ldiag+np.sum(ent_diag[2:,:])+ent_1L+ent_2L\
             ,label='All L')
    plt.legend(loc='lower right')
    plt.show()

    #Get Diss stuff
    diss_diag=np.empty((num_keep,ntime))
    diss_tot=np.empty(ntime)
    diss_12=np.empty(ntime)
    diss_1L=np.empty(ntime)
    diss_2L=np.empty(ntime)
    diss_Lcross=np.empty(ntime)
    diss_res=np.empty(ntime)
    diss_keep_res=np.empty(ntime)
    for t in range(istart,iend+1):
        tind=t-istart
        for evn in range(num_keep):
            diss_diag[evn,tind]=get_energy_single_k(ar_n[evn,tind]*revecs_new[:,evn],\
                                 ar_n[evn,tind]*revecs_new[:,evn],kx,ky,kz,1)                   
        diss_12[tind]=get_energy_single_k(g_1[:,tind],g_2[:,tind],kx,ky,kz,1)
        diss_12[tind]=diss_12[tind]+get_energy_single_k(g_2[:,tind],g_1[:,tind],kx,ky,kz,1)
        diss_1L[tind]=get_energy_single_k(g_1[:,tind],g_landau[:,tind],kx,ky,kz,1)
        diss_1L[tind]=diss_1L[tind]+get_energy_single_k(g_landau[:,tind],g_1[:,tind],kx,ky,kz,1)
        diss_2L[tind]=get_energy_single_k(g_2[:,tind],g_landau[:,tind],kx,ky,kz,1)
        diss_2L[tind]=diss_2L[tind]+get_energy_single_k(g_landau[:,tind],g_2[:,tind],kx,ky,kz,1)
        diss_Lcross[tind]=get_energy_single_k(g_landau[:,tind],g_landau[:,tind],kx,ky,kz,1)-\
                                               np.real(np.sum(diss_diag[2:,tind]))
        diss_keep_res[tind]=get_energy_single_k(g_knk[:,tind],g_res[:,tind],kx,ky,kz,1)+\
                            get_energy_single_k(g_res[:,tind],g_knk[:,tind],kx,ky,kz,1)
        diss_res[tind]=get_energy_single_k(g_res[:,tind],g_res[:,tind],kx,ky,kz,1)
        diss_tot[tind]=get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,1)
        #print np.abs(np.sum(np.conj(g_k[:,tind])*g_k[:,tind])),\
        #      np.abs(np.sum(np.conj(diffr)*diffr)),\
        #      np.abs(np.sum(np.conj(al_n[:,tind])*ar_n[:,tind]))

    plt.title("Dissipation Components")
    plt.plot(time[istart:iend+1],diss_tot,color='black',label='Total Diss.')
    for evn in range(num_keep):
        plt.plot(time[istart:iend+1],diss_diag[evn,:],label='Diss. mode '+str(evn+1))
    #    plt.plot(time[istart:iend+1],ent_diag[1,:],label='Ent mode 2')
    #    plt.plot(time[istart:iend+1],ent_diag[2,:],label='Ent mode 3')
    #    plt.plot(time[istart:iend+1],ent_diag[3,:],label='Ent mode 4')
    #    plt.plot(time[istart:iend+1],ent_diag[4,:],label='Ent mode 5')
    #    plt.plot(time[istart:iend+1],ent_diag[5,:],label='Ent mode 6')
    plt.plot(time[istart:iend+1],(diss_12),marker='+',linestyle='-.',label='Diss. 12 Cross')
    plt.plot(time[istart:iend+1],(diss_1L),marker='+',linestyle='-.',label='Diss. 1L Cross')
    plt.plot(time[istart:iend+1],(diss_2L),marker='+',linestyle='-.',label='Diss. 2L Cross')
    plt.plot(time[istart:iend+1],(diss_Lcross),marker='+',linestyle='-.',label='Diss. L Cross')
    plt.plot(time[istart:iend+1],(diss_keep_res),marker='+',linestyle='-.',label='Diss. KR Cross')
    plt.plot(time[istart:iend+1],(diss_res),marker='+',linestyle='-.',label='Diss. Res.')
    plt.plot(time[istart:iend+1]\
             ,np.sum(diss_diag,axis=0)+diss_12+diss_1L+diss_2L+diss_res+diss_keep_res+diss_Lcross \
             ,marker='+'\
             ,color='black',linestyle='-.',label='Sum of All')
    plt.plot(time[istart:iend+1]\
             ,np.sum(diss_diag,axis=0)+diss_12+diss_1L+diss_2L+diss_res+diss_keep_res+diss_Lcross-diss_tot \
             ,marker='+'\
             ,color='red',linestyle='-.',label='Diff')
    plt.legend(loc='lower right')
    plt.show()

    diss_diff=np.sum(diss_diag,axis=0)+diss_12+diss_1L+diss_2L+diss_res+diss_keep_res+diss_Lcross-diss_tot 
    print diss_diff

    #plt.plot(time[istart:iend+1],diss_Ldiag+np.sum(diss_diag[2:,:])+diss_1L+diss_2L\
    #         ,label='Diss. All L')
    plt.plot(time[istart:iend+1],diss_tot,label='Diss. Total')
    for evn in range(2,num_keep):
        plt.plot(time[istart:iend+1],diss_diag[evn,:],label='Diss. mode '+str(evn+1))
    plt.plot(time[istart:iend+1],diss_Lcross,label='Diss. L Cross')
    plt.plot(time[istart:iend+1],diss_1L,label='Diss. 1L ')
    plt.plot(time[istart:iend+1],diss_2L,label='Diss. 2L ')
    plt.plot(time[istart:iend+1],np.sum(diss_diag[2:,:],axis=0)+diss_Lcross+diss_1L\
             +diss_2L,marker='+',label='Diss. All L ')
    plt.legend(loc='lower right')
    plt.show()

    plt.plot(time[istart:iend+1],diss_tot,label='Diss. Total')
    #for evn in range(2):
    #    plt.plot(time[istart:iend+1],diss_diag[evn,:],label='Diss. mode '+str(evn+1))
    plt.plot(time[istart:iend+1],diss_tot-np.sum(diss_diag[:2,:],axis=0),marker='+',label='Diss. Tot-(1,2) ')
    plt.legend(loc='lower right')
    plt.show()

    diss_tot_tavg=np.sum(diss_tot)/float(ntime)
    diss_ITG_tavg=np.sum(diss_diag[0,:])/float(ntime)
    diss_DW_tavg=np.sum(diss_diag[1,:])/float(ntime)
    print "Total Time Avg dissipation:",diss_tot_tavg 
    print "ITG/Tot Time Avg dissipation:",diss_ITG_tavg/abs(diss_tot_tavg )
    print "DW/Tot Time Avg dissipation:", diss_DW_tavg/abs(diss_tot_tavg)
    print "(ITG+DW)/Tot Time Avg dissipation:", (diss_ITG_tavg+diss_DW_tavg)/abs(diss_tot_tavg)
    print "Res/Tot Time Avg dissipation:", (diss_tot_tavg-diss_ITG_tavg-diss_DW_tavg)/abs(diss_tot_tavg)


    return revecs_new,levecs_new

def ev_proj_12(time,kx,ky,kz,nev=2,start_time=-1.0,end_time=-1.0,show_plots=True):

    num_keep=2
    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    new_dir=par['diagdir'][1:len(par['diagdir'])-1]+'/ev_files/'
    ev_file='ev_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    print ev_file
    evals=np.genfromtxt(new_dir+ev_file)
    #print evals

    ###Read in left and right eigenvectors
    revec_file='revec_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    rf=open(new_dir+revec_file,'rb')
    #evf=open(new_dir+ev_file,'rb')
    revecs=np.zeros((par['nv0'],par['nv0']),dtype='complex')
    ev_num=np.array(range(nev))
    for i in range(nev):
        rf.seek(i*4+i*par['nv0']*16)
        ev_num[i]=np.fromfile(rf,dtype='int32',count=1)
        rf.seek((i+1)*4+i*par['nv0']*16)
        revecs[:,i]=np.fromfile(rf,dtype='complex128',count=par['nv0'])
        revecs[:,i]=revecs[:,i]/\
                np.sqrt(np.sum(np.conj(revecs[:,i])*revecs[:,i]))

    ortho_tol=1.0e-10
    #num_keep=6

    ru,rs,rv=np.linalg.svd(revecs[:,0:num_keep])
    #print "svals of revec",rs
    revecs_new=np.zeros((par['nv0'],par['nv0']),dtype='complex')
    revecs_new[:,0:num_keep]=revecs[:,0:num_keep]
    revecs_new[:,num_keep:]=ru[:,num_keep:]
    ru,rs,rv=np.linalg.svd(revecs_new)
    #print "svals of revec_new",rs
    levecs_new=np.linalg.pinv(revecs_new,rcond=10.0**-8)
    levecs_new=np.conj(levecs_new)

    products=np.empty((par['nv0'],par['nv0']),dtype='float')
    prod_quality=np.empty(par['nv0'])
    cutoff_test=0
    cutoff=0
    for i in range(par['nv0']):
        ###Normalize left eigenvectors
        #prod=abs(sum(np.conj(levecs[:,i])*revecs[:,i]))
        prod=(np.sum(np.conj(levecs_new[i,:])*revecs_new[:,i]))
        levecs_new[i,:]=levecs_new[i,:]/prod
        for j in range(par['nv0']):
            ###Note: Must use one of the following
            products[i,j]=np.real(np.sum(np.conj(levecs_new[i,:])*revecs_new[:,j]))
            #products[i,j]=abs(np.sum(np.conj(levecs[:,i])*revecs[:,j]))
            #if i==j:
            #    print i,j,products[i,j]
        
        prod_quality[i]=np.sum(products[i,:])-products[i,i]
        if prod_quality[i] < ortho_tol and cutoff_test==0:
            cutoff=i
        else:
            cutoff_test=1
        
    ngood=cutoff+1

    #print "Number of good Left EV's",ngood
    nv0=par['nv0']
    #xaxis=np.arange((par['nv0']))+1
    #yaxis=np.arange(cutoff+4)+1
    if show_plots:
        xaxis=np.arange(par['nv0'])+1
        yaxis=np.arange(par['nv0'])+1
        levels=np.arange(0.0,1.01,0.01)
        plt.contour(xaxis,yaxis,products,levels)
        #plt.title('Number of good (err < 10^-3) eigenvectors: '+str(cutoff+1))
        plt.title('Right vs. Left Eigenvectors from pseudo inverse')
        plt.xlabel('Right Eigenvector')
        plt.ylabel('Left Eigenvector')
        plt.colorbar()
        plt.show()
    #Get indices corresponding to k's


    if kx < 0.0:
        print "Error! kx must be positive!"
        stop
    ikx=get_kindex(kx,par['kxmin'],par['nkx0'])
    iky=get_kindex(ky,par['kymin'],par['nky0'])
    ikz=get_kindex(kz,par['kzmin'],par['nkz0'])
    print 'ikx',ikx
    print 'iky',iky
    print 'ikz',ikz

    #g_k: Total distribution function
    g_k=np.empty((par['nv0'],ntime),dtype='complex')
    #g_kr: Total distribution function reproduced from decomposition
    g_kr=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_knk: Sum of all ev's used
    g_knk=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_res: Residual distribution function: g_k-g_knk
    g_res=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_1: dist. func. for ITG mode
    g_1=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_2: dist. func. for Mode 2 - i.e. drift wave
    g_2=np.zeros((par['nv0'],ntime),dtype='complex')
    ar_n=np.empty((nv0,ntime),dtype='complex')
    for t in range(istart,iend+1):
        tind=t-istart
        print 'time=',time[t],' of ',time[iend]
        gt0=read_time_step_g(t)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        g_k[:,tind]=gt0[ikx,iky,ikz,:]
        ###Find coefficients in REV basis
        for evn in range(nv0):
            ar_n[evn,tind]=np.sum(np.conj(levecs_new[evn,:])*g_k[:,tind])
            g_kr[:,tind]=g_kr[:,tind]+ar_n[evn,tind]*revecs_new[:,evn]

        for evn in range(num_keep):        
            #Get sum of the modes we are keeping
            g_knk[:,tind]=g_knk[:,tind]+ar_n[evn,tind]*revecs_new[:,evn]

        g_1[:,tind]=ar_n[0,tind]*revecs_new[:,0]
        g_2[:,tind]=ar_n[1,tind]*revecs_new[:,1]
        #Residual distribution function
        g_res[:,tind]=g_k[:,tind]-g_knk[:,tind]
        #print np.abs(np.sum(np.conj(g_res[:,tind])*g_res[:,tind]))/np.abs(np.sum(np.conj(g_k[:,tind])*g_k[:,tind]))

    g_kr=g_1+g_2+g_res
    diff=g_k-g_kr
    #plt.plot(time[istart:iend+1],np.sum(np.abs(diff),axis=0),color='black',label='Error')
    #plt.show()
    #print np.sum(np.abs(diff),axis=0)
    diff=g_k-g_knk-g_res
    #print np.sum(np.abs(diff),axis=0)

    #Get Entropy stuff
    ent_diag=np.empty((num_keep,ntime))
    ent_tot=np.empty(ntime)
    ent_12=np.empty(ntime)
    ent_res=np.empty(ntime)
    ent_keep_res=np.empty(ntime)
    for t in range(istart,iend+1):
        tind=t-istart
        for evn in range(num_keep):
            #ent_diag[evn,tind]=0.5*np.pi**0.5*np.real(\
            #                   np.sum(np.conj(ar_n[evn,tind]*revecs_new[:,evn])*\
            #                            ar_n[evn,tind]*revecs_new[:,evn]))
            ent_diag[evn,tind]=get_energy_single_k(ar_n[evn,tind]*revecs_new[:,evn],ar_n[evn,tind]*revecs_new[:,evn],kx,ky,kz,-1)
        #ent_12[tind]=0.5*np.pi**0.5*np.real(\
        #                       np.sum(np.conj(g_1[:,tind])*g_2[:,tind])+\
        #                       np.sum(np.conj(g_2[:,tind])*g_1[:,tind]))
        ent_12[tind]=get_energy_single_k(g_1[:,tind],g_2[:,tind],kx,ky,kz,-1)
        ent_12[tind]=ent_12[tind]+get_energy_single_k(g_2[:,tind],g_1[:,tind],kx,ky,kz,-1)
        #ent_keep_res[tind]=0.5*np.pi**0.5*np.real(\
        #                       np.sum(np.conj(g_knk[:,tind])*g_res[:,tind])+\
        #                       np.sum(np.conj(g_res[:,tind])*g_knk[:,tind]))
        ent_keep_res[tind]=                   get_energy_single_k(g_knk[:,tind],g_res[:,tind],kx,ky,kz,-1)
        ent_keep_res[tind]=ent_keep_res[tind]+get_energy_single_k(g_res[:,tind],g_knk[:,tind],kx,ky,kz,-1)
        #ent_res[tind]=0.5*np.pi**0.5*np.real(np.sum(np.conj(g_res[:,tind])*g_res[:,tind]))
        ent_res[tind]=get_energy_single_k(g_res[:,tind],g_res[:,tind],kx,ky,kz,-1)
        ent_tot[tind]=get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,-1)

    grs=evals[0:num_keep,1] 
    frs=evals[0:num_keep,0] 

    if show_plots:
        plt.title("Energy Components")
        plt.plot(time[istart:iend+1],ent_tot,color='black',label='Total Energy')
        for evn in range(num_keep):
            plt.plot(time[istart:iend+1],ent_diag[evn,:],label='En mode '+str(evn+1))
        plt.plot(time[istart:iend+1],(ent_res),label='En Res.')
        plt.plot(time[istart:iend+1],(ent_12),marker='+',linestyle='-.',label='En 12 Cross')
        plt.plot(time[istart:iend+1],(ent_keep_res),marker='+',linestyle='-.',label='En KR Cross')
        plt.plot(time[istart:iend+1]\
                 ,np.sum(ent_diag,axis=0)+ent_12+ent_res+ent_keep_res\
                 ,marker='+'\
                 ,color='black',linestyle='-.',label='Sum of All')
        plt.plot(time[istart:iend+1]\
                 ,np.sum(ent_diag,axis=0)+ent_12+ent_res+ent_keep_res-ent_tot \
                 ,marker='+'\
                 ,color='red',linestyle='-.',label='Diff')
        plt.legend(loc='lower right')
        plt.show()


    #plt.title("ITG, DW, and Res.")
    #for evn in range(2):
    #    plt.plot(time[istart:iend+1],ent_diag[evn,:],label='Ent mode '+str(evn+1))
    #plt.plot(time[istart:iend+1],(ent_res),label='Ent Res.')
    #plt.plot(time[istart:iend+1],ent_res+ent_diag[0,:]+ent_diag[1,:],label='ITG+DW+Res')
    #plt.plot(time[istart:iend+1],ent_tot,label='Total')
    #plt.legend()
    #plt.show()

    #plt.title("Entropy Components")
    #plt.semilogy(time[istart:iend+1],ent_tot,color='black',label='Total Entropy')
    #for evn in range(num_keep):
    #    plt.semilogy(time[istart:iend+1],ent_diag[evn,:],label='Ent mode '+str(evn+1))
    #plt.semilogy(time[istart:iend+1],(ent_12),marker='+',linestyle='-.',label='Ent 12 Cross')
    #plt.semilogy(time[istart:iend+1],(ent_keep_res),marker='+',linestyle='-.',label='Ent KR Cross')
    #plt.semilogy(time[istart:iend+1],(ent_res),label='Ent Res.')
    #plt.semilogy(time[istart:iend+1]\
    #         ,np.sum(ent_diag,axis=0)+ent_12+ent_res+ent_keep_res\
    #         ,marker='+'\
    #         ,color='black',linestyle='-.',label='Sum of All')
    #plt.semilogy(time[istart:iend+1]\
    #         ,np.sum(ent_diag,axis=0)+ent_12+ent_res+ent_keep_res-ent_tot \
    #         ,marker='+'\
    #         ,color='red',linestyle='-.',label='Diff')
    #plt.legend(loc='lower right')
    #plt.show()


    #for evn in range(2,num_keep):
    #    plt.plot(time[istart:iend+1],ent_diag[evn,:],label='Ent mode '+str(evn+1))
    #plt.legend(loc='lower right')
    #plt.show()

    #Get QC stuff
    QC_diag=np.empty((num_keep,ntime))
    QC_tot=np.empty(ntime)
    QC_12=np.empty(ntime)
    QC_res=np.empty(ntime)
    QC_keep_res=np.empty(ntime)
    for t in range(istart,iend+1):
        tind=t-istart
        for evn in range(num_keep):
            QC_diag[evn,tind]=get_energy_single_k(ar_n[evn,tind]*revecs_new[:,evn],\
                                 ar_n[evn,tind]*revecs_new[:,evn],kx,ky,kz,1)                   
            QC_diag[evn,tind]=QC_diag[evn,tind]+get_energy_single_k(ar_n[evn,tind]*revecs_new[:,evn],\
                                 ar_n[evn,tind]*revecs_new[:,evn],kx,ky,kz,2)                   
            QC_diag[evn,tind]=QC_diag[evn,tind]+get_energy_single_k(ar_n[evn,tind]*revecs_new[:,evn],\
                                 ar_n[evn,tind]*revecs_new[:,evn],kx,ky,kz,6)                   
            QC_diag[evn,tind]=QC_diag[evn,tind]+get_energy_single_k(ar_n[evn,tind]*revecs_new[:,evn],\
                                 ar_n[evn,tind]*revecs_new[:,evn],kx,ky,kz,8)                   

        QC_12[tind]=get_energy_single_k(g_1[:,tind],g_2[:,tind],kx,ky,kz,1)
        QC_12[tind]=QC_12[tind]+get_energy_single_k(g_1[:,tind],g_2[:,tind],kx,ky,kz,2)
        QC_12[tind]=QC_12[tind]+get_energy_single_k(g_1[:,tind],g_2[:,tind],kx,ky,kz,6)
        QC_12[tind]=QC_12[tind]+get_energy_single_k(g_1[:,tind],g_2[:,tind],kx,ky,kz,8)

        QC_12[tind]=QC_12[tind]+get_energy_single_k(g_2[:,tind],g_1[:,tind],kx,ky,kz,1)
        QC_12[tind]=QC_12[tind]+get_energy_single_k(g_2[:,tind],g_1[:,tind],kx,ky,kz,2)
        QC_12[tind]=QC_12[tind]+get_energy_single_k(g_2[:,tind],g_1[:,tind],kx,ky,kz,6)
        QC_12[tind]=QC_12[tind]+get_energy_single_k(g_2[:,tind],g_1[:,tind],kx,ky,kz,8)

        QC_keep_res[tind]=get_energy_single_k(g_knk[:,tind],g_res[:,tind],kx,ky,kz,1)+\
                            get_energy_single_k(g_res[:,tind],g_knk[:,tind],kx,ky,kz,1)
        QC_keep_res[tind]=QC_keep_res[tind]+get_energy_single_k(g_knk[:,tind],g_res[:,tind],kx,ky,kz,2)+\
                            get_energy_single_k(g_res[:,tind],g_knk[:,tind],kx,ky,kz,2)
        QC_keep_res[tind]=QC_keep_res[tind]+get_energy_single_k(g_knk[:,tind],g_res[:,tind],kx,ky,kz,6)+\
                            get_energy_single_k(g_res[:,tind],g_knk[:,tind],kx,ky,kz,6)
        QC_keep_res[tind]=QC_keep_res[tind]+get_energy_single_k(g_knk[:,tind],g_res[:,tind],kx,ky,kz,8)+\
                            get_energy_single_k(g_res[:,tind],g_knk[:,tind],kx,ky,kz,8)

        QC_res[tind]=get_energy_single_k(g_res[:,tind],g_res[:,tind],kx,ky,kz,1)
        QC_res[tind]=QC_res[tind]+get_energy_single_k(g_res[:,tind],g_res[:,tind],kx,ky,kz,2)
        QC_res[tind]=QC_res[tind]+get_energy_single_k(g_res[:,tind],g_res[:,tind],kx,ky,kz,6)
        QC_res[tind]=QC_res[tind]+get_energy_single_k(g_res[:,tind],g_res[:,tind],kx,ky,kz,8)

        QC_tot[tind]=get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,1)
        QC_tot[tind]=QC_tot[tind]+get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,2)
        QC_tot[tind]=QC_tot[tind]+get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,6)
        QC_tot[tind]=QC_tot[tind]+get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,8)

        #print np.abs(np.sum(np.conj(g_k[:,tind])*g_k[:,tind])),\
        #      np.abs(np.sum(np.conj(diffr)*diffr)),\
        #      np.abs(np.sum(np.conj(al_n[:,tind])*ar_n[:,tind]))

    if show_plots:
        plt.title("QC Components")
        plt.plot(time[istart:iend+1],QC_tot,color='black',label='Total QC.')
        for evn in range(num_keep):
            plt.plot(time[istart:iend+1],QC_diag[evn,:],label='QC. mode '+str(evn+1))
        plt.plot(time[istart:iend+1],(QC_res),marker='+',label='QC. Res.')
        plt.plot(time[istart:iend+1],(QC_12),marker='+',linestyle='-.',label='QC. 12 Cross')
        plt.plot(time[istart:iend+1],(QC_keep_res),marker='+',linestyle='-.',label='QC. KR Cross')
        plt.plot(time[istart:iend+1]\
                 ,np.sum(QC_diag,axis=0)+QC_12+QC_res+QC_keep_res\
                 ,marker='+'\
                 ,color='black',linestyle='-.',label='Sum of All')
        plt.plot(time[istart:iend+1]\
                 ,np.sum(QC_diag,axis=0)+QC_12+QC_res+QC_keep_res-QC_tot \
                 ,marker='+'\
                 ,color='red',linestyle='-.',label='Diff')
        plt.legend(loc='lower right')
        plt.show()

    QC_diff=np.sum(QC_diag,axis=0)+QC_12+QC_res+QC_keep_res-QC_tot 

    if show_plots:
        plt.plot(time[istart:iend+1],QC_tot,label='QC. Total')
        plt.plot(time[istart:iend+1],QC_tot-QC_diag[0,:]-QC_diag[1,:]-QC_12,marker='+',label='QC. Tot-(1,2,diag) ')
        plt.legend(loc='lower right')
        plt.show()

    if show_plots:
        plt.plot(time[istart:iend+1],QC_tot,label='QC. Total')
        plt.plot(time[istart:iend+1],QC_tot-QC_diag[0,:],marker='+',label='QC. Tot-QC 1 ')
        plt.legend(loc='lower right')
        plt.show()

    #print "Calculating diss averages."
    QC_tot_tavg=np.sum(QC_tot)/float(ntime)
    QC_ITG_tavg=np.sum(QC_diag[0,:])/float(ntime)
    QC_DW_tavg=np.sum(QC_diag[1,:])/float(ntime)
    QC_12_tavg=np.sum(QC_12)/float(ntime)
    QC_res_tavg=np.sum(QC_res)/float(ntime)
    QC_keep_res_tavg=np.sum(QC_keep_res)/float(ntime)

    #print "Calculating ent averages."
    ent_tot_tavg=np.sum(ent_tot)/float(ntime)
    ent_ITG_tavg=np.sum(ent_diag[0,:])/float(ntime)
    ent_DW_tavg=np.sum(ent_diag[1,:])/float(ntime)
    ent_12_tavg=np.sum(ent_12)/float(ntime)
    ent_res_tavg=np.sum(ent_res)/float(ntime)
    ent_keep_res_tavg=np.sum(ent_keep_res)/float(ntime)

    return ent_tot_tavg,ent_ITG_tavg,ent_DW_tavg,ent_12_tavg,ent_res_tavg,ent_keep_res_tavg,\
          QC_tot_tavg,QC_ITG_tavg,QC_DW_tavg,QC_12_tavg,QC_res_tavg,QC_keep_res_tavg

def ev_proj_scan_12(start_time=-1.0,end_time=-1.0):
    time=get_time_from_gout()
    #First ky scan
    ev_proj_ky=np.zeros((par['nky0']/4-1,12))
    kygrid=np.zeros(par['nky0']/4-1)
    for i in range(par['nky0']/4-1):
        ky=par['kymin']*2*(i+1)
        print "ky",ky
        kygrid[i]=ky
        kz=0.2
        ev_proj_ky[i,:]=ev_proj_12(time,0.0,ky,kz,start_time=start_time,end_time=end_time,show_plots=0)

    ev_proj_kz=np.zeros((par['nkz0']/4-1,12))
    kzgrid=np.zeros(par['nkz0']/4-1)
    for i in range(par['nkz0']/4-1):
        kz=par['kzmin']*2*(i+1)
        print "kz",kz
        kzgrid[i]=kz
        ky=0.48
        ev_proj_kz[i,:]=ev_proj_12(time,0.0,ky,kz,start_time=start_time,end_time=end_time,show_plots=0)

    ev_proj12_scan_plot(kygrid,ev_proj_ky,kzgrid,ev_proj_kz)

    return kygrid,ev_proj_ky,kzgrid,ev_proj_kz

def QC_scan(start_time=-1.0,end_time=-1.0):
    time=get_time_from_gout()
    #First ky scan
    ev_proj_ky=np.zeros((par['nky0']/4-1,6))
    kygrid=np.zeros(par['nky0']/4-1)
    for i in range(par['nky0']/4-1):
        ky=par['kymin']*2*(i+1)
        print "ky",ky
        kygrid[i]=ky
        kz=0.2
        ev_proj_ky[i,:]=QC_de(time,0.0,ky,kz,start_time=start_time,end_time=end_time,show_plots=0)

    ev_proj_kz=np.zeros((par['nkz0']/4-1,6))
    kzgrid=np.zeros(par['nkz0']/4-1)
    for i in range(par['nkz0']/4-1):
        kz=par['kzmin']*2*(i+1)
        print "kz",kz
        kzgrid[i]=kz
        ky=0.48
        ev_proj_kz[i,:]=QC_de(time,0.0,ky,kz,start_time=start_time,end_time=end_time,show_plots=0)

    #QC_scan_plot(kygrid,ev_proj_ky,kzgrid,ev_proj_kz)

    return kygrid,ev_proj_ky,kzgrid,ev_proj_kz

def QC_scan_plot(kygrid,ev_proj_ky,kzgrid,ev_proj_kz):

    plt.plot(kygrid,ev_proj_ky[:,0],label='QC_tot')
    plt.plot(kygrid,ev_proj_ky[:,1],label='QC1')
    plt.plot(kygrid,ev_proj_ky[:,2],label='QC_res')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.show()

    plt.plot(kygrid,ev_proj_ky[:,3],label='gamma_nl')
    plt.plot(kygrid,ev_proj_ky[:,4],label='gamma_lin')
    plt.plot(kygrid,ev_proj_ky[:,5],label='gamma_mode')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.show()

    plt.plot(kzgrid,ev_proj_kz[:,0],label='QC_tot')
    plt.plot(kzgrid,ev_proj_kz[:,1],label='QC1')
    plt.plot(kzgrid,ev_proj_kz[:,2],label='QC_res')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_z \rho_i$',size=18)
    plt.show()

    plt.plot(kzgrid,ev_proj_kz[:,3],label='gamma_nl')
    plt.plot(kzgrid,ev_proj_kz[:,4],label='gamma_lin')
    plt.plot(kzgrid,ev_proj_kz[:,5],label='gamma_mode')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_z \rho_i$',size=18)
    plt.show()

def ev_proj12_scan_plot(kygrid,ev_proj_ky,kzgrid,ev_proj_kz):

    plt.plot(kygrid,ev_proj_ky[:,0],label='Total Energy')
    plt.plot(kygrid,ev_proj_ky[:,1],label='Energy ITG')
    plt.plot(kygrid,ev_proj_ky[:,2]+ev_proj_ky[:,3],label='Energy DW+C1')
    plt.plot(kygrid,ev_proj_ky[:,4]+ev_proj_ky[:,5],label='Energy Landau+C12')
    plt.plot(kygrid,np.sum(ev_proj_ky[:,1:6],axis=1),'+',color='black',label='Sum')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.show()

    plt.semilogy(kygrid,ev_proj_ky[:,0],label='Total Energy')
    plt.semilogy(kygrid,ev_proj_ky[:,1],label='Energy ITG')
    plt.semilogy(kygrid,ev_proj_ky[:,2]+ev_proj_ky[:,3],label='Energy DW+C1')
    plt.semilogy(kygrid,ev_proj_ky[:,4]+ev_proj_ky[:,5],label='Energy Landau+C12')
    plt.semilogy(kygrid,np.sum(ev_proj_ky[:,1:6],axis=1),'+',color='black',label='Sum')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.show()

    plt.plot(kygrid,ev_proj_ky[:,6],label='Total QC')
    plt.plot(kygrid,ev_proj_ky[:,7],label='QC ITG')
    plt.plot(kygrid,ev_proj_ky[:,8]+ev_proj_ky[:,9],label='QC DW+C1')
    plt.plot(kygrid,ev_proj_ky[:,10]+ev_proj_ky[:,11],label='QC Landau+C12')
    plt.plot(kygrid,np.sum(ev_proj_ky[:,7:],axis=1),'+',color='black',label='Sum')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.show()

    plt.plot(kygrid,ev_proj_ky[:,6]-ev_proj_ky[:,7],label='Total QC - QC ITG')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_y \rho_i$',size=18)
    plt.show()

    plt.plot(kygrid,ev_proj_ky[:,6]-ev_proj_ky[:,7],label='Total non-ITG QC')
    plt.plot(kygrid,ev_proj_ky[:,10]+ev_proj_ky[:,11],label='QC Landau+C12')
    plt.plot(kygrid,ev_proj_ky[:,8]+ev_proj_ky[:,9],label='QC DW+C1')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_y R$',size=18)
    plt.show()


    ######Now kz
    ######Now kz
    ######Now kz

    plt.plot(kzgrid,ev_proj_kz[:,0],label='Total Energy')
    plt.plot(kzgrid,ev_proj_kz[:,1],label='Energy ITG')
    plt.plot(kzgrid,ev_proj_kz[:,2]+ev_proj_kz[:,3],label='Energy DW+C1')
    plt.plot(kzgrid,ev_proj_kz[:,4]+ev_proj_kz[:,5],label='Energy Landau+C12')
    plt.plot(kzgrid,np.sum(ev_proj_kz[:,1:6],axis=1),'+',color='black',label='Sum')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_z R$',size=18)
    plt.show()

    plt.semilogy(kzgrid,ev_proj_kz[:,0],label='Total Energy')
    plt.semilogy(kzgrid,ev_proj_kz[:,1],label='Energy ITG')
    plt.semilogy(kzgrid,ev_proj_kz[:,2]+ev_proj_kz[:,3],label='Energy DW+C1')
    plt.semilogy(kzgrid,ev_proj_kz[:,4]+ev_proj_kz[:,5],label='Energy Landau+C12')
    plt.semilogy(kzgrid,np.sum(ev_proj_kz[:,1:6],axis=1),'+',color='black',label='Sum')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_z R$',size=18)
    plt.show()

    plt.plot(kzgrid,ev_proj_kz[:,6],label='Total QC')
    plt.plot(kzgrid,ev_proj_kz[:,7],label='QC ITG')
    plt.plot(kzgrid,ev_proj_kz[:,8]+ev_proj_kz[:,9],label='QC DW+C1')
    plt.plot(kzgrid,ev_proj_kz[:,10]+ev_proj_kz[:,11],label='QC Landau+C12')
    plt.plot(kzgrid,np.sum(ev_proj_kz[:,7:],axis=1),'+',color='black',label='Sum')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_z R$',size=18)
    plt.show()

    plt.plot(kzgrid,ev_proj_kz[:,6]-ev_proj_kz[:,7],label='Total non-ITG QC')
    plt.plot(kzgrid,ev_proj_kz[:,10]+ev_proj_kz[:,11],label='QC Landau+C12')
    plt.plot(kzgrid,ev_proj_kz[:,8]+ev_proj_kz[:,9],label='QC DW+C1')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_z R$',size=18)
    plt.show()

    plt.semilogy(kzgrid,-ev_proj_kz[:,6],label='Total QC')
    plt.semilogy(kzgrid,-ev_proj_kz[:,7],label='QC ITG')
    plt.semilogy(kzgrid,-ev_proj_kz[:,8]+ev_proj_kz[:,9],label='QC DW+C1')
    plt.semilogy(kzgrid,-ev_proj_kz[:,10]+ev_proj_kz[:,11],label='QC Landau+C12')
    plt.semilogy(kzgrid,-np.sum(ev_proj_kz[:,7:],axis=1),'+',color='black',label='Sum')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_z R$',size=18)
    plt.show()

    plt.plot(kzgrid,ev_proj_kz[:,6]-ev_proj_kz[:,7],label='Total QC - QC ITG')
    plt.legend(loc='upper right')
    plt.xlabel(r'$k_z \rho_i$',size=18)
    plt.show()

def ev_tests(kx,ky,kz):

    ###Read in left and right eigenvectors
    levec_file='levec_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    revec_file='revec_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    new_dir=par['diagdir'][1:len(par['diagdir'])-1]+'/ev_files/'
    lf=open(new_dir+levec_file,'rb')
    rf=open(new_dir+revec_file,'rb')
    revecs=np.empty((par['nv0'],par['nv0']),dtype='complex')
    levecs=np.empty((par['nv0'],par['nv0']),dtype='complex')
    ev_num=np.array(range(par['nv0']))
    for i in range(par['nv0']):
        lf.seek(i*4+i*par['nv0']*16)
        ev_num[i]=np.fromfile(lf,dtype='int32',count=1)
        #print ev_num[i]
        rf.seek((i+1)*4+i*par['nv0']*16)
        lf.seek((i+1)*4+i*par['nv0']*16)
        input=np.fromfile(lf,dtype='complex128',count=par['nv0'])
        #print type(input)
        #print len(input)
        levecs[:,i]=input
        revecs[:,i]=np.fromfile(rf,dtype='complex128',count=par['nv0'])
        #normalize revecs
        revecs[:,i]=revecs[:,i]/np.sqrt(np.sum(np.conj(revecs[:,i])*revecs[:,i]))

    ortho_tol=1.0e-10
    num_keep=6

    ru,rs,rv=np.linalg.svd(revecs[:,0:num_keep])
    #print "svals of revec",rs
    revecs_new=np.empty((par['nv0'],par['nv0']))
    revecs_new[:,0:num_keep]=revecs[:,0:num_keep]
    revecs_new[:,num_keep:]=ru[:,num_keep:]
    ru,rs,rv=np.linalg.svd(revecs_new)
    print "svals of revec_new",rs
    levecs_new=np.linalg.pinv(revecs_new,rcond=10.0**-8)

    products=np.empty((par['nv0'],par['nv0']),dtype='float')
    prod_quality=np.empty(par['nv0'])
    cutoff_test=0
    cutoff=0
    for i in range(par['nv0']):
        ###Normalize left eigenvectors
        #prod=abs(sum(np.conj(levecs[:,i])*revecs[:,i]))
        prod=np.abs(np.sum((levecs_new[i,:])*revecs_new[:,i]))
        levecs_new[i,:]=levecs_new[i,:]/prod
        for j in range(par['nv0']):
            ###Note: Must use one of the following
            products[i,j]=np.abs(np.sum(levecs_new[i,:]*revecs_new[:,j]))
            #products[i,j]=abs(np.sum(np.conj(levecs[:,i])*revecs[:,j]))
        
        prod_quality[i]=np.sum(products[i,:])-products[i,i]
        if prod_quality[i] < ortho_tol and cutoff_test==0:
            cutoff=i
        else:
            cutoff_test=1
        
    num_good_evs=cutoff+1

    print "Number of good Left EV's",num_good_evs
    #xaxis=np.arange((par['nv0']))+1
    #yaxis=np.arange(cutoff+4)+1
    xaxis=np.arange(par['nv0'])+1
    yaxis=np.arange(par['nv0'])+1
    levels=np.arange(0.0,1.01,0.01)
    plt.contour(xaxis,yaxis,products,levels)
    #plt.title('Number of good (err < 10^-3) eigenvectors: '+str(cutoff+1))
    plt.title('Right vs. Left Eigenvectors from pseudo inverse')
    plt.xlabel('Right Eigenvector')
    plt.ylabel('Left Eigenvector')
    plt.colorbar()
    plt.show()
    return products
    stop

        
    #plt.plot(abs(revecs[:,0]))     
    #plt.show()
    
    products=np.empty((par['nv0'],par['nv0']),dtype='float')
    prod_quality=np.empty(par['nv0'])
    cutoff_test=0
    cutoff=0
    for i in range(par['nv0']):
        ###Normalize left eigenvectors
        #prod=abs(sum(np.conj(levecs[:,i])*revecs[:,i]))
        prod=abs(sum((levecs[:,i])*revecs[:,i]))
        levecs[:,i]=levecs[:,i]/prod
        for j in range(par['nv0']):
            ###Note: Must use one of the following
            products[i,j]=abs(np.sum(levecs[:,i]*revecs[:,j]))
            #products[i,j]=abs(np.sum(np.conj(levecs[:,i])*revecs[:,j]))
        
        prod_quality[i]=np.sum(products[i,:])-products[i,i]
        if prod_quality[i] < ortho_tol and cutoff_test==0:
            cutoff=i
        else:
            cutoff_test=1
        
    num_good_evs=cutoff+1

    print "Number of good Left EV's",num_good_evs
    #xaxis=np.arange((par['nv0']))+1
    #yaxis=np.arange(cutoff+4)+1
    xaxis=np.arange(par['nv0'])+1
    yaxis=np.arange(par['nv0'])+1
    levels=np.arange(0.0,1.01,0.01)
    plt.contour(xaxis,yaxis,products,levels)
    #plt.title('Number of good (err < 10^-3) eigenvectors: '+str(cutoff+1))
    plt.title('Right vs. Left Eigenvectors from SLEPc')
    plt.xlabel('Right Eigenvector')
    plt.ylabel('Left Eigenvector')
    plt.colorbar()
    plt.show()

    ###Orthogonalize right evecs
    ###Note: Doesn't work - not high enough precision
    gs_orthogonalize(revecs,gs_test='T')

    ###Get orthogonal basis for remaining degrees of freedom
    ru,rs,rv=np.linalg.svd(revecs[:,:num_good_evs])
    print "From svd of good evs:"
    print "svals",rs
    print np.shape(ru)
    print np.shape(rv)
    ortho_test=np.empty((par['nv0'],par['nv0']))
    for i in range(np.shape(ru)[0]):
        for j in range(np.shape(ru)[1]):
            ortho_test[i,j]=np.abs(np.sum(np.conj(ru[:,i])*ru[:,j]))
            ###Note: orthogonality with conj (as above)
            ###Note: Two scalar products below don't work
            ###ortho_test[i,j]=np.abs(np.dot(rv[:,i],rv[:,j]))
            ###ortho_test[i,j]=np.abs(np.sum((rv[:,i])*rv[:,j]))

    #print ortho_test[:3,:3]
    #plt.contour(ortho_test,200)
    #plt.colorbar()
    #plt.show()

    #Now fill out bases with svd basis
    revecs_new=np.empty((par['nv0'],par['nv0']),dtype='complex')
    levecs_new=np.empty((par['nv0'],par['nv0']),dtype='complex')
    revecs_new[:,0:num_good_evs]=revecs[:,0:num_good_evs]
    levecs_new[:,0:num_good_evs]=levecs[:,0:num_good_evs]
    revecs_new[:,num_good_evs:par['nv0']]=ru[:,num_good_evs:par['nv0']]
    levecs_new[:,num_good_evs:par['nv0']]=np.conj(ru[:,num_good_evs:par['nv0']])

    ###############Test##############
    ###############Test##############
    ###############Test##############
    products=np.empty((par['nv0'],par['nv0']),dtype='float')
    prod_quality=np.empty(par['nv0'])
    cutoff_test=0
    for i in range(par['nv0']):
        for j in range(par['nv0']):
            products[i,j]=abs(np.sum(levecs_new[:,i]*revecs_new[:,j]))
            ###Note: using conj does not produce orthogonality
            ###products[i,j]=abs(np.sum(np.conj(levecs[:,i])*revecs[:,j]))
        
        prod_quality[i]=np.sum(products[i,:])-products[i,i]
        if prod_quality[i] < ortho_tol and cutoff_test==0:
            cutoff=i
        else:
            cutoff_test=1
        
    num_good_evs=cutoff+1

    print "Number of good Left EV's",num_good_evs
    xaxis=np.arange((par['nv0']))+1
    yaxis=np.arange(cutoff+4)+1
    levels=np.arange(0.0,1.01,0.01)
    #plt.contour(xaxis,yaxis,products[0:cutoff+4,:],levels)
    plt.contour(products,levels)
    #plt.title('Number of good (err < 10^-3) eigenvectors: '+str(cutoff+1))
    plt.title('Good left and right evecs, plus orthogonal complement')
    plt.xlabel('Right Eigenvector')
    plt.ylabel('Left Eigenvector')
    plt.colorbar()
    plt.show()

    ###############Test##############
    ###############Test##############
    ###############Test##############

    #Now try pseudo-inverse
    ru,rs,rv=np.linalg.svd(revecs) 
    print "svals",rs
    print "svals first/last",rs[0]/rs[par['nv0']-1]
    keep_going=1
    i=0
    while keep_going:
        if rs[i]/rs[0] < 1.0e-10:
            rcond0=rs[i]
            num_lin_ind=i
            keep_going=0
        else:
            pass
        i=i+1

    print

    ###Get pseudo inverse of right eigenvectors
    levecs_pi=np.linalg.pinv(revecs,rcond=rcond0)
    #levecs_pi=np.linalg.pinv(revecs)
    print np.shape(levecs_pi)
    print np.info(levecs_pi)

    ###############Test##############
    ###############Test##############
    ###############Test##############
    products=np.empty((par['nv0'],par['nv0']),dtype='float')
    prod_quality_rdl=np.empty(par['nv0'])
    prod_quality_ldr=np.empty(par['nv0'])
    cutoff_test=0
    cutoff_ldr=0
    cutoff_rdl=0
    for i in range(par['nv0']):
        for j in range(par['nv0']):
            products[i,j]=abs(np.sum(levecs_pi[i,:]*revecs[:,j]))
            #products[i,j]=abs(np.sum(np.conj(levecs_pi[i,:])*revecs[:,j]))
        prod_quality_ldr[i]=np.sum(products[i,:])-products[i,i]
        if prod_quality_ldr[i] < ortho_tol and cutoff_test==0:
            cutoff_ldr=i
        else:
            cutoff_test=1

    cutoff_test=0
    for i in range(par['nv0']):
        prod_quality_rdl[i]=np.sum(products[:,i])-products[i,i]
        if prod_quality_rdl[i] < ortho_tol and cutoff_test==0:
            cutoff_rdl=i
        else:
            cutoff_test=1
        
    num_good_ldr=cutoff_ldr+1
    num_good_rdl=cutoff_rdl+1
    print "Number of good Left EV's",num_good_ldr
    print "Number of good Right EV's",num_good_rdl
    print "prod_quality_rdl",prod_quality_rdl[0:num_good_rdl+1]
    print "prod_quality_ldr",prod_quality_ldr[0:num_good_ldr+1]

    xaxis=np.arange((par['nv0']))+1
    yaxis=np.arange((par['nv0']))+1
    levels=np.arange(0.0,1.01,0.01)
    #plt.contour(xaxis,yaxis,products[0:cutoff+4,:],levels)
    plt.contour(xaxis,yaxis,products,levels)
    #plt.title('Number of good (err < 10^-3) eigenvectors: '+str(num_good_ldr)\
    #          +","+str(num_good_rdl))
    plt.title('Left eigenvectors with pseudo-inverse')
    plt.xlabel('Right Eigenvector')
    plt.ylabel('Left Eigenvector')
    plt.colorbar()
    plt.show()

    ###############Test##############
    ###############Test##############
    ###############Test##############

def gs_orthogonalize(vecs_in,normalize=1.0,gs_test=0.0):
   """Function for diagonalizing a set of vectors. Each column must be a vector"""
   if np.ndim(vecs_in) != 2:
       print "Error gs_orthogonalize: must input a 2 dimensional array."
       stop
   if np.shape(vecs_in)[0] < np.shape(vecs_in)[1]:
       print "Error gs_orthogonalize: Must have more elements than vectors."

   nev=np.shape(vecs_in)[1]
   sp=np.empty(nev,dtype='complex')
   #Check for divide by zero
   for i in range(nev):
       sp[i]=np.sqrt(np.sum((np.conj(vecs_in[:,i])*\
                                            vecs_in[:,i])))

   gs_vecs=np.zeros(np.shape(vecs_in),dtype='complex')
   if np.abs(sp[0]) > 1.0e-14:
       gs_vecs[:,0]=vecs_in[:,0]/sp[0]
   else:
       print "Error!"
       stop
   #print "magnitude of normalized vec 1:"
   #print np.sum((np.conj(gs_vecs[:,0])* gs_vecs[:,0]))

   for i in range(nev):
       gs_vecs[:,i]=vecs_in[:,i]/sp[i]
       #print i, sp[i],np.sqrt(np.sum((np.conj(gs_vecs[:,i])*\
                                            #gs_vecs[:,i])))
       
   for i in range(nev):
       for j in range(i):
           ip=np.sum(np.conj(gs_vecs[:,j])*gs_vecs[:,i])
           #print i,j,ip
           gs_vecs[:,i]=gs_vecs[:,i]-gs_vecs[:,j]*ip
       if normalize:
           gs_vecs[:,i]=gs_vecs[:,i]/\
                   np.sqrt(np.sum(np.conj(gs_vecs[:,i])*gs_vecs[:,i]))
     
   orth_test=np.empty((nev,nev))
   if gs_test:
       for i in range(nev):
           for j in range(nev):
               orth_test[i,j]=np.abs(np.sum(np.conj(gs_vecs[:,i])*gs_vecs[:,j]))
       #print orth_test
       pltaxis=np.arange(nev)+1
       plt.contour(pltaxis,pltaxis,orth_test,200)
       plt.title('Orthogonalized Vectors')
       #plt.contour(orth_test[:,:20],200)
       plt.colorbar()
       plt.show()

def ev_proj_1(time,kx,ky,kz,start_time=-1.0,end_time=-1.0,show_plots=True):

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    new_dir=par['diagdir'][1:len(par['diagdir'])-1]+'/ev_files/'
    ev_file='ev_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    print ev_file
    evals=np.genfromtxt(new_dir+ev_file)
    #print evals

    ###Read in left and right eigenvectors
    revec_file='revec_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    rf=open(new_dir+revec_file,'rb')
    #evf=open(new_dir+ev_file,'rb')
    revecs=np.zeros((par['nv0'],par['nv0']),dtype='complex')
    ev_num=np.array(range(par['nv0']))
    for i in range(2):
        rf.seek(i*4+i*par['nv0']*16)
        ev_num[i]=np.fromfile(rf,dtype='int32',count=1)
        rf.seek((i+1)*4+i*par['nv0']*16)
        revecs[:,i]=np.fromfile(rf,dtype='complex128',count=par['nv0'])
        revecs[:,i]=revecs[:,i]/\
                np.sqrt(np.sum(np.conj(revecs[:,i])*revecs[:,i]))

    ortho_tol=1.0e-10

    #num_keep = 2 for getting left eigenvectors
    num_keep=2
    ru,rs,rv=np.linalg.svd(revecs[:,0:num_keep])
    #print "svals of revec",rs
    revecs_new=np.zeros((par['nv0'],par['nv0']),dtype='complex')
    revecs_new[:,0:num_keep]=revecs[:,0:num_keep]
    revecs_new[:,num_keep:]=ru[:,num_keep:]
    ru,rs,rv=np.linalg.svd(revecs_new)
    #print "svals of revec_new",rs
    levecs_new=np.linalg.pinv(revecs_new,rcond=10.0**-8)
    levecs_new=np.conj(levecs_new)

    products=np.empty((par['nv0'],par['nv0']),dtype='float')
    prod_quality=np.empty(par['nv0'])
    cutoff_test=0
    cutoff=0
    for i in range(par['nv0']):
        ###Normalize left eigenvectors
        #prod=abs(sum(np.conj(levecs[:,i])*revecs[:,i]))
        prod=(np.sum(np.conj(levecs_new[i,:])*revecs_new[:,i]))
        levecs_new[i,:]=levecs_new[i,:]/prod
        for j in range(par['nv0']):
            ###Note: Must use one of the following
            products[i,j]=np.real(np.sum(np.conj(levecs_new[i,:])*revecs_new[:,j]))
            #products[i,j]=abs(np.sum(np.conj(levecs[:,i])*revecs[:,j]))
        
        prod_quality[i]=np.sum(products[i,:])-products[i,i]
        if prod_quality[i] < ortho_tol and cutoff_test==0:
            cutoff=i
        else:
            cutoff_test=1
        
    ngood=cutoff+1

    #print "Number of good Left EV's",ngood
    nv0=par['nv0']
    #xaxis=np.arange((par['nv0']))+1
    #yaxis=np.arange(cutoff+4)+1
    if show_plots:
        xaxis=np.arange(par['nv0'])+1
        yaxis=np.arange(par['nv0'])+1
        levels=np.arange(0.0,1.01,0.01)
        plt.contour(xaxis,yaxis,products,levels)
        #plt.title('Number of good (err < 10^-3) eigenvectors: '+str(cutoff+1))
        plt.title('Right vs. Left Eigenvectors from pseudo inverse')
        plt.xlabel('Right Eigenvector')
        plt.ylabel('Left Eigenvector')
        plt.colorbar()
        plt.show()
    #Get indices corresponding to k's

    if kx < 0.0:
        print "Error! kx must be positive!"
        stop
    ikx=get_kindex(kx,par['kxmin'],par['nkx0'])
    iky=get_kindex(ky,par['kymin'],par['nky0'])
    ikz=get_kindex(kz,par['kzmin'],par['nkz0'])
    print 'ikx',ikx
    print 'iky',iky
    print 'ikz',ikz

    #Now change num_keep to 1 to only get most unstable mode
    num_keep=1
    #g_k: Total distribution function
    g_k=np.empty((par['nv0'],ntime),dtype='complex')
    #g_kr: Total distribution function reproduced from decomposition
    g_kr=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_res: Residual distribution function: g_k-g_knk
    g_res=np.zeros((par['nv0'],ntime),dtype='complex')
    #g_1: dist. func. for ITG mode
    g_1=np.zeros((par['nv0'],ntime),dtype='complex')
    ar_n=np.empty((nv0,ntime),dtype='complex')
    for t in range(istart,iend+1):
        tind=t-istart
        print 'time=',time[t],' of ',time[iend]
        gt0=read_time_step_g(t)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        g_k[:,tind]=gt0[ikx,iky,ikz,:]
        ###Find coefficients in REV basis
        for evn in range(nv0):
            ar_n[evn,tind]=np.sum(np.conj(levecs_new[evn,:])*g_k[:,tind])
            g_kr[:,tind]=g_kr[:,tind]+ar_n[evn,tind]*revecs_new[:,evn]


        g_1[:,tind]=ar_n[0,tind]*revecs_new[:,0]
        #Residual distribution function
        g_res[:,tind]=g_k[:,tind]-g_1[:,tind]
        #print np.abs(np.sum(np.conj(g_res[:,tind])*g_res[:,tind]))/np.abs(np.sum(np.conj(g_k[:,tind])*g_k[:,tind]))

    g_kr=g_1+g_res
    diff=g_k-g_1-g_res
    print "Diff",np.sum(np.abs(diff),axis=0)

    #Get Energy stuff
    en_tot=np.empty(ntime)
    en_1=np.empty(ntime)
    en_res=np.empty(ntime)
    for t in range(istart,iend+1):
        tind=t-istart
        en_1[tind]=get_energy_single_k(g_1[:,tind],g_1[:,tind],kx,ky,kz,-1)
        en_tot[tind]=get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,-1)
        en_res[tind]=en_tot[tind]-en_1[tind]

    grs=evals[0:num_keep,1] 
    frs=evals[0:num_keep,0] 

    if show_plots:
        plt.title("Energy")
        plt.plot(time[istart:iend+1],en_tot,color='black',label='Total Energy')
        plt.plot(time[istart:iend+1],(en_res),label='Energy Res.')
        plt.plot(time[istart:iend+1],(en_1),marker='+',linestyle='-.',label='Energy 1')
        plt.plot(time[istart:iend+1]\
                 ,en_1+en_res,marker='+'\
                 ,linestyle='-.',label='Sum of All')
        plt.legend(loc='lower right')
        plt.show()

    #Get Q+C stuff
    QC_tot=np.empty(ntime)
    QC_1=np.empty(ntime)
    QC_res=np.empty(ntime)
    for t in range(istart,iend+1):
        tind=t-istart
        QC_1[tind]=get_energy_single_k(g_1[:,tind],g_1[:,tind],kx,ky,kz,6)
        QC_1[tind]=QC_1[tind]+get_energy_single_k(g_1[:,tind],g_1[:,tind],kx,ky,kz,1)
        QC_1[tind]=QC_1[tind]+get_energy_single_k(g_1[:,tind],g_1[:,tind],kx,ky,kz,2)
        QC_1[tind]=QC_1[tind]+get_energy_single_k(g_1[:,tind],g_1[:,tind],kx,ky,kz,8)
        QC_tot[tind]=get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,6)
        QC_tot[tind]=QC_tot[tind]+get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,1)
        QC_tot[tind]=QC_tot[tind]+get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,2)
        QC_tot[tind]=QC_tot[tind]+get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,8)
        QC_res[tind]=QC_tot[tind]-QC_1[tind]

    #if show_plots:
    #    plt.plot(time[istart:iend+1],QC_tot,color='black',label='Total Q+C')
    #    plt.legend()
    #    plt.show()
    #    plt.plot(time[istart:iend+1],QC_1,color='black',label='Q+C one')
    #    plt.legend()
    #    plt.show()
    #    plt.plot(time[istart:iend+1],QC_res,color='black',label='Q+C res')
    #    plt.legend()
    #    plt.show()

    if show_plots:
        plt.title("Q+C")
        plt.plot(time[istart:iend+1],QC_tot,color='black',label='Total Q+C')
        plt.plot(time[istart:iend+1],QC_res,label='Q+C Res.')
        plt.plot(time[istart:iend+1],QC_1,marker='+',linestyle='-.',label='Q+C One')
        plt.plot(time[istart:iend+1]\
                 ,QC_1+QC_res,marker='+'\
                 ,linestyle='-.',label='Sum of All')
        plt.legend(loc='upper left')
        plt.show()

    gamlin=np.empty(ntime)
    gamlin[:]=2.0*grs[0]
    if show_plots:
        plt.title("Gammas")
        plt.plot(time[istart:iend+1],QC_tot/en_tot,color='black',label='Total Gam')
        plt.plot(time[istart:iend+1],(QC_res/en_res),label='Gam Res.')
        plt.plot(time[istart:iend+1],(QC_1/en_1),marker='+',linestyle='-.',label='Gam 1')
        plt.plot(time[istart:iend+1],gamlin,marker='+',linestyle='-.',label='Gam lin.')
        plt.legend(loc='lower right')
        plt.show()


    #print "Calculating ent averages."
    #ent_tot_tavg=np.sum(ent_tot)/float(ntime)
    #ent_ITG_tavg=np.sum(ent_diag[0,:])/float(ntime)
    #ent_DW_tavg=np.sum(ent_diag[1,:])/float(ntime)
    #ent_12_tavg=np.sum(ent_12)/float(ntime)
    #ent_res_tavg=np.sum(ent_res)/float(ntime)
    #ent_keep_res_tavg=np.sum(ent_keep_res)/float(ntime)

    #return ent_tot_tavg,ent_ITG_tavg,ent_DW_tavg,ent_12_tavg,ent_res_tavg,ent_keep_res_tavg,\
    #      diss_tot_tavg,diss_ITG_tavg,diss_DW_tavg,diss_12_tavg,diss_res_tavg,diss_keep_res_tavg

def QC_de(time,kx,ky,kz,start_time=-1.0,end_time=-1.0,show_plots=True):

    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    new_dir=par['diagdir'][1:len(par['diagdir'])-1]+'/ev_files/'
    ev_file='ev_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    print ev_file
    evals=np.genfromtxt(new_dir+ev_file)
    #print evals

    ###Read in left and right eigenvectors
    revec_file='revec_kx'+str(kx)+'ky'+str(ky)+'kz'+str(kz)
    rf=open(new_dir+revec_file,'rb')
    #evf=open(new_dir+ev_file,'rb')
    revecs=np.zeros((par['nv0'],par['nv0']),dtype='complex')
    ev_num=np.array(range(par['nv0']))
    for i in range(1):
        rf.seek(i*4+i*par['nv0']*16)
        ev_num[i]=np.fromfile(rf,dtype='int32',count=1)
        rf.seek((i+1)*4+i*par['nv0']*16)
        revecs[:,i]=np.fromfile(rf,dtype='complex128',count=par['nv0'])
        revecs[:,i]=revecs[:,i]/\
                np.sqrt(np.sum(np.conj(revecs[:,i])*revecs[:,i]))

    ortho_tol=1.0e-10

    if kx < 0.0:
        print "Error! kx must be positive!"
        stop
    ikx=get_kindex(kx,par['kxmin'],par['nkx0'])
    iky=get_kindex(ky,par['kymin'],par['nky0'])
    ikz=get_kindex(kz,par['kzmin'],par['nkz0'])
    print 'ikx',ikx
    print 'iky',iky
    print 'ikz',ikz

    g_k=np.empty((par['nv0'],ntime),dtype='complex')
    for t in range(istart,iend+1):
        tind=t-istart
        print 'time=',time[t],' of ',time[iend]
        gt0=read_time_step_g(t)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        g_k[:,tind]=gt0[ikx,iky,ikz,:]

    #Get Energy stuff
    en_tot=np.empty(ntime)
    for t in range(istart,iend+1):
        tind=t-istart
        en_tot[tind]=get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,-1)

    grs=evals[0,1] 
    frs=evals[0,0] 
    gamlin=np.empty(ntime)
    gamlin[:]=2.0*grs
    gamlmode=np.empty(ntime)
    gamtemp=get_energy_single_k(revecs[:,0],revecs[:,0],kx,ky,kz,6)
    gamtemp=gamtemp+get_energy_single_k(revecs[:,0],revecs[:,0],kx,ky,kz,1)
    gamtemp=gamtemp+get_energy_single_k(revecs[:,0],revecs[:,0],kx,ky,kz,2)
    gamtemp=gamtemp+get_energy_single_k(revecs[:,0],revecs[:,0],kx,ky,kz,8)
    gamtemp=gamtemp/get_energy_single_k(revecs[:,0],revecs[:,0],kx,ky,kz,-1)
    gamlmode[:]=gamtemp

    if show_plots:
        plt.title("Energy")
        plt.plot(time[istart:iend+1],en_tot,color='black',label='Total Energy')
        plt.legend(loc='lower right')
        plt.show()

    #Get Q+C stuff
    QC_tot=np.empty(ntime)
    QC_res=np.empty(ntime)
    for t in range(istart,iend+1):
        tind=t-istart
        QC_tot[tind]=get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,6)
        QC_tot[tind]=QC_tot[tind]+get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,1)
        QC_tot[tind]=QC_tot[tind]+get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,2)
        QC_tot[tind]=QC_tot[tind]+get_energy_single_k(g_k[:,tind],g_k[:,tind],kx,ky,kz,8)

    QC_res=QC_tot-gamlmode*en_tot

    if show_plots:
        plt.title("Q+C")
        plt.plot(time[istart:iend+1],QC_tot,color='black',label='Total Q+C')
        plt.plot(time[istart:iend+1],gamlmode*en_tot,label='Q+C one')
        plt.plot(time[istart:iend+1],QC_res,label='Q+C Res.')
        plt.legend(loc='upper left')
        plt.show()

    if show_plots:
        plt.title("Gammas")
        plt.plot(time[istart:iend+1],QC_tot/en_tot,color='black',label='Total Gam')
        plt.plot(time[istart:iend+1],gamlin,marker='+',linestyle='-.',label='Gam lin.')
        plt.plot(time[istart:iend+1],gamlmode,marker='+',linestyle='-.',label='Gam From 1')
        plt.legend(loc='lower right')
        plt.show()

    QC_tota=np.sum(QC_tot)/float(ntime)
    QC_1=np.sum(gamlmode*en_tot)/float(ntime)
    QC_resa=np.sum(QC_res)/float(ntime)
    gam_nla=np.sum(QC_tot/en_tot)/float(ntime)

    return QC_tota,QC_1,QC_resa,gam_nla,gamlin[0],gamlmode[0]




