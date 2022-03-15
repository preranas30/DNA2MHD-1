import numpy as np
import matplotlib.pyplot as plt
from config import *
from dna_diags import *

def nlt_test_old(time,start_time=-1.0,end_time=-1.0):
    kx_nlt,ky_nlt,kz_nlt=get_grids_nlt()
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    print "kxgrid_nlt",kx_nlt
    print "kygrid_nlt",ky_nlt
    print "kzgrid_nlt",kz_nlt


    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time >= end_time:
        stop

    nkx1=(2*par['nkx0']-1)
    nky1=(par['nky0']-1)
    nkz1=(par['nkz0']-1)
    NLT=np.zeros((nkx1,nky1,nkz1))
    #NLT_herm=np.zeros(par['nv0'])
    #rhs=np.zeros(par['nv0'])
 
    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1
    print "ntime",ntime
    nlt_tot=np.zeros(ntime)
    rhs_tot=np.zeros(ntime)
    en_tot=np.zeros(ntime)
    dedt=np.zeros(ntime)

    ind_in=raw_input("Enter k indices:")
    ind_split=ind_in.split()
    ikx=int(float(ind_split[0]))
    iky=int(float(ind_split[1]))
    ikz=int(float(ind_split[2]))

    print "indices:",ikx,iky,ikz    
    print "k's:",kxgrid[ikx],kygrid[iky],kzgrid[ikz]    

    kxstring=str(kxgrid[ikx])[0:6]
    kystring=str(kygrid[iky])[0:6]
    kzstring=str(kzgrid[ikz])[0:6]
    print "kxstring",kxstring
    print "kystring",kystring
    print "kzstring",kzstring

    for i in range(istart,iend+1):
        print 'time=',time[i],' of ',time[iend]
        gt0=read_time_step_g(i)
        gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        NLT0=get_nlt_singlek(gt0,ikx,iky,ikz)
        NLT=NLT+NLT0
        nlt_tot[i-istart]=np.sum(np.sum(np.sum(NLT0,axis=0),axis=0),axis=0)
        rhs0=get_energy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],\
                                 kxgrid[ikx],kygrid[iky],kzgrid[ikz],0)
        rhs_tot[i-istart]=rhs0
        en0=get_energy_single_k(gt0[ikx,iky,ikz,:],gt0[ikx,iky,ikz,:],\
                                 kxgrid[ikx],kygrid[iky],kzgrid[ikz],-1)
        en_tot[i-istart]=en0
        #rhs=rhs+rhs0

    NLT=NLT/(float(ntime))
    #rhs=rhs/(float(ntime))
    NLT_kzavg=np.sum(NLT,axis=2) 
    plt.contourf(kx_nlt,ky_nlt,np.transpose(NLT_kzavg),200)
    plt.colorbar()
    plt.xlabel(r'$k_x^\prime \rho_s$',size=18)
    plt.ylabel(r'$k_y^\prime \rho_s$',size=18)
    title_str=(r'$T_{k,k^\prime}(k_x \rho_s =$'+kxstring+r'$,k_y \rho_s =$'+ kystring+'$,k_z R=$'+kzstring+')')
    plt.title(title_str)
    plt.show()

    #print "sum NLT",np.sum(np.sum(np.sum(NLT,axis=0),axis=0),axis=0)
    #print "sum rhs",np.real(np.sum(rhs))
    for i in range(0,ntime-1):
        dedt[i]=(en_tot[i+1]-en_tot[i])/(time[i+istart+1]-time[i+istart])
    print "nlt_tot",nlt_tot
    print "rhs_tot",rhs_tot
    plt.plot(time[istart:iend+1],nlt_tot,label='nlt')
    plt.plot(time[istart:iend+1],rhs_tot,'-x',label='rhs')
    plt.plot(time[istart:iend+1],en_tot,'-*',label='energy')
    plt.plot(time[istart:iend+1],rhs_tot+nlt_tot,'-+',label='rhs+nlt')
    plt.plot(time[istart+1:iend+1],dedt[1:],'-o',label='dedt')
    plt.legend()
    plt.show()
    plt.plot(time[istart:iend+1],rhs_tot+nlt_tot,'-',label='rhs+nlt')
    plt.plot(time[istart+1:iend+1],dedt[1:],'-',label='dedt')
    plt.legend()
    plt.show()

def get_nlt_singlek(g_in,ikx,iky,ikz):
    phib=get_phi_bar(g_in)
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    kxgrid_nlt,kygrid_nlt,kzgrid_nlt=get_grids_nlt()
    #print "kxgrid_nlt",kxgrid_nlt
    #print "kygrid_nlt",kygrid_nlt
    #print "kzgrid_nlt",kzgrid_nlt
    nkx1=(2*par['nkx0']-1)
    nky1=(par['nky0']-1)
    nkz1=(par['nkz0']-1)
    NLT=np.zeros((nkx1,nky1,nkz1))
    kx=kxgrid[ikx]
    ky=kygrid[iky]
    kz=kzgrid[ikz]
    phibk=phib[ikx,iky,ikz]
    gk=g_in[ikx,iky,ikz,:]
    #The following loops over k prime
    for i in range(nkx1):
        print i, " of ", nkx1-1
        for j in range(nky1):
            for k in range(nkz1):
                kxp=kxgrid_nlt[i]
                kyp=kygrid_nlt[j]
                kzp=kzgrid_nlt[k]
                kxpp=kx-kxp
                kypp=ky-kyp
                kzpp=kz-kzp
                if np.abs(kxp) < par['kxmax0'] and \
                   np.abs(kyp) < par['kymax0'] and \
                   np.abs(kzp) < par['kzmax0'] and \
                   np.abs(kxpp) < par['kxmax0'] and \
                   np.abs(kypp) < par['kymax0'] and \
                   np.abs(kzpp) < par['kzmax0']:

                    ikxp,ikyp,ikzp,take_conjg_p=get_all_indices(kxp,kyp,kzp)
                    #print "kxp,kyp,kzp",kxp,kyp,kzp
                    #print "ikxp,ikyp,ikzp,take_conjg",ikxp,ikyp,ikzp,take_conjg_p
                    ikxpp,ikypp,ikzpp,take_conjg_pp=get_all_indices(kxpp,kypp,kzpp)
                    #print "kxpp,kypp,kzpp",kxpp,kypp,kzpp
                    #print "ikxpp,ikypp,ikzpp,take_conjg",ikxpp,ikypp,ikzpp,take_conjg_pp
                    ckkp=kxp*ky-kx*kyp
                    if take_conjg_p:
                        phib_p=np.conj(phib[ikxp,ikyp,ikzp])
                        g_p=np.conj(g_in[ikxp,ikyp,ikzp,:])
                    else:
                        phib_p=(phib[ikxp,ikyp,ikzp])
                        g_p=(g_in[ikxp,ikyp,ikzp,:])
                    if take_conjg_pp:
                        phib_pp=np.conj(phib[ikxpp,ikypp,ikzpp])
                        g_pp=np.conj(g_in[ikxpp,ikypp,ikzpp,:])
                    else:
                        phib_pp=(phib[ikxpp,ikypp,ikzpp])
                        g_pp=(g_in[ikxpp,ikypp,ikzpp,:])


                    NLT[i,j,k]=np.real(np.pi**0.25*ckkp*np.conj(phibk)*\
                        phib_p*g_pp[0])
                    for n in range(par['nv0']):
                        NLT[i,j,k]=NLT[i,j,k]-np.real(np.pi**0.5*ckkp*np.conj(gk[n])*\
                            phib_pp*g_p[n])
                    #NLT[i,j,k]=np.pi**0.25*ckkp*np.conj(phibk)*\
                    #    phib[ikxp,ikyp,ikzp]*g_in[ikxpp,ikypp,ikzpp,0]
                    #for n in range(par['nv0']):
                    #    NLT[i,j,k]=NLT[i,j,k]-np.pi**0.5*ckkp*np.conj(gk[n])*\
                    #        phib[ikxpp,ikypp,ikzpp]*g_in[ikxp,ikyp,ikzp,n]
    return NLT

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

def get_grids_nlt():
    nkx1=(2*par['nkx0']-1)
    nky1=(par['nky0']-1)
    nkz1=(par['nkz0']-1)
    kxgrid=np.arange(nkx1)/(float(nkx1-1))*par['kxmax0']*2-par['kxmax0']
    kygrid=np.arange(nky1)/(float(nky1-1))*par['kymax0']*2-par['kymax0']
    kzgrid=np.arange(nkz1)/(float(nkz1-1))*par['kzmax0']*2-par['kzmax0']
    #print par['kxmax']
    #print par['kymax']
    #print par['kzmax']
    #print kxgrid
    #print kygrid
    #print kzgrid
    return kxgrid,kygrid,kzgrid

def get_index_from_kx_nlt(kx):
    nkx=kx+par['kxmax0']
    ikx=np.rint(nkx/par['kxmin'])
    return ikx

def get_index_from_ky_nlt(ky):
    nky=ky+par['kymax0']
    iky=np.rint(nky/par['kymin'])
    return iky

def get_index_from_kz_nlt(kz):
    nkz=kz+par['kzmax0']
    ikz=np.rint(nkz/par['kzmin'])
    return ikz

def grid_tests():
    kxgrid,kygrid,kzgrid,herm_grid=get_grids()
    kxgrid_nlt,kygrid_nlt,kzgrid_nlt=get_grids_nlt()
    print "kx test"
    for i in range(len(kxgrid)):
        ikx=get_index_from_kx(kxgrid[i])
        if i-ikx != 0.0:
            print "Error!", i,ikx
    print "ky test"
    for i in range(len(kygrid)):
        iky=get_index_from_ky(kygrid[i])
        if i-iky != 0.0:
            print "Error!", i,iky
    print "kz test"
    for i in range(len(kzgrid)):
        ikz=get_index_from_kz(kzgrid[i])
        if i-ikz != 0.0:
            print "Error!", i,ikz
    print "kx nlt test"
    print kxgrid_nlt
    for i in range(len(kxgrid_nlt)):
        ikx=get_index_from_kx_nlt(kxgrid_nlt[i])
        if i-ikx != 0.0:
            print "Error!", i,ikx
    print "ky nlt test"
    print kygrid_nlt
    for i in range(len(kygrid_nlt)):
        iky=get_index_from_ky_nlt(kygrid_nlt[i])
        if i-iky != 0.0:
            print "Error!", i,iky
    print "kz nlt test"
    print kzgrid_nlt
    for i in range(len(kzgrid_nlt)):
        ikz=get_index_from_kz_nlt(kzgrid_nlt[i])
        if i-ikz != 0.0:
            print "Error!", i,ikz


