import numpy as np
import matplotlib.pyplot as plt
from dna_diags import *

def get_time_from_tempfile(ikx,iky):
    cikx=str(int(ikx))
    if len(cikx)==1:
        cikx='0'+cikx
    elif len(cikx) > 2:
        print cikx
        print "Error in get_time_from_temp"
        stop
    
    ciky=str(int(iky))
    if len(ciky)==1:
        ciky='0'+ciky
    elif len(ciky) > 2:
        print "Error in get_time_from_gk"
        stop

    diagdir=par['diagdir'][1:-1]
    print diagdir
        
    file_name='temp_kx'+cikx+'ky'+ciky+'.dat'

    print "Reading file",diagdir+'/'+file_name
    file_exists=os.path.isfile(diagdir+'/'+file_name)
    if file_exists:
        pass
    else:
        print "File does not exist:",diagdir+'/'+file_name
        stop

    f=open(diagdir+'/'+file_name,'r')
    ntot=3*par['nkz0']
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
      else:
          continue_read=0

    f.close()
    return time


def test_gknl(ikx,iky,start_time=-1.0,end_time=-1.0,calc_from_gout=False):

    time=get_time_from_gkfile(ikx,iky,read_nl=True)
    time0=get_time_from_gout()


    if start_time==-1.0:
        start_time=time[0]
    if end_time==-1.0:
        end_time=time[len(time)-1]
    if start_time > end_time:
        stop

    istart=np.argmin(abs(time-start_time))
    iend=np.argmin(abs(time-end_time))
    ntime=iend-istart+1

    istart0=np.argmin(abs(time0-start_time))
    iend0=np.argmin(abs(time0-end_time))
    ntime0=iend0-istart0+1

    if ntime != ntime0:
        print "Error in test_gknl!"
        stop

    kxgrid,kygrid,kzgrid,herm_grid=get_grids()

    #Get from gk file
    etot=np.zeros(ntime,dtype='float')
    #Get from time derivative of above 
    detot=np.zeros(ntime,dtype='float')
    #Get from gkfile
    erhslin=np.zeros(ntime,dtype='float')
    #Get from calculated nonlinearity from g_out
    erhsnl0=np.zeros(ntime,dtype='float')
    #Get from nl file
    erhsnl=np.zeros(ntime,dtype='float')
    #Get from temp file
    terhsnl=np.zeros(ntime,dtype='float')
    #Get from temp file
    terhs=np.zeros(ntime,dtype='float')
    #Get from temp file
    tetot=np.zeros(ntime,dtype='float')
    #Get from calculated nonlinearity from g_out
    erhsnl0=np.zeros(ntime0,dtype='float')

    #Get from gk file
    etemp=np.zeros(par['nkz0'],dtype='float')
    #Get from gk file
    rhstemp=np.zeros(par['nkz0'],dtype='float')
    #Get from nl file
    rhsnltemp=np.zeros(par['nkz0'],dtype='float')
    #Get from gout file
    rhsnltemp0=np.zeros(par['nkz0'],dtype='float')
    #Get from temp file
    tetemp=np.zeros(par['nkz0'],dtype='float')
    trhstemp=np.zeros(par['nkz0'],dtype='float')
    trhsnltemp=np.zeros(par['nkz0'],dtype='float')

    ################
    ################
    #etotk=np.zeros((par['nkz0'],ntime0),dtype='float')
    #erhslink=np.zeros((par['nkz0'],ntime0),dtype='float')
    #erhsnlk=np.zeros((par['nkz0'],ntime0),dtype='float')
    ################
    ################

    for i in range(istart,iend+1):
        it=i-istart
        nl=read_time_step_gkfile(ikx,iky,i,read_nl=True)
        nl=np.reshape(nl,(par['nkz0'],par['nv0']),order='F')
        print it+1, " of ", ntime0
        gt0=read_time_step_gkfile(ikx,iky,i)
        temp=read_time_step_tempfile(ikx,iky,i)
        e_from_temp=temp[:,0]
        gt0=np.reshape(gt0,(par['nkz0'],par['nv0']),order='F')
        #print np.info(gt0)
        g0=read_time_step_g(i)
        g0=np.reshape(g0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
        for k in range(par['nkz0']):
        #for k in range(1,2):
            #print k, " of ", par['nkz0']
            if calc_from_gout:
                rhsnltemp0[k]=get_energy_single_k(g0,g0,kxgrid[ikx],kygrid[iky],kzgrid[k],9)
            rhstemp[k]=get_energy_single_k(gt0[k,:],gt0[k,:],kxgrid[ikx],kygrid[iky],kzgrid[k],0)
            etemp[k]=get_energy_single_k(gt0[k,:],gt0[k,:],kxgrid[ikx],kygrid[iky],kzgrid[k],-1)
            eop=energy_operator_single_k(gt0[k,:],kxgrid[ikx],kygrid[iky]) 
            rhsnltemp[k]=np.real(np.sum(eop*nl[k,:]))
            #temp file stuff
            tetemp[k]=temp[k,0] 
            trhstemp[k]=temp[k,1] 
            trhsnltemp[k]=temp[k,2] 
    ################
    ################
    #        erhslink[k,it]=rhstemp[k]
    #        erhsnlk[k,it]=rhsnltemp[k]
    #        etotk[k,it]=etemp[k]
    ################
    ################
        erhslin[it]=np.sum(rhstemp)
        erhsnl0[it]=np.sum(rhsnltemp0)
        erhsnl[it]=np.sum(rhsnltemp)
        etot[it]=np.sum(etemp)
        #temp file stuff
        terhs[it]=np.sum(trhstemp)
        terhsnl[it]=np.sum(trhsnltemp)
        tetot[it]=np.sum(tetemp)
    
    for i in range(ntime-1):
        detot[i]=(etot[i+1]-etot[i])/(time[i+1]-time[i])

    plt.plot(time[istart:iend+1],etot,label='FE (gk)')
    plt.plot(time[istart:iend+1],tetot,'x-',label='Etot temp')
    plt.legend()
    plt.show()

    plt.plot(time[istart:iend+1],erhslin,label='RHS lin (gk)')
    plt.plot(time[istart:iend+1],terhs-terhsnl,'x-',label='RHS lin temp')
    plt.legend()
    plt.show()

    plt.plot(time[istart:iend+1],erhslin-(terhs-terhsnl),label='RHS lin - RHS lin (gk/temp)')
    plt.legend()
    plt.show()

    plt.plot(time[istart:iend+1],erhsnl,label='RHS nl (nl)')
    if calc_from_gout:
        plt.plot(time[istart:iend+1],erhsnl0,'+-',label='RHS nl0 (gout)')
    plt.plot(time[istart:iend+1],terhsnl,'x-',label='RHS nl temp')
    plt.legend()
    plt.show()

    plt.plot(time[istart:iend+1],terhs,'x-',label='RHS temp')
    plt.plot(time[istart:iend+1],terhsnl,'x-',label='RHSnl temp')
    plt.legend()
    plt.show()


    plt.plot(time[istart:iend+1],erhslin+erhsnl,label='RHS tot (gk+nl)')
    #plt.plot(time[istart:iend+1],erhslin+erhsnl0,label='RHS tot0 (gk_gout)')
    plt.plot(time[istart:iend+1],terhs,'x',label='RHS temp')
    plt.plot(time[istart:iend+1],detot,label='detot')
    plt.legend()
    plt.show()

    ################
    ################
    #for k in range(par['nkz0']):
    #    for i in range(istart,iend):
    #        detot[i]=2.0*(etotk[k,i+1]-etotk[k,i])/(time[i+1]-time[i])
    #    plt.plot(time[istart:iend+1],erhslink[k,:]+erhsnlk[k,:],label='RHS tot')
    #    plt.plot(time[istart:iend+1],detot,label='detot')
    #    plt.legend()
    #    plt.show()
    ################
    ################

    #for i in range(istart0,iend0+1):
    #    it=i-istart0
    #    print it+1, " of ", ntime0
    #    gt0=read_time_step_g(i)
    #    gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
    #    for k in range(par['nkz0']):
    #        print k, " of ", par['nkz0']
    #        rhs[:,k]=get_rhs_nl_single_k(gt0,kxgrid[ikx],kygrid[iky],kzgrid[k])
    #        nl_tot0[k,it]=np.real(np.sum(rhs[:,k]))
    #        print nl_tot[k,it],nl_tot0[k,it]
           
      
 
def read_time_step_tempfile(ikx,iky,which_itime):
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
        
    file_name='temp_kx'+cikx+'ky'+ciky+'.dat'
    #print "Reading file",diagdir+'/'+file_name
    file_exists=os.path.isfile(diagdir+'/'+file_name)
    if file_exists:
        pass
    else:
        print "File does not exist:",diagdir+'/'+file_name
        stop


    f = open(diagdir+'/'+file_name,'rb')
    ntot=3*par['nkz0']
    mem_tot=ntot*8
    gt0=np.empty((par['nkz0'],par['nv0']))
    f.seek(8+which_itime*(8+mem_tot))
    gt0=np.fromfile(f,dtype='float64',count=ntot)
    gt0=np.reshape(gt0,(par['nkz0'],3),order='F')
    #print sum(gt0)
    f.close()
    return gt0


