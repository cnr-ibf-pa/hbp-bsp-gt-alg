from django.shortcuts import render
from django.http import HttpResponse
from django.http import QueryDict

import sys
import json
import csv
import pprint
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.io as io
import random
import logging
import datetime

from scipy.stats import binom
from scipy.special import beta

# set logging up
logging.basicConfig(stream=sys.stdout)
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# create logger if not in DEBUG mode
accesslogger = logging.getLogger('gtalg_access.log')
accesslogger.addHandler(logging.FileHandler('/var/log/gtalg/gtalg_access.log'))
accesslogger.setLevel(logging.DEBUG)

def index(request):
    return render(request, 'model/model.html')

def log_user(request, user_id, pagename):

    accesslogger.info("At " + str(datetime.datetime.now()) + " user " + user_id + " entered page " + pagename)

    return HttpResponse(json.dumps({"response":"ok"}), content_type="application/json")

def read_csv(request,field):
    try:
        stp = str(request.FILES[field].read()).replace("'", "").replace("b", "")
        ret = [list(map(int,row)) for row in csv.reader(stp.split("\\n"), delimiter=',') if row != []]
        return np.array(ret).flatten('C')
    except:
        return []

def ba_dir2(n, m, a, gamma, rho,file,offset):
    T=(np.random.rand(m,m)<rho)*(np.ones((m,m))-np.eye(m))
    for i in range(0,m):
        for j in range(0,m):
            #G[i][0:m]=T[i][0:m]
            if T[i][j]==1:
                file.write(str(offset+i)+"\t"+str(offset+j)+"\r\n")
    targets=np.arange(0,m)
    dout=np.zeros(n)
    din=np.zeros(n)
    dout[0:m]=np.sum(T,axis=1)
    din[0:m]=np.sum(T,axis=0)
    source=int(m)
    D = np.cumsum(gamma)
    while source<n:
        cc=0
        for i in targets:
            #G[int(i)][source]=1
            file.write(str(offset+i)+"\t"+str(offset+source)+"\r\n")
            dout[i]=dout[i]+1
            cc=cc+1
        dout[source]=0
        din[source]=cc
        tr=dout[0:source]+a
        tstr=np.sum(tr)
        tr=tr/tstr
        bins=np.zeros(source+1)
        bins[1:source+1]=np.cumsum(tr)
        mt=source
        while mt>source-1 or mt<1:
            r=random.random()
            mt=min(np.argwhere(r<=D))[0]-1
        j=0
        targets=np.zeros(mt,dtype=int)
        while j<mt:
            r=random.random()
            cb=np.searchsorted(bins,r)-1
            if max(targets==cb)==0:
                targets[j]=cb
                j=j+1
        source=source+1
        if source % 100 == 0:
            pass
            #pprint.pprint('source')
            #pprint.pprint(source)

def impl_inf_deconv(f,h):
    l=min(len(f),len(h))
    A=np.zeros((l,l))
    A[0,0]=f[0]
    for eq in range(1,l):
        for i in range(0,eq+1):
            A[eq,i]=f[eq-i]
    g=np.linalg.solve(A,h)
    return g
	
def test_sim(data_indeg,data_outdeg,Ek,l,rho,m0,pup,pdw,N,tsp):
    binvi=list(range(0,max(data_indeg),5))
    binvo=list(range(0,max(data_outdeg),5))
    #pprint.pprint(data_indeg.shape)
    fig2 = plt.figure(num=2)
    hdo, binsdo, patchesdo = plt.hist(data_outdeg, binvi, density=True, facecolor='b', alpha=0.5, label='data')
    fig1 = plt.figure(num=1)
    hdi, binsdi, patchesdi = plt.hist(data_indeg, binvo, density=True, facecolor='b', alpha=0.5, label='data')
    #plt.show()
	
    ed=binsdi
    data_vedi=[0]
    data_vhi=[0]
    for i in range(0,len(ed)-1):
        if hdi[i]>0:
            data_vedi.append(int((ed[i]+ed[i+1])*0.5))
            data_vhi.append(hdi[i])
    
    ed=binsdo
    data_vedo=[0]
    data_vho=[0]
    for i in range(0,len(ed)-1):
        if hdo[i]>0:
            data_vedo.append(int((ed[i]+ed[i+1])*0.5))
            data_vho.append(hdo[i])
			
    m=np.mean(data_indeg)
	
    n=np.arange(1,np.floor(data_vedo[-1]))
    nz=np.arange(0,np.floor(data_vedo[-1]))
    ni=np.arange(1,np.floor(data_vedi[-1]))
    binsi=np.arange(-0.5,max(data_indeg)+0.5)
    binso=np.arange(-0.5,max(data_outdeg)+0.5)

    N1=int(N/2)	
    N2=int(N/2) 
    M=int(np.floor(N1/l))
    cutoff=0.001
    nf=np.floor(data_vedi[-1])
    hi2=np.interp(np.arange(0,nf+1),data_vedi,data_vhi)
    T=hi2*np.arange(0,nf+1)
	
    Es=np.sum(T)
    E=np.cumsum(T)
    s=np.sum(hi2)
    S=np.cumsum(hi2)
    V=((Es-E)-np.arange(1,nf+2)*(s-S))/(s-S)
    dl=min(np.argwhere(V<(m-Ek)))[0]-1
    c=(m-Ek-m0/N1*(m0-1)*rho)*N1/(N1-m0)
    a=300*c
	
	
    f=beta(nz+c,3)/(beta(c,2));
    ho1=np.interp(n,data_vedo,data_vho)
    hoz1=np.interp(nz,data_vedo,data_vho)
    go1=impl_inf_deconv(f,hoz1)
    p=Ek/N1
	
    hi1=np.zeros(int(np.floor(data_vedi[-1])))
    hi1[0:data_vedi[-1]-1]=np.interp(ni,data_vedi,data_vhi)
	
    #pprint.pprint(str(ni[-1]))
    alpha_true=np.zeros(int(ni[-1]))
    alpha_true[0:int(ni[-1])-dl]=hi1[dl:int(ni[-1])]
    alpha_true=alpha_true/np.sum(alpha_true)
    sa=sum(alpha_true*np.arange(0,int(ni[-1])))
    la=len(alpha_true)

    iota=binom.pmf(np.arange(0,la),m0-1,rho)
    gamma_true=N1/(N1-m0)*alpha_true-m0/(N1-m0)*iota
    gamma_true=gamma_true/sum(gamma_true)
    
    clc=open('/mnt/bsp-gtalg-storage/connectivity/Connectivity'+tsp+'.txt', 'w').close()
    file=open('/mnt/bsp-gtalg-storage/connectivity/Connectivity'+tsp+'.txt', 'w+')
    file.write('This file was created by the Human Brain Project Platform (From the paper Giacopelli et al 2020)\r\n')
    file.write('This is the connectivity data for a network with parameters N='+str(N)+', Ek='+str(Ek)+', l='+str(l)+', rho='+str(rho)+', m0='+str(m0)+', phi_up='+str(pup)+' and phi_down='+str(pdw)+'\r\n' )
    file.write("start\tend\r\n")
	
    ba_dir2(N1, int(m0), a, gamma_true, rho, file, 0)
	
    ba_dir2(N1, int(m0), a, gamma_true, rho, file, N1)
	
    upper_graph_in=(np.random.rand(M,M)<p)
    upper_graph_out=(np.random.rand(M,M)<p)
    part1=np.random.permutation(N1)
    part2=np.random.permutation(N2)
	
    for li in range(0,M):
        for lj in range(0,M):
            if upper_graph_in[li][lj]==1:
                for i in part1[(li-1)*l+1:li*l]:
                    for j in part2[(lj-1)*l+1:lj*l]:
                        if (np.random.rand()<pup):
                            file.write(str(N1+i)+"\t"+str(j)+"\r\n")
            else:
                for i in part1[(li-1)*l+1:li*l]:
                    for j in part2[(lj-1)*l+1:lj*l]:
                        if (np.random.rand()<pdw):
                            file.write(str(N1+i)+"\t"+str(j)+"\r\n")
							
        #pprint.pprint((1.0*li)/M)
						
    for li in range(0,M):
        for lj in range(0,M):
            if upper_graph_out[li][lj]==1:
                for i in part1[(li-1)*l+1:li*l]:
                    for j in part2[(lj-1)*l+1:lj*l]:
                        if (np.random.rand()<pup):
                            file.write(str(i)+"\t"+str(N1+j)+"\r\n")
            else:
                for i in part1[(li-1)*l+1:li*l]:
                    for j in part2[(lj-1)*l+1:lj*l]:
                        if (np.random.rand()<pdw):
                            file.write(str(i)+"\t"+str(N1+j)+"\r\n")
							
    #pprint.pprint((1.0*li)/M)
    file.close()
    #pprint.pprint('File saved!')
	
    

def actionpage(request):
    #pprint.pprint(request.POST)
	
    try:
        N=int(request.POST['N'])
        Ek=int(request.POST['Ek'])
        l=int(request.POST['l'])
        rho=float(request.POST['rho'])
        m0=int(request.POST['m0'])
        pup=float(request.POST['phiu'])
        pdw=float(request.POST['phid'])
        step=float(request.POST['step'])
    except:
        pprint.pprint('Error in loading parameters!')
        return render(request, 'model/error2.html')
	
    flag_data=0
	
    indeg=read_csv(request,'in_file')
    outdeg=read_csv(request,'out_file')

	
    if (len(indeg)+len(outdeg))==0:
        Conn=read_csv(request,'conn_file')
        lC=len(Conn)
        if lC>0:
            Conn.shape=(int(lC/2),2)
            #pprint.pprint(Conn)
            OutConn=np.array(Conn[:,0])
            InConn=np.array(Conn[:,1])
	    #pprint.pprint(InConn)
            UNQ=np.concatenate((OutConn,InConn))
            GIDs=np.unique(UNQ)
            indeg=[]
            outdeg=[]
            #pprint.pprint('Converting!')
            cg=0
            for gid in GIDs:
                indeg.append(list(InConn).count(gid))
                outdeg.append(list(OutConn).count(gid))
                if cg%1000==0:
                    pass
                    #pprint.pprint(cg)
                cg=cg+1
            indeg=np.array(indeg)
            outdeg=np.array(outdeg)
	    #pprint.pprint('Data loaded!')
		
	
	
    if (len(indeg)+len(outdeg))>0:
        pass
	#pprint.pprint('Data loaded!')
    else:
        indeg=np.genfromtxt('static/dinII.csv',delimiter=',').astype(int)
        outdeg=np.genfromtxt('static/doutII.csv',delimiter=',').astype(int)
	
    t = time.time()
    if N<1:
        N=len(indeg)
        N=int(np.floor(N/(2*l))*2*l)
    tsp=str(int(10**4*time.time()))
    try:
        test_sim(indeg,outdeg,Ek,l,rho,m0,pup,pdw,N,tsp)
        sim_indeg=np.zeros(N)
        sim_outdeg=np.zeros(N)
        
        lineList = [line.rstrip('\n') for line in open('/mnt/bsp-gtalg-storage/connectivity/Connectivity'+tsp+'.txt')]

        lln=len(lineList)
        #pprint.pprint(lineList)
        for ll in range(10,lln):
            tl=lineList[ll]
            if len(tl)>0:
		#pprint.pprint(tl)
                tp=tl.split('\t')
                inidx=int(tp[1])
                outidx=int(tp[0])
                #pprint.pprint(str(outidx)+' '+str(inidx))
                sim_indeg[inidx]=sim_indeg[inidx]+1
                sim_outdeg[outidx]=sim_outdeg[outidx]+1
        elapsed = time.time() - t
        #pprint.pprint(sim_indeg)
        #pprint.pprint('Elapsed: '+str(elapsed))
        
        instart=0
        outstart=0
        inend=max(indeg)
        outend=max(outdeg)
        instep=step
        outstep=step
        instat=[instart,inend,instep]
        outstat=[outstart,outend,outstep]
        binsi=np.arange(instart,inend,instep)
        binso=np.arange(outstart,outend,outstep)
        plotin=np.append(np.array([0]), np.array((binsi[0:-1]+binsi[1:])*0.5))
        plotout=np.append(np.array([0]), np.array((binso[0:-1]+binso[1:])*0.5))
        temp1=np.histogram(indeg, bins=binsi, density=True)
        hdi=np.append(np.array([0]), temp1[0])
        temp2=np.histogram(outdeg, bins=binso, density=True)
        hdo=np.append(np.array([0]), temp2[0])
        temp3=np.histogram(sim_indeg, bins=binsi, density=True)
        hsi=np.append(np.array([0]), temp3[0])
        temp4=np.histogram(sim_outdeg, bins=binso, density=True)
        hso=np.append(np.array([0]), temp4[0])
        dp=dict(request.POST)
        st='aaa'
        dv={'indeg':indeg.tolist(),'outdeg':outdeg.tolist(),'sim_indeg':sim_indeg.tolist(),'sim_outdeg':sim_outdeg.tolist(),'instat':instat,'outstat':outstat,'timestamp':tsp,'filename':'Connectivity'+tsp+'.txt','plotin':plotin.tolist(),'plotout':plotout.tolist(),'hdi':hdi.tolist(),'hdo':hdo.tolist(),'hsi':hsi.tolist(),'hso':hso.tolist()}
        df={**dp,**dv}
        qdict = QueryDict('', mutable=True)
        qdict.update(df)
        context = qdict
        return render(request, 'model/model2.html', context)
    except:
        return render(request, 'model/error.html')
