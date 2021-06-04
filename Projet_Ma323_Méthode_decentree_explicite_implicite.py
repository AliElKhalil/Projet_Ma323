# -*- coding: utf-8 -*-
"""
Created on Mon May 24 17:40:52 2021

@author: elena
"""

import numpy as np
import matplotlib.pyplot as plt

#Construction matrice explicite décentré amont
def matrice_eda(N,c,d):
    A = (1-d-2*c)*np.eye(N)+(d+c)*np.eye(N,N,-1)+c*np.eye(N,N,1)
    return(A)

#Construction matrice implicite décentré amont
def matrice_ida(N,c,d): 
    A=(1+d+2*c)*np.eye(N)+(-d-c)*np.eye(N,N,-1)-c*np.eye(N,N,1)
    return(A)

def construitX(N):
    X = np.linspace(-10,10,N+1)
    return(X)

def U0(x):
    U = np.exp(-x**2)
    return(U)
def U02(x):
    U=np.exp(-(x-2)**2)+np.exp(-(x+2)**2)
    return(U)
#résolution
def SolDA(h,N,tau,c,d):
    ntfinal=int(Tmax/tau)+1
    ADA = matrice_eda(N,c,d)
    UT=np.zeros((ntfinal,N))
    T=np.arange(ntfinal)*tau
    X=construitX(N)
    a=int(input("Voulez-vous réaliser les test pour U0 question8 ou U0 question9 ?"))
    if (a==8):
        U=U0(X[0:N])
    else:
        U=U02(X[0:N])
    UT[0,0:N]=U
    for n in range(ntfinal-1) :
        U=ADA@U
        UT[n+1,0:N]=U
    return X,T,UT

def SolDAimplicite(h,N,tau,c,d):
    ntfinal=int(Tmax/tau)+1
    ADAimp = matrice_ida(N,c,d)
    print(ADAimp)
    print(h)
    print(N)
    print(c)
    print(d)
    UT=np.zeros((ntfinal,N))
    T=np.arange(ntfinal)*tau
    X=construitX(N)
    U=U0(X[0:N])
    UT[0,0:N]=U
    for n in range(ntfinal-1) :
        U=np.linalg.solve(ADAimp,U)
        UT[n+1,0:N]=U
    return X,T,UT

def UnTest(h,N,c,d,tau,Tmax,solveur,nomsolv,ax,listet): 
    X,T,U=solveur(h,N,tau,c,d)
    ntfinal = int(Tmax/tau)+1
    b = np.zeros((ntfinal,1))
    b[:,0] = U[:,0]
    U = np.concatenate((U,b),axis=1)
    for t in listet :
        n=int(t/tau)
        ax.plot(X,U[n,:],label='t='+str(n*tau))
    ax.legend()
    ax.set(title=nomsolv+' h={:.3f}\n tau={:.4f}'.format(20/(N+1),tau))

Tmax = 4. #Initialisation de Tmax
liste_t=[0.,1.,2.,3.,4.] #Initialisation du vecteur t avec 0<=t<=4
Solveurs=[SolDA,SolDAimplicite] #Initialisation du veteur solveurs reprenant toutes les méthodes que nous allons utiliser
Nomssolveurs=['Méthode explicite décentrée amont','Méthode implicite décentrée amont'] #Initialisation du vecetur reprenan toute sles chaines de caractères de chaque méthode

#Définition des différents cas à traiter, [h,tau]
Cas=[[0.2,0.5],[0.1,0.05],[0.05,0.0025]]

#définition de la fenêtre d'affichage des différents graphes que nous allons représenter
fig,liste_ax=plt.subplots(len(Cas),len(Solveurs),sharex='all',figsize=(25,25))

#boucle de la taille du nombre de cas traités permettant de traiter tous les cas avec toutes les méthodes
for cm in range(len(Cas)) : 
    m=Cas[cm]
    print("\n\nCas :", m ,"\n")
    c = (1/20)*m[1]/(m[0]**2)
    d=m[1]/m[0]
    N = int(20/(m[0])-1)
    for k in range(len(Solveurs)):
        print('----------------------------------------------------\n','h=',m[0],'tau=',m[1],Nomssolveurs[k])
        UnTest(m[0],N,c,d,m[1],Tmax,Solveurs[k],Nomssolveurs[k],liste_ax[cm][k],liste_t)
    print('----------------------------------------------------\n')        

plt.show()
