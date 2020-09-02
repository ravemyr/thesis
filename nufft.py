# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 19:17:31 2020

@author: Emanuel
"""

import numpy as np
import matplotlib.pyplot as plt


def g(x,t):
    l = np.arange(-1000,1000)
    g = np.sum(np.exp(-(x-2*np.pi*l)**2/(4*t)))
    return g

def f(x):
    return np.sin(x) + np.multiply(np.cos(x),x)

def validate(fk,x, M):
    # Add check that fj is np.array
    F = np.empty([1,M],dtype=np.complex128)
    for k in range(int(-M/2),int(M/2)):
        F[0,k+int(M/2)] = np.sum(np.multiply(fk.astype(complex),np.exp(-1j*k*x)))
    return F

N = 16
M = 2*N
Mr = M
t_vec = [2/(Mr**2), 6/(Mr**2), 12/(Mr**2), 24/(Mr**2)]

data = np.random.rand(N)*2*np.pi
data = np.sort(data)
fj = f(data)


val = validate(fj,data,M)
ft = np.empty([1,Mr])
res = np.empty([1,len(t_vec)])
p = 0
for t in t_vec:
    for m in range(0,Mr):
        for n in range(0,N):
            ft[0,m] = ft[0,m] + g(2*np.pi*m/Mr - data[n],t)*fj[n]
          
    ff = np.fft.fft(np.multiply(ft,np.exp(1.j*M*np.pi/Mr)))/Mr
    fff = np.empty([1,Mr],dtype=complex)
    for k in range(int(-M/2),int(M/2)):
        fff[0,k+int(M/2)] = ff[0,k+int(M/2)]*np.exp(t*k**2)*np.sqrt(np.pi/t)
    res[0,p] = np.linalg.norm(fff-val)
    p = p+1
print(res)
plt.plot(t_vec,np.transpose(res))