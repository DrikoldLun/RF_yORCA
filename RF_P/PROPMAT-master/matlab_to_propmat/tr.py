#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:22:04 2020

@author: lun
"""

import numpy as np
import sys
from numpy import sqrt, power

def eta(c,p):
    return sqrt(complex(power(1/c,2)-power(p,2)))

# from top to bottom (for P-SV system)
def tr_coef(a1,b1,a2,b2,rho1,rho2,p):
    etaa1, etaa2, etab1, etab2 = eta(a1,p), eta(a2,p), eta(b1,p), eta(b2,p)
    a = rho2*(1-2*power(b2*p,2))-rho1*(1-2*power(b1*p,2))
    b = rho2*(1-2*power(b2*p,2))+2*rho1*power(b1*p,2)
    c = rho1*(1-2*power(b1*p,2))+2*rho2*power(b2*p,2)
    d = 2*(rho2*power(b2,2)-rho1*power(b1,2))
    E = b*etaa1+c*etaa2
    F = b*etab1+c*etab2
    G = a-d*etaa1*etab2
    H = a-d*etaa2*etab1
    D = E*F+G*H*power(p,2)
    A = power(power(1/b1,2)-2*power(p,2),2)+4*power(p,2)*etaa1*etab1
    coef = {}
    
    # incident P
    coef['Rpp'] = ((b*etaa1-c*etaa2)*F-(a+d*etaa1*etab2)*H*power(p,2))/D
    coef['Rpsv'] = -(2*etaa1*(a*b+c*d*etaa2*etab2)*p*(a1/b1))/D
    coef['Tpp'] = (2*rho1*etaa1*F*(a1/a2))/D
    coef['Tpsv'] = (2*rho1*etaa1*H*p*(a1/b2))/D
    
    # incident SV
    coef['Rsvsv'] = -((b*etab1-c*etab2)*E-(a+d*etaa2*etab1)*G*power(p,2))/D
    coef['Rsvp'] = -(2*etab1*(a*b+c*d*etaa2*etab2)*p*(b1/a1))/D
    coef['Tsvsv'] = (2*rho1*etab1*E*b1/b2)/D
    coef['Tsvp'] = -(2*rho1*etab1*G*p*b1/a2)/D

    # free surface only related to medium 1
    coef['Rpp_f'] = (-power((power(1/b1,2)-2*power(p,2)),2)+4*power(p,2)*etaa1*etab1)/A
    coef['Rpsv_f'] = (4*(a1/b1)*p*etaa1*(power(1/b1,2)-2*power(p,2)))/A
    coef['Rsvp_f'] = (4*(b1/a1)*p*etab1*(power(1/b1,2)-2*power(p,2)))/A
    
    return coef

def t_layer(c,h,rayp):
    return h/(c*sqrt(1-power(rayp*c,2)))

# Ps conv_layer_index from top to the bottom (from 0), this index <= n_layer-2
def Ps_t_amp(a_l,b_l,h_l,rho_l,inci,conv_layer_index):
    conv_layer_index = int(conv_layer_index)
    p = np.sin(inci)/a_l[-1] # inci is defined at the bottom interface
    inci_fa = np.arcsin(a_l[0]*p)
    inci_fb = np.arcsin(b_l[0]*p)
    t_P = 0
    amp_P = 1
    t_Ps = 0
    amp_Ps = 1
    n_l = len(a_l) # number of layers
    if conv_layer_index > n_l-2:
        print('conversion index is too large!')
        sys.exit(1)
    result = {}
    for i in range(n_l):
        if i < n_l-1:
            coef = tr_coef(a_l[i+1],b_l[i+1],a_l[i],b_l[i],rho_l[i+1],rho_l[i],p)
            #print(coef['Tsvsv'],coef['Tpp'],coef['Tpsv'])
            if i <= conv_layer_index:
                t_Ps += t_layer(b_l[i],h_l[i],p)
                if i < conv_layer_index:
                    amp_Ps *= coef['Tsvsv']
                else:
                    amp_Ps *= coef['Tpsv']
            else:
                t_Ps += t_layer(a_l[i],h_l[i],p)
            amp_P *= coef['Tpp']
        else:
            t_Ps += t_layer(a_l[-1],h_l[-1],p)
        t_P += t_layer(a_l[i],h_l[i],p)
    
    coef_surf = tr_coef(a_l[0],b_l[0],a_l[0],b_l[0],rho_l[0],rho_l[0],p)
    coef_conv = tr_coef(a_l[conv_layer_index],b_l[conv_layer_index],a_l[conv_layer_index+1],b_l[conv_layer_index+1]\
                        ,rho_l[conv_layer_index],rho_l[conv_layer_index+1],p)
    # multiple reflection
    amp_PpPs = amp_P*coef_surf['Rpp_f']*coef_conv['Rpsv']
    amp_PsPs = amp_Ps*coef_surf['Rsvp_f']*coef_conv['Rpsv']
    amp_PpSs = amp_P*coef_surf['Rpsv_f']*coef_conv['Rsvsv']
    t_PpPs = t_P
    t_PsPs = t_Ps
    t_PpSs = t_P
    # multiple downward
    for i in range(conv_layer_index+1):
        coef = tr_coef(a_l[i],b_l[i],a_l[i+1],b_l[i+1],rho_l[i],rho_l[i+1],p)
        amp_PpPs *= coef['Tpp']
        amp_PsPs *= coef['Tpp']
        amp_PpSs *= coef['Tsvsv']
        t_PpPs += t_layer(a_l[i],h_l[i],p)
        t_PsPs += t_layer(a_l[i],h_l[i],p)
        t_PpSs += t_layer(b_l[i],h_l[i],p)
    # multiple upward
    for i in range(conv_layer_index+1):
        coef = tr_coef(a_l[i+1],b_l[i+1],a_l[i],b_l[i],rho_l[i+1],rho_l[i],p)
        amp_PpPs *= coef['Tsvsv']
        amp_PsPs *= coef['Tsvsv']
        amp_PpSs *= coef['Tsvsv']
        t_PpPs += t_layer(b_l[i],h_l[i],p)
        t_PsPs += t_layer(b_l[i],h_l[i],p)
        t_PpSs += t_layer(b_l[i],h_l[i],p)
    
    amp_2p2s = amp_PsPs+amp_PpSs

    # Z,R component
    Z_amp = np.hstack((np.array([amp_P*np.cos(inci_fa)]),np.array([amp_Ps,amp_PpPs,amp_2p2s])*np.sin(inci_fb)))
    R_amp = np.hstack((np.array([amp_P*np.sin(inci_fa)]),np.array([amp_Ps,amp_PpPs,amp_2p2s])*np.cos(inci_fb)))
    t = [t_P,t_Ps,t_PpPs,t_PpSs]
    phase = ['P','Ps','PpPs','PpSs+PsPs']
    
    for key, value in zip(['t','Zamp','Ramp','phase'],[t,Z_amp.tolist(),R_amp.tolist(),phase]):
        result[key] = value
        
    return result

'''
# Zhu & Kanamori 2000
inci = np.deg2rad(25)
Vp_l = [5.6744,6.3008,6.3528,7.7985]
Vs_l = [3.28,3.64,3.67,4.5]
rho_l = [2.6238,2.7754,2.7894,3.3223]
h_l = [5.5,10.5,16,68]
'''
inci = np.deg2rad(15)
#Vp_l = [6.,7.2]
#Vs_l = [3.5,4.1]
#rho_l = [2.72,3.03]
#h_l = [6.,24.]

#result = Ps_t_amp(Vp_l,Vs_l,h_l,rho_l,inci,0)

inci = np.deg2rad(15)
Vs_l = [3.7, 4.4, 4.1]
Vp_l = [6.5, 8, 7.2]
rho_l = [2.85, 3.4, 3.2]
h_l = [5., 65., 30.]

result = Ps_t_amp(Vp_l,Vs_l,h_l,rho_l,inci,0)