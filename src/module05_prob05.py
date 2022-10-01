# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 16:30:23 2022

@author: angel
"""

import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

files= glob.glob("../data/*EL.csv")
files_old = glob.glob("../data/*XV_old.csv")

for i in range(0, len(files)):
    df = pd.read_csv(files[i])
    r = df['a']*(1.0-df['e']**2)/(1.0+df['e']*np.cos(df['true_anom']))
    h = (df['a']*(1.0-df['e']**2))**1/2
    
    x = r*(np.cos(df['lon_asc_node'])*np.cos(df['arg_peri']+df['true_anom'])) - np.sin(df['lon_asc_node']*np.sin(df['arg_peri']+df['true_anom'])*np.cos(df['I']))
    y = r*(np.sin(df['lon_asc_node'])*np.cos(df['arg_peri']+df['true_anom'])) + np.cos(df['lon_asc_node']*np.sin(df['arg_peri']+df['true_anom'])*np.cos(df['I']))
    z = r*(np.sin(df['I'])*np.sin(df['arg_peri']+df['true_anom']))
    
    vx = x*h*df['e']/r*(df['a']*(1.0-df['e']**2))*np.sin(df['true_anom']) - (h/r)*(np.cos(df['lon_asc_node'])*np.sin(df['arg_peri']+df['true_anom'])+np.sin(df['lon_asc_node'])*np.cos(df['arg_peri']+df['true_anom'])*np.cos(df['I']))
    vy = y*h*df['e']/r*(df['a']*(1.0-df['e']**2))*np.sin(df['true_anom']) - (h/r)*(np.sin(df['lon_asc_node'])*np.sin(df['arg_peri']+df['true_anom'])-np.cos(df['lon_asc_node'])*np.cos(df['arg_peri']+df['true_anom'])*np.cos(df['I']))
    vz = z*h*df['e']/r*(df['a']*(1.0-df['e']**2))*np.sin(df['true_anom']) + (h/r)*(np.sin(df['I']*np.cos(df['arg_peri']+df['true_anom'])))
    
    
    cart= {'t':df['t'],'xh':x, 'yh':y, 'zh':z,'vxh':vx, 'vyh':vy, 'vz':vz}
    
    df_cart = pd.DataFrame(data=cart)
    
    df_cart.to_csv(files[i].replace("EL", "XV"))
    
    mag_pos = np.sqrt(x**2+y**2+z**2)
    mag_vel = np.sqrt(vx**2+vy**2+vz**2)
    
    df_old = pd.read_csv(files_old[i])
    mag_pos_old = np.sqrt(df_old['xh']**2 + df_old['yh']**2 + df_old['zh']**2 )
    mag_vel_old = np.sqrt(df_old['vxh']**2 + df_old['vyh']**2 + df_old['vz']**2)
    
    pos_diff = mag_pos - mag_pos_old
    vel_diff = mag_vel - mag_vel_old
    
    plt.plot(df['t'],pos_diff)
    plt.title("Difference in Magnitude of Position Vector vs t")
    plt.savefig('../plots/pos_vs_t_'+files[i].split('\\')[-1][:-4]+'.png',bbox='tight',dpi=300)
    plt.clf()
    plt.plot(df['t'],vel_diff)
    plt.title("Difference in Magnitude of Velocity Vector vs t")
    plt.savefig('../plots/vel_vs_t_'+files[i].split('\\')[-1][:-4]+'.png',bbox='tight',dpi=300)
    plt.clf()