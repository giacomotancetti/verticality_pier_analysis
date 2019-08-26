#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 24 21:24:14 2019

@author: giacomo tancetti
"""
import os
import csv
import numpy as np
import math
from scipy.optimize import leastsq

path='./coordinates'
files=os.listdir(path)

def f_min(X,p):
    plane_xyz = p[0:3]
    distance = (plane_xyz*X.T).sum(axis=1) + p[3]
    return distance / np.linalg.norm(plane_xyz)

def residuals(params, signal, X):
    return f_min(X, params)

incl_angle={}

for filename in files:   
    pier_name=filename[0:-4]
    type_surf=filename[5:9]  # check if pier face is longitudinal or transverse
    
    coord_str=[]
    with open(path+'/'+filename,'r') as csvfile:
        csv_cont=csv.reader(csvfile, delimiter=';')
        next(csv_cont, None)
        for row in csv_cont:
            coord_str.append(row)
      
    x=[]
    y=[]
    z=[]
    for row in coord_str:
        x.append(float(row[0]))
        y.append(float(row[1]))
        z.append(float(row[2]))

    XYZ=np.array([x,y,z])   
    # Inital guess of the plane
    # calculate coordinates of mean point of the plane
    xp=XYZ[0].mean()
    yp=XYZ[1].mean()
    zp=XYZ[2].mean()
    P=[xp,yp,zp]

    # calcolo parametri direttori del piano
    alpha_deg=9.7407           # angolo fra sezione trasversale e asse E (+ antiorario)
    alpha_v_long=alpha_deg+270      # angolo fra versore piano e asse E (+ antiorario)
    alpha_v_tras=alpha_v_long+90
    alpha_v_long_rad=math.radians(alpha_v_long)
    alpha_v_tras_rad=math.radians(alpha_v_tras)

    if type_surf == 'long':
        v=(math.cos(alpha_v_long_rad),math.sin(alpha_v_long_rad),0)
        d= -(v[0]*P[0]+v[1]*P[1]+v[2]*P[2])
    elif type_surf == 'tras':
        v=(math.cos(alpha_v_tras_rad),math.sin(alpha_v_tras_rad),0)
        d= -(v[0]*P[0]+v[1]*P[1]+v[2]*P[2])
    
    p0 = [v[0], v[1], v[2], d]

    sol = leastsq(residuals, p0, args=(None, XYZ))[0]
    incl=math.degrees(math.acos(sol[2]))

    print("Solution: ", sol)
    print("Old Error: ", (f_min(XYZ, p0)**2).sum())
    print("New Error: ", (f_min(XYZ, sol)**2).sum())
    
    incl_angle[pier_name]=[-(90-incl),(f_min(XYZ, sol)**2).sum()]

