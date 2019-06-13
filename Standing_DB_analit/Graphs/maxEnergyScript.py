# -*- coding: utf-8 -*-
"""
Created on Tue May 21 19:31:00 2019

@author: mohit
"""

import os
#import numpy as np
import pandas as pd

n = 3
dx = 0.25
i = 0.5
path = 'C:\\Users\\mohit\\Desktop\\Mohit\\DB\\Moving_DB_analit\\Moving_DB_analit\\Hard'
os.chdir(path)
d = {}

while i<=n:
    if(i - int(i) != 0):
        os.chdir('.\\' + str(i)[0] + '_' + str(i)[2:])
    else:
        os.chdir('.\\' + str(i)[0])
    
    dataset = pd.read_csv('newResults_l.txt',sep=' ',header=None)
    d[dataset.iloc[0,1]] = dataset.iloc[0,4]
    
    os.chdir('..')
    i = i + dx
    
d = pd.DataFrame(list(d.items()) ,columns=['ADB', 'MaxEnergy'])