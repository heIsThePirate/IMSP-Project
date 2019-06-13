# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:36:20 2019

@author: mohit
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter

path = 'C:\\Users\\mohit\\Desktop\\Mohit\\DB\\Moving_DB_analit\\Soft'
os.chdir(path)

i = 0.3
dx = 0.05
n = 0.6
l = ['newEnergies.csv']
d = {}

while i<=n:
    os.chdir('.\\' + str(i)[0] + '_' + str(i)[2:])
    for filename in l:
        dataset = pd.read_csv(filename)
        X = dataset.iloc[:, 1].values
        y = dataset.iloc[:, 6].values
        yhat = savgol_filter(y, window_length = 201, polyorder = 3)
        yhat = savgol_filter(yhat, window_length = 51, polyorder = 3)
        plt.plot(X, savgol_filter(y, window_length = 201, polyorder = 3))
        #plt.plot(X, yhat)
        #m = max(y)
        #d[i] = m
    os.chdir('..')
    i = i + dx
    i = round(i, 2)
