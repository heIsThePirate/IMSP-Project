# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 13:33:03 2019

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

while i<=n:
    os.chdir('.\\' + str(i)[0] + '_' + str(i)[2:])
    for filename in l:
        dataset = pd.read_csv(filename)
        X = dataset.iloc[:, 1].values
        y = dataset.iloc[:, 4].values
        yhat = savgol_filter(y, window_length = 201, polyorder = 3)
        yhat = savgol_filter(y, window_length = 51, polyorder = 3)
        plt.plot(X, yhat)
        dataset.iloc[:, 4] = yhat
        dataset.to_csv(filename[: -4] + '4.csv')
    os.chdir('..')
    i = i + dx
    i = round(i, 2)
