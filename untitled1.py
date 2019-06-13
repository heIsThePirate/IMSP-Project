# -*- coding: utf-8 -*-
"""
Created on Wed May 29 17:04:04 2019

@author: mohit
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 14 16:57:19 2019

@author: mohit
"""

import os
import re
import pandas as pd
l = ['DB_list.txt']
path = 'C:\\Users\\mohit\\Desktop\\Mohit\\DB\\Moving_DB_analit\\1_2'
os.chdir(path)
for filename in l:
    filename_2 = 'new' + filename
    file = open(filename, 'r')
    newFile = open(filename_2, 'w')
    for line in file:
        newLine = re.sub(' +', ' ', line)
        newFile.write(newLine)
    file.close()
    newFile.close()
dataset = pd.read_csv(filename_2, sep = " ", header = None)
dataset = dataset.iloc[:, 1:]
dataset.to_csv(filename_2 + '.csv')