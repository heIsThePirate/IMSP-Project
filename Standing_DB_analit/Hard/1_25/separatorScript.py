# -*- coding: utf-8 -*-
"""
Created on Tue May 14 16:57:19 2019

@author: mohit
"""

import re
l = ['Energies_l.txt', 'AvgKE_l.txt']

for filename in l:
    filename_2 = 'new' + filename
    file = open(filename, 'r')
    newFile = open(filename_2, 'w')

    for line in file:
        newLine = re.sub(' +', ' ', line)
        newFile.write(newLine)
    file.close()
    newFile.close()