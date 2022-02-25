# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 17:01:52 2022

@author: ryanw

not really sure what this will do yet :)
"""

from numpy import *
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))

for i in ["Back", "Down", "Front", "Left", "Right", "Up"]:
    for j in ["A", "B", "C", "D", "E", "F"]:
        for k in ["01", "02", "03", "04", "05", "06"]:
            folder = dir_path + "/" + i + "/" + j + k
            fuzzy = open(folder+"/fuzzy.txt", 'r')
            with fuzzy as f:
                contents = f.readlines()
            