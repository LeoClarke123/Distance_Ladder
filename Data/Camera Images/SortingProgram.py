# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 17:01:52 2022

@author: ryanw

not really sure what this will do yet :)

polar is 0 at 'due north' [(0,0) in Up image].
    Positive polar direction is down 
equatorial [equat] is 0 at the 'front' [(0,0) in the front image, or (0, 90) in polar coords]. 
    Positive equatorial direction is clockwise when viewed from polar north looking down. 
    That is, rotating right is increasing equatorial angle. 
"""

from numpy import *
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))  #finds the path of this program to use later
totalfuzzy = open(dir_path+"/"+"total fuzzy data.txt", "w")         #opens and overwrites existing data. Creates this file if not there
totalpoints = open(dir_path+"/"+"total point-like data.txt", "w")       #opens and overwrites existing data. Creates this file if not there

totalfuzzy.write("Name X Y BlueFlux GreenFlux RedFlux Size RadialVelocity \n")      #
totalpoints.write("Name X Y BlueFlux GreenFlux RedFlux Parallax RadialVelocity \n")

for i in ["Back", "Down", "Front", "Left", "Right", "Up"]:      #each of the six image directions
    for j in ["A", "B", "C", "D", "E", "F"]:                    #each of the image column subdivisions
        for k in ["01", "02", "03", "04", "05", "06"]:          #each of the image row subdivisions
            folder = dir_path + "/" + i + "/" + j + k           #sorts through the data to find data for each subdivision
            fuzzy = open(folder+"/fuzzy.txt", 'r')
            points = open(folder+"/points.txt", 'r')
            with fuzzy as f:
                fuzzycontents = f.readlines()[1:]                    #reads the file and puts all content (except line 1) into the list 'contents'
                for row in fuzzycontents:
                    [name, xpos, ypos, bluef, greenf, redf, size, veloc] = row.split()      #splits the row into strings of each data type instead of one long string
                    xpos, ypos, bluef, greenf, redf, size, veloc = float(xpos), float(ypos), float(bluef), float(greenf), float(redf), float(size), float(veloc)
                    #below code converts all data coordinates into a unified, spherical coordinate system
                    if i == "Up":
                        polar = abs(ypos)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Down":
                        polar = 180 - abs(ypos)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Front":
                        polar = abs(ypos - 90)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Back":
                        polar = abs(ypos - 90)
                        equat = 180 + xpos
                    elif i == "Right":
                        polar = abs(ypos - 90)
                        equat = 90 + xpos
                    else:
                        polar = abs(ypos - 90)
                        equat = 270 + xpos
                    #below code arranges data into a list string, and then adds it to the end of the document. 
                    addendum = str([name, round(equat, 3), round(polar, 3), bluef, greenf, redf, size, veloc]) + "\n"
                    totalfuzzy.write(addendum)
            with points as p:       #this section is identical to the fuzzy section, but replaces size variable with parallax angle variable
                pointcontents = p.readlines()[1:]                    #reads the file and puts all content (except line 1) into the list 'contents'
                for row in pointcontents:
                    [name, xpos, ypos, bluef, greenf, redf, parallax, veloc] = row.split()      #splits the row into strings of each data type instead of one long string
                    xpos, ypos, bluef, greenf, redf, parallax, veloc = float(xpos), float(ypos), float(bluef), float(greenf), float(redf), float(parallax), float(veloc)
                    if i == "Up":
                        polar = abs(ypos)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Down":
                        polar = 180 - abs(ypos)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Front":
                        polar = abs(ypos - 90)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Back":
                        polar = abs(ypos - 90)
                        equat = 180 + xpos
                    elif i == "Right":
                        polar = abs(ypos - 90)
                        equat = 90 + xpos
                    else:
                        polar = abs(ypos - 90)
                        equat = 270 + xpos
                    addendum = str([name, round(equat, 3), round(polar, 3), bluef, greenf, redf, parallax, veloc]) + "\n"
                    totalpoints.write(addendum)