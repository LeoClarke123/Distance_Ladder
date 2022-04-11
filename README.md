# PHYS3080 Distance Ladder Project
This GitHub repo serves as the home for all of the data and subsequent analysis for the mock-universe supplied as part of the PHYS3080 Distance Ladder project (c. Sem 1 2022). Ryan White, Ciaran Komarakul-Greene and Leo Clarke contributed to this project primarily through Spyder and Jupyter python interpreters. 
# If you're a visitor...
- Take a look at the Cluster Analysis folder for graphs and programs related to nearby galaxies, and some for the galaxies in the distant universe. 
- Take a look at HR_Diagram to look at the benchmark HR diagram from nearby stars, or the vertically shifted HR diagrams for nearby galaxies.
- Head to Period_Lum to look at graphs relating the luminosity of variable stars to their periodicity. 
- Go to \Data for some graphs of stellar temperatures / galaxy velocities vs distance 
- Go to \Data\Camera Images to see all-sky maps of stars and galaxies (redshifted or true colour). 
# Program Running Order
Since many of the programs rely on data previously generated from other programs, the .py and .ipynb files must be run in a certain order (if being run for the first time). It's my impression that this order should be:
1. SortingProgram.py  (auth. Ryan)
2. starmap.py (or starmap.ipynb)  (auth. Ryan)
3. Star Cluster Population Program.py   and   Fuzzy Cluster Identifier.py   (auth. Ryan)
4. clusterDistances.py  and   period_lum_analysis.ipynb   (auth. Ciaran and Leo resp.)
5. StarProperties.py  and   X-Ray Origin Finder.py  (auth. Ryan)
6. mainAnalysis.py  (auth. Ryan)
This order is intended as a guideline mostly and follows the rough order in which the programs were written. 
