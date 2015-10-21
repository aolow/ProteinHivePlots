# Kinase-Substrate Connection Hive Plots

This is a function with small example dataset for visualizing kinase-peptide-substrate connections in the human kinome.

This script relies on HiveR package:
  Bryan A. Hanson (2015) HiveR: 2D and 3D Hive Plots for R. R package version 0.2.44,
  academic.depauw.edu/~hanson/HiveR/HiveR.htm
  
  
# Example AKT1 hive plot

AKT1 can act as a kinase or as a substrate. Edges connect from Kinase axis to distinct peptides that the protein phosphorylates; which belong to substrate proteins (edge from Peptide to Substrate). In the absence of known target peptides, edges connect directly from kinase to substrate - yellow edges

![alt tag](https://github.com/aolow/RProteinHivePlots/blob/master/AKT1_hive.jpg)
