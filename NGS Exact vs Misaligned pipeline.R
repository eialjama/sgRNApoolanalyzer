library("plyr")
library("dplyr")
library("readxl")
library("stringr")
library("writexl")
library("Rsamtools")

setwd("C:/Users/matthew.woods/Desktop/Ion Torrent Read Counter and Aligner/Ion Read Counter and Aligner")


source('Ion Torrent Read Counter.R') #initial filtering to determine any matching sequences to reference file
source("Density plot Final.R") #Density and distribution plot generation
source("NGS Exact vs Misaligned CrisprID Counter.R") #Takes count based off gene
source("NGS Exact vs Misaligned Sequence Counter.R") #Takes count based off sequence 
source("NGS CRISPR and Sequence Exact and Misaligned comparison.R") #takes threshold of mutation rate to perfect alignment rate and compile results from two previous scripts onto one spreadsheet

