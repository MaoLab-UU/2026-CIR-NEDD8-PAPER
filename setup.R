library(dplyr)
library(Seurat)
library(readxl)
library(readr)
library(stringr)
library(scales)

lapply(list.files("C:/github 2026.2.9 submission/formal/R tools",full.names = T,pattern = '.R'),source) %>% invisible()


options(future.globals.maxSize= 32000*1024^2)

options(timeout = 6000)

setwd('R tools')

easyBuild()
