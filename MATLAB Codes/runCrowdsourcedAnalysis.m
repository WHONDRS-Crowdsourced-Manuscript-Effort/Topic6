function runCrowdsourcedAnalysis

clear all; close all; clc;

filename = 'Processed_S19S_Sediments_Water_2-2_newcode.csv';
tbl = readtable(filename);
dataDescrp = filename;
tblOut = runLambda_v0p31_crowdsourced(tbl,7,'y',dataDescrp);
assignin("base","tblOut",tblOut)

end