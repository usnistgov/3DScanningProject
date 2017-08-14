close all; clear all;

%Read the data file 
file1 = fullfile(pwd(),'testSphereData.xyz')
data1 = dlmread(file1);

knownRadius = 0.05; %50 mm
tic;
resultsCon1 = sphereFitLSQ1_conR(data1,knownRadius)
toc
tic
resultsUnc1 = sphereFitLSQ1_uncR(data1)
toc;

tic;
resultsCon2 = sphereFitLSQ1_constrained(data1,knownRadius)
toc
tic
resultsUnc2 = sphereFitLSQ1(data1)
toc;