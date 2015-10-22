#!/bin/bash 


#Generate the correlation matrix for all SNPs
plink --bfile example/dat/example --r --matrix --out example/result/cor

#Calculate the loading of SNPs. The parameters are 1)the number of SNPs, 2) number of principal components and 3) input directory
nSNP=`wc -l ./example/dat/example.bim | awk '{print $1}'`
./getPCA $nSNP 100 ./example/result/
	
#Project the samples to the eigen vectors
plink --bfile example/dat/example --recodeA --out example/result/geno
nSample=`wc -l ./example/dat/example.fam  | awk '{print $1}'`
./getMult $nSNP 100 $nSample ./example/result/


