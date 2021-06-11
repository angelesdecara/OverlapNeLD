#OverlapNeLD

Here we provide the codes to replicate the results of Garcia-Cortes et al. (2021). Our methods are developed in Fortran90, and thus can be compiled with gfortran or ifort, depending on the OS. We provide binaries for linux which end in .ix or in Mac (Power PC) which end in .mac. Should you encounter difficulties when running them, please let us know.

##Estimating Ne from r2 and life history data

Firstly, we need our data in a format suitable for [Ne Estimator](http://www.molecularfisherieslaboratory.com.au/neestimator-software/).
In the case of the sparrows, (available from https://datadryad.org/stash/dataset/doi:10 the data were in VCF format.
We used [PGDSpider] (http://www.cmpg.unibe.ch/software/PGDSpider/) (Lishcher and Excoffier 2012), and converted the VCF into GENO format.

For this dataset, no recombination map was available, thus we assumed all markers were independent. Otherwise, markers at similar recombination distances would have to be analysed independently. [Ne Estimator](http://www.molecularfisherieslaboratory.com.au/neestimator-software/) provides the observed and expected r2, where the latter is r2 due to finite sample size.

Using N0NeG.f90 (in folder r2cNeOverlap), we can obtain the estimated num- ber of newborns and the effective population size. This program can be compiled with gfortran, but we provide two binaries, one for Linux (N0NeG.ix) and one for Mac OSX (N0NeG.mac not intel). It is run on the command line

```bash
./N0NeG.ix life history.txt
```

where the file life history.txt includes the following:
1. Name of output file
2. 1 for detailed output, 0 otherwise
3. number of parental cohorts
4. survival for each cohort
5. contribution of each cohort but last, which is calculated as $1-\sum_{1}{n-1}p_i$
6. name of file with recombination distancec and the observed minus expected $r^2$ (two examples are provided, c_r2.txt and c0.5_r2.txt, the latter obtained with the sparrow data).

##Standard error via jackknife

Using the vcf file, we separated each sample in an indexed bcf file using

```
mkdir singlesample
bcftools +split allpug.recode.vcf.gz -Ob -o singlesample/
cd singlesample
for i in `(ls *bcf)`; do bcftools index $i; done
```

Given that there were 102 pugetensis sparrows, we created 102 “synthetic” populations using the following bash script:

```bash
#!/bin/bash
bcfs=(‘ls ../singlesample/*bcf‘) for p in {0..101};do
    newpop=()
    i=0
    for i in {0..101};do
if [ $i != $p ]; then newpop+=(${bcfs[$i]} )
fi done
    echo ${newpop[@]}
bcftools merge ${newpop[@]} > pop${p}.vcf done
```
We converted them in a format suitable for [Ne estimator](http://www.molecularfisherieslaboratory.com.au/neestimator-software/) using one of the scripts provided with [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/).

```bash
#!/bin/bash
for i in {0..101};do
    ~/Downloads/jdk1.7.0_80/bin/java -jar ~/Downloads/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputformat VCF -outputformat GENEPOP -inputfile pop${i}.vcf -outputfile pop${i}.geno -spid ../../ciNuttalli/template_VCF_GENEPOP.spid
done
```

Note the path will probably vary, as well as your java version. Finally, we run in this directory Ne estimator (linux version) for each population in command line mode. We created one input file for population 0 which we modified accordingly
```bash
#!/bin/bash
./Ne2-1L c:inputNeEstpop0.txt 
for i in {1..101};do
    j=$((i-1))
    sed -i "s/pop${j}/pop${i}/g" inputNeEstpop0.txt
    ./Ne2-1L c:inputNeEstpop0.txt
done
```

We selected the 2 lines where observed and expected $r^2$ are using
```
grep -A 2 Over outputsfromNeEstimator >> r2neNeEst.txt
```

Finally, our code CalcNfromR2sparrow.f90 provided standard errors from this jackknife for our estimates of $N_e$ and for those obtained from Ne estimator for three MAF. We provide the required files for this calculation for pugetensis sparrows in folder Sparrow/ciPugetensis.

##Multinomial simulations, independent replicates of one pair of loci
The code is provided in folder InferNfromL and is called InferNfromL.f90. Two binaries are provided, one for Mac OSX (Power PC), (infer.mac) and another for linux (infer.ix).

##Multinomial simulations, independent replicates with genomic correlations
The codes for these simulations are in folder ManyLociChrom. Program CoeffVarFast.f90 provides results for 1000 replicates which are then analysed with program sd.f90, which provides the mean and the variance between replicates. Binaries for mac and linux are provided (.mac and .ix, respectively).

##Plots

Most of the plots in the manuscript have been produced with [gnuplot](http://www.gnuplot.info/). Codes and binaries to produce fig. 3 from the manuscript are in folder triangle. First, run Losdos_triangulo.f90, then haceplot.f90, and lastly the script plot will provide triangle.png, as follows in Linux

```bash
./losdos.ix
./haceplot.ix
gnuplot plot
```
