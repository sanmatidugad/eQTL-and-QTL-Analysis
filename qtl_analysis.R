rm(list = ls())    # clear all variables in ym environment. 
#install.packages("qtl")
library(qtl)
setwd("/home/data/4.QTL_Analysis/")
#listeria = read.cross("csv", ".", file = "listeria.csv",genotypes = c('A', 'H', 'B', 'C'), na.strings = "-")
data("listeria")

## knowing your data
  #jittermap(listeria)
  summary(listeria)
  names(listeria$pheno)   # what phenotypes I have
  names(listeria$geno)    # what chromosomes I have
  names(listeria)
  
  nphe(listeria)    # how many phenotypes
  nchr(listeria)    # how many chromosomes
  
  table(listeria$pheno$sex)
  #table(listeria$pheno$pgm)
  table(listeria$pheno$T264)

plot(listeria, pheno.col = c(1))
plot.map(listeria)
#plot(listeria, pheno.col = c(1, 2))

## Explanation of plots
## genetic map shows markers, Each tick is a marker. For example chromosome 1 has 13 markers (13 columns) 
## missing genotypes shows a lot of missing data on each marker, which is okay. 
## the histogram shows frequency of values for column T264. There are 35 values which are 264 exactly. 

# To save graphs
jpeg("genetic_map.jpeg")
plot.map(listeria)
graphics.off()

#install.packages("vioplot")
library(vioplot)
vioplot((listeria$pheno$T264), xlab = "sex", ylab = "T264", names = "female")
lines(c(0.75, 1.25), c(100,100), col = "blue", lwd = "5")
lines(c(0.75, 1.25), c(264,264), col = "salmon", lwd = "5")
# Violin plot is just a histogram but turned 90 degrees. 

estimated_map = est.map(listeria)
plot.map(listeria, estimated_map, main = "Comparison of Genetic Maps") ## It shows what the actual map is compared with the estimated map. 
## The above plot shows statistical similarity. Showing % error calculation. Based on Linear Regression
## Its kind of a diagnostic tool

## Pairwise recombination fractions and LOD scores. 
listeria = est.rf(listeria)
plotRF(listeria)
## The X-axis is number of markers v/s number of markers.
## The bottom trianlgle shows lod score. If the bottom triangle has hotspots there are errors in data.


## IN TERMS OF DIAGNOSTICS - CALCULATE LOG SCORES. -- Read more about lod scores. 
listeria = calc.errorlod(listeria, error.prob = 0.01) ## look at genotypes and ask if are there any ones that look like they might have statistical error of probabiltiy of 0.01 error 
top.errorlod(listeria, cutoff = 4)

#plotGeno(listeria, chr = 1)
plotGeno(listeria, chr = 3, include.xo = FALSE, ind = 1:30)

  ## Marker Genotypes
  geno.image(listeria, chr=1)
  geno.image(subset(listeria, ind = 1:100), chr = 1)    
  # vertically check what is the probability of having a geneotype near that loci.
  # Each color is a genotype. (We have A,B,C,H -- four colors.)
#  plot(listeria$pheno$T264)    # scatter plot

## Main Scan - Looking at analysis of each chromosomes. -- Genome Scan
listeria = clean(listeria)    # remove the statistics done on listeria. Just to clean any intermediate calculations done

## Calculating genotype probability
listeria = calc.genoprob(listeria, step = 2.0, off.end = 0.0, error.prob = 1.0e-4,
                         map.function = "haldane", stepwidth = "fixed")    # step = 2.0 Centimorgans (cM) are a unit of genetic distance used to measure the relative positions of genetic loci along a chromosome
## build a simulated genotype, to compare data -- not necessary to do this step.
listeria = sim.geno(listeria, n.draws = 32, step = 2, off.end = 0.0, error.prob = 1.0e-04,
                    map.function = "haldane", stepwidth = "fixed")
## run your first main scan
listeria.scan1 = scanone(listeria, pheno.col = 1, model = "normal", method = "em")
listeria.scan1.perm = scanone(listeria, pheno.col = 1, model = "normal", method = "em", n.perm = 1000)   # basically generate a map that we can look at

threshold1 = summary(listeria.scan1.perm, alpha = c(0.63, 0.20, 0.10))    # set up some threshold to compare data, 37%, 80% and 90%
threshold1

options(repr.plot.width = 16)
plot(listeria.scan1, main = "Mainscan plot of T264") 
abline(h = threshold1[1], lty = 2, lwd = 1, col = "blue")
abline(h = threshold1[2], lty = "dotted", lwd = 1, col = "red")    # number or name anything works
abline(h = threshold1[3], lty = "dotted", lwd = 1, col = "green")
## So what does the above plot say -- We have a pretty big QTL on chromosome 5, chromosome 13 and chromosome 15, that is causing effects on T264 phenotype

## Get a summary of all this data
summary(listeria.scan1, perms = listeria.scan1.perm, lodcolumn = 1, alpha = 0.68)
summary(listeria.scan1, perms = listeria.scan1.perm, lodcolumn = 1, alpha = 0.20)
summary(listeria.scan1, perms = listeria.scan1.perm, lodcolumn = 1, alpha = 0.10)
# Chromosome 1,5,6,12,13,15 have some information, And the high lod scores are present on c5.loc28 , D13M147 markers

## PLOTTING QTL Effects
# Plot of trait by three possible genotypes on the chromosome 5 locus. It appers to be a dominant locus)
plotPXG(listeria, marker = "D5M357", pheno.col = log2(listeria$pheno$T264))

marker.name1 = find.marker(listeria, chr = 5, pos = 28.0)
  marker.name1
marker.name2 = find.marker(listeria, chr = 13, pos = 26.2)
  marker.name2
marker.name3 = find.marker(listeria, chr = 6, pos = 59.4)
  marker.name1

effectplot(listeria, pheno.col = 1, mname1 = marker.name1)
effectplot(listeria, pheno.col = 1, mname1 = marker.name2)
effectplot(listeria, pheno.col = 1, mname1 = marker.name3)

## Confidence Interval
CIchr5 = bayesint(listeria.scan1, chr = 5, prob = 0.95)
CIchr5

plot(listeria.scan1, chr = 5, lodcolumn = 1, main = "Confidence Interval for Chromosome 5")
lines(x=CIchr5[c(1,3),2], y = c(0,0), type = "l", col = "green", lwd=4)
## The above plot hence shows that something is happening in Chr 5 from regions 18 to 36 that is causing T264 to increase.

## Histogram of 
hist(as.vector(listeria.scan1.perm), breaks = seq(0,30, by=0.5),
     freq = FALSE, ylim = c(0,0.5), xlab = "LOD Score", main = "Histogram of permutation max LOD score vs Chi-square distribution")
lines(seq(0,30,by=0.1), dchisq(seq(0,30,by=0.1)/log(10)*2, df = 2))
##distribution of maximum lod score from the shuffled dataset, because we have a plain chisq, everything beyond 5 LOD score we consider unusual. 

#Another QTL Analysis
## Adjusting for Covariates.
## We can adjust for sex analysis if there are differences by sex in the data.
# We can then perform the permutations stratified by sex.

## Strengths and Limitations for QTL Mapping
#Strengths
#1. Identifies genetic regions linked to traits.
#2. Helps understand complex traits influenced by multiple genes and the environment.
#3. Allows cross-species comparisons for evolutionary insights.
#4. Useful in breeding programs for crop and livestock improvement.

#Limitations
#1. Limited to identifying genetic regions, not understanding gene interactions.
#2. Requires large sample sizes to avoid false positives/negatives.
#3. Findings may not generalize universally across populations or species.
#4. Time-consuming, expensive, and resource-intensive.







