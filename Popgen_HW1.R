#Popgen HW 1
#1B
#The genotype frequences can be represented as by the following, where the allele frequencies of S, I and G are represented by fs, fi, and fg, repectively. 
#GSS = fs^2
#GFF =ff^2
#GII = fi^2
#GSF = 2fs*ff
#GSI = 2fs*fi
#GFI = 2ff*fi
#1B 
#tot represents the total population number
tot <- 141 + 111 + 15 + 28 + 32 + 5
tot
#The population total is 332 individuals 
#allelefreq is a function to calculate the allele frequencies where x represents homozygotes where y and z represent the two different heterozygotes.
allelefreq <- function(x,y,z) {(x + (y + z)/2)/tot}
#fs represents the allele frequency of S
fs <- allelefreq(141,111,15)
fs
#fs is 0.6144578
#ff is the allele frequency for the F allele 
ff <- allelefreq(28,32,111)
ff
#ff is 0.2996988
#fi is the allele frequency of I
fi <- allelefreq(5,32,15)
fi
#fi is 0.08584337
#fs, fi, and ff need to add to 1
fs + fi + ff
#and they do :) 
#And now to calculate genotype frequencies. 
GSS = fs^2
GSS
#GSS is 0.3775584
GFF =ff^2
GFF
#GFF is 0.08981937
GII = fi^2
GII
#GII is 0.007369085
GSF = 2*fs*ff
GSF
#GSF is 0.3683045
GSI = 2*fs*fi
GSI
#GSI is 0.1054943
GFI = 2*ff*fi
GFI
#GFI is 0.05145431
#These genotype frequences should add to 1
GSS + GFF + GII + GFI +GSF +GSI
#AND THEY DO!!
#ge is the expected genotype totals of all individuals in the population
ge <- tot*c(GSS,GFF,GII,GSF,GSI,GFI)
ge
#egenotypenames is an object I'm going to use just to symbolically represent each expected genotype.
egenotypenames <- c("SS","FF","II","SF","SI","FI")
egenotypenames
egnames <- as.list(egenotypenames)
#Now I am going to represent the total amount of individuals per genotype with the following array: expectedgenotypes
expectedgenotypes <- data.frame(ge, row.names = egnames)
expectedgenotypes
#  Expected Genotypes
#SS 125.349398
#FF  29.820030
#II   2.446536
#SF 122.277108
#SI  35.024096
#FI  17.082831
#go is the observed genotypes
go <- c(141,28,5,111,15,32)
go
#dif is the difference between and expected genotypes
dif <- go-ge
dif
#dif 2 is the square values of the differences between observed and expected genotypes
dif2 <- dif^2
#prechi is the division for genotypes of (observed - expected) by the corresponding expected values
prechi <- dif2/ge
prechi
#chi2 is the sum of all the values in prechi (aka the chi squared value)
chi2 <- sum(prechi)
chi2
#The chi squared value is 30.24456
#There are 2 non-independent variables and 6 possible genotypes. So k, the degrees of freedome is:
k <- 6-1-2
k
#There are 3 degrees of freedom. 
#This chi squared value at a k of 3 statistically deviates from Hardy Weinberg. 




#Question 2 
#homo represents the frequency of homozygotes (based on the Hardy-Weinberg model) and hete represents the heterozygote frequency
homo<- function(p,q){(p^2)+(q^2)}
hete <- function(p,q){2*p*q}
homop0 <- homo(0,1)
homop0
homop0.25 <- homo(0.25,0.75)
homop0.25
homop0.5 <- homo(0.5,0.5)
homop0.5
homop0.75 <- homo(0.75,0.25)
homop0.75
homop1 <- homo(1,0)
homozygotes <-c(homop0,homop0.25,homop0.5,homop0.75,homop1)
homozygotes
heteq1 <- hete(0,1)
heteq0.75 <- hete(0.25,0.75)
heteq0.5 <- hete(0.5,0.5)
heteq0.25 <- hete(0.75,0.25)
heteq0 <- hete(1,0)
heterozygotes <- c(heteq0,heteq0.25,heteq0.5,heteq0.75,heteq1)
heterozygotes
pfreq <- list("0","0.25","0.5","0.75","1")
pfreq
genotype_frequencies <- data.frame(homozygotes, heterozygotes, row.names = pfreq)
genotype_frequencies
#p       homozygotes heterozygotes
#0          1.000         0.000
#0.25       0.625         0.375
#0.5        0.500         0.500
#0.75       0.625         0.375
#1          1.000         0.000
#Heterozygote frequencies are maximized when p (and also q) is 0.5. 

#Question 3 
#D = gOD-pOpD where D is the non-d allele. 
#gOD = 0.1 +(0.67)(0.4)
gOD <- 0.1 +0.67*0.4
gOD
#gOD = 0.368. O- individuals are homozygotes for O and the non-d allele (D). And gOD is 0.368 so:
Omindividuals <- gOD^2
Omindividuals
#There is a frequency of 0.135424 O- individuals in the population. 
#Question 4A
#total individuals will be represented with itot
itot <- 63+7+12+58
itot
#There are 140 individuals. gtot it the total number of gametes that went into making this population. 
gtot <- 2*itot
gtot
#There are 280 gametes. 
#mdmd = gmd^2, 
#mDmd = 2*(gmD)*(gmd),
#Mdmd = 2*(gMd)*(gmd), 
#MmDd (without recombination) = (1/2)*2*(MD)*(md)(1-r) + (1/2)*(Md)*(mD)*r

#gmd^2 = 58 , gmd = sqrt(58)
gmd <- sqrt(58/gtot)
gmd
#gmd = 7.615773
#2*(gMd)*(gmd) = 7
#7/(2*gmd) = gMd
gMd <- (7/gtot)/(2*gmd)
gMd
#gMd = 0.4595725
#12 = 2*(gmD)*(gmd), gmD = 12/(2*gmd)
gmD <- (12/gtot)/(2*gmd)
gmD
#gmD = 0.7878386 
gamfreqs <- c(gmd, gMd, gmD, "unknown")
gamnames <- list("gmd","gMd","gmD", "gMD")
gametes <- data.frame(gamfreqs, row.names = gamnames)
gametes
gMD <- 1 - (gmd+gMd+gmD)
gMD
#gMD = 0.4703234
#If gMD is 0.4703234, then MmDd expected (MmDdexp) should be 2*gMD*gmd. 
#MmDdobs is the observed frequency of MmDd 
MmDdobs <- 63/itot
MmDdobs
#MmDdobs is 0.45 
MmDdexp <- 2*gMD*gmd
MmDdexp
#The result is 0.4281161. 
#MDmdobs = gMD*gmd*(1-r) + gMd*gmD*r. I did some algebra ... 
r <- (MmDdobs - MmDdexp)/(gMd*gmD - MmDdexp)
r
# r is 0.05127152
#4B 
#Dt = Do((1-r)^t)
#Algebra --- t = ln(Dt/Do)/ln(1-r)
#t = ln(0.05/0.23)/ln(1-r)








