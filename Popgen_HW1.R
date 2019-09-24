#Popgen HW 1


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
#ge is all the expected genotypes
ge <- c(123.5,29.88,2.69,121.51,36.45,17.93)
ge
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
#The chi squared value is 29.155




