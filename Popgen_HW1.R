#Popgen HW 1


#1B 
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




