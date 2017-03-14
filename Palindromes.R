data <- read.table("hcmv.txt",header=TRUE)
#input library
library(ggplot2)
library(MASS)
## begin ##
index <- as.integer(data$location)
N <- 229354
n <- 296
gene <- seq(1, N)
###################    First Approach    ####################
# -------------------- random scatter ----------------------#
#############################################################
K<-40
tab1<- table(cut(data$location,breaks=seq(1,N,length.out=K+1),include.lowest=TRUE))
counts<- as.vector(tab1)
E_i<-n/K
chi_2<-sum((counts-E_i)^2/E_i)
p_value<-1-pchisq(chi_2,df=59)
p_value

K<-50
tab2<- table(cut(data$location,breaks=seq(1,N,length.out=K+1),include.lowest=TRUE))
counts<- as.vector(tab2)
E_i<-n/K
chi_2<-sum((counts-E_i)^2/E_i)
p_value<-1-pchisq(chi_2,df=59)
p_value

K<-60
tab3<- table(cut(data$location,breaks=seq(1,N,length.out=K+1),include.lowest=TRUE))
counts<- as.vector(tab3)
E_i<-n/K
chi_2<-sum((counts-E_i)^2/E_i)
p_value<-1-pchisq(chi_2,df=59)
p_value

K<-70
tab4<- table(cut(data$location,breaks=seq(1,N,length.out=K+1),include.lowest=TRUE))
counts<- as.vector(tab4)
E_i<-n/K
chi_2<-sum((counts-E_i)^2/E_i)
p_value<-1-pchisq(chi_2,df=59)
p_value

qplot(data$location,data=data,color = "blue",geom="histogram",
      xlab="Location on chain of DNA base pairs",
      ylab="Number of palindromes in the interval",binwidth=N/70,
      main="n.region = 70")

qplot(data$location,data=data,color = "red",geom="histogram",
      xlab="Location on chain of DNA base pairs",
      ylab="Number of palindromes in the interval",binwidth=N/40,
      main="n.region = 40")

qplot(data$location,data=data,color = "blue",geom="histogram",
      xlab="Location on chain of DNA base pairs",
      ylab="Number of palindromes in the interval",binwidth=N/50,
      main="n.region = 50")

qplot(data$location,data=data,color = "red",geom="histogram",
      xlab="Location on chain of DNA base pairs",
      ylab="Number of palindromes in the interval",binwidth=N/60,
      main="n.region = 60")

par(new=FALSE)
### results saved! 
## ------------------------END-------------------------------------

####################  Second Approach ##########################
## --------------------- Counts -----------------------------##
###############################################################

# First, we draw the histgram to compare the sample distribution and 
# the poisson distribution, we consider interval to be
# 4000 interval length, interval counts around 57
# 4587 interval length, interval counts around 50
# 5734 interval length, interval counts around 40
#---- Draw the histgram ---------------------------------------
#---- Let interval size be 4000, find out # of hits inside each interval
data.count <- table(cut(data$location, seq(min(data$location),max(data$location),
                                           4000),include.lowest = TRUE))
# Graphical comparison between sample and poisson distribution
hist(data.count, breaks = 15, col = rgb(1,0,0,0.5), probability = TRUE, 
     xlab = "number of points inside an interval", ylim = c(0,0.25))
lines(density(data.count, adjust = 2), col = rgb(1,0,0,0.5))
Pois <- rpois(1000, lambda = mean(data.count))
hist(Pois, breaks = 15, col = rgb(0,0,1,0.5), probability = TRUE, add = TRUE)
lines(density(Pois, adjust = 2), col = rgb(0,0,1,0.5))
legend(x = 10, y = 0.15, legend = c("Sample", "Poisson"), lty = c(1,1), 
       col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

#---- Let interval size be 4587, find out # of hits inside each interval
size <- 4587
data.count <- table(cut(data$location, seq(min(data$location),max(data$location), size),
                        include.lowest = TRUE))
# Graphical comparison between sample and poisson distribution
hist(data.count, breaks = 15, col = rgb(1,0,0,0.5), probability = TRUE, 
     xlab = "number of points inside an interval", ylim = c(0,0.25))
lines(density(data.count, adjust = 2), col = rgb(1,0,0,0.5))
Pois <- rpois(1000, lambda = mean(data.count))
hist(Pois, breaks = 15, col = rgb(0,0,1,0.5), probability = TRUE, add = TRUE)
lines(density(Pois, adjust = 2), col = rgb(0,0,1,0.5))
legend(x = 10, y = 0.15, legend = c("Sample", "Poisson"), 
       lty = c(1,1), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

#---- Let interval size be 5733, find out # of hits inside each interval
size <- 5733
data.count <- table(cut(data$location, seq(min(data$location),max(data$location), size),include.lowest = TRUE))
# Graphical comparison between sample and poisson distribution
hist(data.count, breaks = 15, col = rgb(1,0,0,0.5), probability = TRUE, xlab = "number of points inside an interval", ylim = c(0,0.25))
lines(density(data.count, adjust = 2), col = rgb(1,0,0,0.5))
Pois <- rpois(5733, lambda = mean(data.count))
hist(Pois, breaks = 15, col = rgb(0,0,1,0.5), probability = TRUE, add = TRUE)
lines(density(Pois, adjust = 2), col = rgb(0,0,1,0.5))
legend(x = 10, y = 0.15, legend = c("Sample", "Poisson"), lty = c(1,1), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

## function for count##
chisqtable <- function(n.region, site, N){
  n <- length(site)
  lambda.est <- n/n.region
  count.int <- table(cut(site, breaks = seq(1, length(gene), length.out=n.region+1), include.lowest=TRUE))
  count.vector <- as.vector(count.int)
  count.range <- max(count.vector) - min(count.vector) + 1
  table <- matrix(rep(NA, count.range*3), count.range, 3)
  for (i in 1:count.range){
    offset <- min(count.vector) - 1
    table[i, 1] <- i + offset
    table[i, 2] <- sum(count.vector == i + offset)
    if ((i + offset == min(count.vector)) && (min(count.vector) != 0))
      table[i, 3] <- ppois(i+offset, lambda.est)*n.region
    else if (i + offset == max(count.vector))
      table[i, 3] <- 1 - ppois(i + offset - 1, lambda.est)
    else
      table[i, 3] <- (ppois(i+offset, lambda.est) - ppois(i + offset - 1, lambda.est))*n.region
  }
  colnames(table) <- c("  Palindrome count","  Number of Observed","  Intervals expected")
  row
  table <- as.table(table)
  return (table)
}

#-------------- first try, intevel length = n.region = 50 ----------------------------
set.seed(2017)
n.region <- 50 # the interval number
site.random.tabtemp <- chisqtable(n.region, index, N)
site.random.tabtemp
lambda.est <- n/n.region
lambda.est
## read the table:
# lambda.set = 5.92
# [0,2] forming a group, [8,9] forming a group,[10+,...] forming a group
site.random.tab <- matrix(rep(NA, 8*2), 8, 2)
site.random.tab[1,] <- colSums(site.random.tabtemp[1:3, 2:3])
site.random.tab[2:6,] <- site.random.tabtemp[4:8, 2:3]
site.random.tab[7,] <- colSums(site.random.tabtemp[9:10, 2:3])
site.random.tab[8,] <- colSums(site.random.tabtemp[11:16, 2:3])
colnames(site.random.tab) <- c("  Number of Observed", "  Intervals expected")
rownames(site.random.tab) <- c("[0,2]","3","4","5","6","7","[8,9]","[10,15]")
site.random.tab <- as.table(site.random.tab)
site.random.tab
# Pearson's chi square-test (goodness of fit)
# Select a desired level of confidence. We will select alpha=0.05.
site.random.stats <- sum((site.random.tab[,2] - site.random.tab[,1])^2/site.random.tab[,2])
site.random.stats
# calculate the p-value
p_value <- 1 - pchisq(site.random.stats, 8 - 2, lower.tail=FALSE)
p_value
# p value = 0.1364952

Residuals <- (site.random.tab[,1] - site.random.tab[,2]) / sqrt(site.random.tab[,2])
plot(Residuals, type = 'h', ylab = "standardized residuals", xlab = "interval index",xaxt = "n",
     main = "n.region = 50, test based on poisson")
axis(1, at=c(1:8), labels = c("Counts [0,2]","Counts 3",
                              "Counts 4","Counts 5","Counts 6",
                              "Counts 7","Counts [8,9]",
                              "Counts [10,15]"))
#-------------------------- first try, END-----------------------------------------

#-------------- second try, intevel length = n.region = 57 ----------------------------
set.seed(2018)
n.region <- 57 # the interval number
site.random.tabtemp <- chisqtable(n.region, index, N)
site.random.tabtemp
lambda.est <- n/n.region
lambda.est
## read the table:
# lambda.set = 5.192982
# [9,13] forming a group
site.random.tab <- matrix(rep(NA, 9*2), 9, 2)
site.random.tab[1:8,] <- site.random.tabtemp[1:8, 2:3]
site.random.tab[9,] <- colSums(site.random.tabtemp[9:13, 2:3])
colnames(site.random.tab) <- c("  Number of Observed", "  Intervals expected")
rownames(site.random.tab) <- c("Counts 1","Counts 2",
                               "Counts 3","Counts 4",
                               "Counts 5","Counts 6",
                               "Counts 7","Counts 8",
                               "Interval Counts [9,13]")
site.random.tab <- as.table(site.random.tab)
site.random.tab
# Pearson's chi square-test (goodness of fit)
# Select a desired level of confidence. We will select alpha=0.05.
site.random.stats <- sum((site.random.tab[,2] - site.random.tab[,1])^2/site.random.tab[,2])
site.random.stats
# calculate the p-value
p_value <- 1 - pchisq(site.random.stats, 9 - 2, lower.tail=FALSE)
p_value
# p value = 0.9540093

Residuals <- (site.random.tab[,1] - site.random.tab[,2]) / sqrt(site.random.tab[,2])
plot(Residuals, type = 'h', ylab = "standardized residuals", xlab = "interval index",xaxt = "n",
     main = "n.region = 57, test based on poisson")
axis(1, at=c(1:9), labels = c("Counts 1","Counts 2",
                              "Counts 3","Counts 4",
                              "Counts 5","Counts 6",
                              "Counts 7","Counts 8",
                              "Counts [9,13]"))
#-------------------------- second try, END-----------------------------------------

#-------------- third try, intevel length = n.region = 40 ----------------------------
set.seed(2019)
n.region <- 40 # the interval number
site.random.tabtemp <- chisqtable(n.region, index, N)
site.random.tabtemp
lambda.est <- n/n.region
lambda.est
## read the table:
# lambda.set = 7.4
# [3,5] forming a group, [10,15] forming a group
site.random.tab <- matrix(rep(NA, 6*2), 6, 2)
site.random.tab[1,] <- colSums(site.random.tabtemp[1:3, 2:3])
site.random.tab[2:5,] <- site.random.tabtemp[4:7, 2:3]
site.random.tab[6,] <- colSums(site.random.tabtemp[8:13, 2:3])
colnames(site.random.tab) <- c("  Number of Observed", "  Intervals expected")
rownames(site.random.tab) <- c("[3,5]","6","7","8","9","[10,15]")
site.random.tab <- as.table(site.random.tab)
site.random.tab
# Pearson's chi square-test (goodness of fit)
# Select a desired level of confidence. We will select alpha=0.05.
site.random.stats <- sum((site.random.tab[,2] - site.random.tab[,1])^2/site.random.tab[,2])
site.random.stats
# calculate the p-value
p_value <- 1 - pchisq(site.random.stats, 6 - 2, lower.tail=FALSE)
p_value
# p value = 0.7969956

Residuals <- (site.random.tab[,1] - site.random.tab[,2]) / sqrt(site.random.tab[,2])
plot(Residuals, type = 'h', ylab = "standardized residuals", xlab = "interval index",xaxt = "n",
     main = "n.region = 40, test based on poisson")
axis(1, at=c(1:7), labels = c("Counts [0,3]","Counts 4",
                              "Counts 5","Counts 6",
                              "Counts 7","Counts [8,9]",
                              "Counts [10,15]"))
#-------------------------- third try, END-----------------------------------------

####################  Third Approach ##########################
## ------------------- Spacing Test--------------------------##
###############################################################
setwd("C:/Users/hao/Dropbox/courses/Math289C/homework 3")
getwd()
data <- read.table("hcmv.txt",header=TRUE)
spacing <- matrix(0, nrow = 295, ncol = 1)

#----------------------------spacing between consecutive palindromes----------------------------------------------------------------
for (i in 2:296)
{
  spacing[i-1,1] = data[i,1] - data[i-1,1]
}
lamda = 1/(mean(spacing))
hist(spacing, breaks = 30, col = 'lightgreen', xlim = c(0,6000), ylim = c(0, 100),
     main = "Distribution of spacings between pairs of palindromes")
#perform chi-square test
chisqtable <- function(spacing, b, rowname){
  
  count.int <- table(cut(spacing, breaks = b, include.lowest = TRUE))
  count.vector <- as.vector(count.int)
  count.range <- max(count.vector) - min(count.vector) + 1
  
  lambda.est <- 1/mean(spacing)
  
  table <- matrix(0, 12, 2)
  for (i in 1:12){
    table[i, 1] <- count.vector[i]
    table[i,2] <- 295*(pexp(b[i+1], lamda) - pexp(b[i], lamda) )
  }
  colnames(table) <- c("Observed counts","Expected counts")
  rownames(table) <- rowname
  table <- as.table(table)
  return (table)
}

#-----------first set of intervals------------------------------
b <- c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000)
rowname = c("0-500", "500-1000", "1000-1500", "1500-2000", "2000-2500", "2500-3000", "3000-3500", "3500-4000", "4000-4500", "4500-5000", "5000-5500", "5500-6000")
table1 = chisqtable(spacing, b, rowname)
table1

newTable1 <- matrix(rep(NA, 8*2), 8, 2)
newTable1[1:7,] <- table1[1:7, 1:2]
newTable1[8,] <- colSums(table1[8:12,1:2])
colnames(newTable1) <- c("Observed counts","Expected counts")
rowname = c("0-500", "500-1000", "1000-1500", "1500-2000", "2000-2500", "2500-3000", "3000-3500", "3500-6000")
rownames(newTable1) <- rowname
newTable1
stats = sum((newTable1[,1] - newTable1[,2])^2/newTable1[,2])
p_value1 = pchisq(stats, df = 6, lower.tail = FALSE)
p_value1
residuals <- (newTable1[,1] - newTable1[,2])/sqrt(newTable1[,2])
plot(residuals, type = 'h', ylab = "standardized residuals", xlab = "intervals",xaxt = "n", main = "Residual plot")
axis(1, at=c(1:8), labels = c("1st", "2ed", "3rd", "4th", "5th", "6th", "7th", "8th"))

##################################################################################
setwd("C:/Users/hao/Dropbox/courses/Math289C/homework 3")
getwd()
data <- read.table("hcmv.txt",header=TRUE)
spacing <- matrix(0, nrow = 295, ncol = 1)

#--------------------------------------------------spacing between consecutive palindromes----------------------------------------------------------------
for (i in 1:147)
{
  spacing[i,1] = data[i*2+1,1] - data[i*2,1]
}
lamda = 1/(mean(spacing))
hist(spacing, breaks = 30, col = 'lightgreen', xlim = c(0,6000), ylim = c(0, 100), main = "Distribution of spacings between consecutive palindromes")
#perform chi-square test
chisqtable <- function(spacing, b, rowname){
  
  count.int <- table(cut(spacing, breaks = b, include.lowest = TRUE))
  count.vector <- as.vector(count.int)
  count.range <- max(count.vector) - min(count.vector) + 1
  
  lambda.est <- 1/mean(spacing)
  
  table <- matrix(0, 12, 2)
  for (i in 1:8){
    table[i, 1] <- count.vector[i]
    table[i,2] <- 147*(pgamma(b[i+1], 2, lamda) - pgamma(b[i], 2, lamda))
  }
  colnames(table) <- c("Observed counts","Expected counts")
  rownames(table) <- rowname
  table <- as.table(table)
  return (table)
}

#------------------------------------------- GAMMA --------------------------------
b <- c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000)
rowname = c("0-500", "500-1000", "1000-1500", "1500-2000", "2000-2500", "2500-3000", "3000-3500", "3500-4000", "4000-4500", "4500-5000", "5000-5500", "5500-6000")
table1 = chisqtable(spacing, b, rowname)
table1

newTable1 <- matrix(rep(NA, 7*2), 7, 2)
newTable1[1:7,] <- table1[1:7, 1:2]
colnames(newTable1) <- c("Observed counts","Expected counts")
rowname = c("0-500", "500-1000", "1000-1500", "1500-2000", "2000-2500", "2500-3000", "3000-6000")
rownames(newTable1) <- rowname
newTable1

stats = sum((newTable1[,1] - newTable1[,2])^2/newTable1[,2])
p_value1 = pchisq(stats, df = 6, lower.tail = FALSE)
p_value1

residuals <- (newTable1[,1] - newTable1[,2])/sqrt(newTable1[,2])
plot(residuals, type = 'h', ylab = "standardized residuals", xlab = "intervals",xaxt = "n", main = "Residual plot")
axis(1, at=c(1:7), labels = c("1st", "2ed", "3rd", "4th", "5th", "6th", "7th"))










