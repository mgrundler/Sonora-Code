# I made most of the chunks into indpendent functions so that we can run them on as many
# populations as we want
rm(list=ls())

# this is the recombination frequency for the linked/nuetral allele
LF <- 0.3
# in the absence of frequency dependent selection, the percentage of the time that the predator will attack each morph upon an encounter
baseAttack <- c(.5, .5, .5, .5)
# number of offspring per female per breeding (can be interpreted as the # that survive to adulthood)
n.off <- 4

# "switching similarity" matrix (see Van Leeuwen et al. 2013 for explanation) - similarity of each morph to all other morph is 10% here
s1=c(1,.1,.1,.1)
s2=c(.1,1,.1,.1)
s3=c(.1,.1,1,.1)
s4=c(.1,.1,.1,1)

sim <- rbind(s1,s2,s3,s4)

# this function takes in the coding alleles and returns phenotype, accounting for dominance

phenotype=function(offspring.phenotype){
offspring.phenotype1=ifelse(offspring.phenotype==0, 1, offspring.phenotype)

offspring.phenotype2=ifelse(offspring.phenotype==1, 2, offspring.phenotype1)

offspring.phenotype3=ifelse(offspring.phenotype==2, 2, offspring.phenotype2)

offspring.phenotype4=ifelse(offspring.phenotype==3, 3, offspring.phenotype3)

offspring.phenotype5=ifelse(offspring.phenotype==4, 4, offspring.phenotype4)

offspring.phenotype6=ifelse(offspring.phenotype==5, 4, offspring.phenotype5)

offspring.phenotype7=ifelse(offspring.phenotype==6, 3, offspring.phenotype6)

offspring.phenotype8=ifelse(offspring.phenotype==7, 4, offspring.phenotype7)

offspring.phenotype9=ifelse(offspring.phenotype==8, 4, offspring.phenotype8)

offspring.phenotype10=ifelse(offspring.phenotype==9, 1, offspring.phenotype9)

offspring.phenotype11=ifelse(offspring.phenotype==10, 2, offspring.phenotype10)

offspring.phenotype12=ifelse(offspring.phenotype==11, 3, offspring.phenotype11)

offspring.phenotype13=ifelse(offspring.phenotype==12, 4, offspring.phenotype12)

return(offspring.phenotype13)
}

# the function to get offspring from parental genotypes
# inputs are number of offspring, a matrix of parental genotypes, the number of parental individuals, and the percentage of the parental individuals that breed

make.off <- function(n.off, mat, start.pop, percent.breed){

# select the individuals that will breed - a random selection	
lucky <- sample(start.pop, percent.breed*start.pop)

# reduce the matrix to just the lucky breeders
pairs <- mat[lucky,]

# divide the lucky ones randomly into "male" and "female" matrices, which have equal numbers of individuals
pair1 <- pairs[1:(nrow(pairs)/2),]
pair2 <- pairs[(1+nrow(pairs)/2):nrow(pairs),]

# make the offspring - create matrices in which each column represents a parent, and each row in that column is a random draw from that parent's diploid genotype

# make the matrices
bands.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
red.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
linked.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
neutral.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))

# loop over all of the first set of parental individuals
for(i in 1:nrow(pair1)){
	# randomly select either the first or second allele - with a length of the number of offspring
	which.bands <- rbinom(n.off, 1, 0.5)+1
	# fill the corresponding column of the matrix
	bands.off1[,i] <- pair1[i,1:2][which.bands]
	# do the same for the second coding allele
	which.allele <- rbinom(n.off, 1, 0.5)+1
	red.off1[,i] <- pair1[i,3:4][which.allele]
	# use the same selection for the linked, noncoding allele
	linked.off1[,i] <- pair1[i,5:6][which.allele]
	# now for a neutral allele
	which.neu <- rbinom(n.off, 1, 0.5)+1
	neutral.off1[,i] <- pair1[i,7:8][which.neu]
}

# do the same with the second set of parents - this means that the mated pairs are determined by the position of the parents in the two matrices. This is fine as long as the matrices are randomized.

bands.off2 <- matrix(nrow=n.off, ncol=nrow(pair1))
red.off2 <- matrix(nrow=n.off, ncol=nrow(pair1))
linked.off2 <- matrix(nrow=n.off, ncol=nrow(pair1))
neutral.off2 <- matrix(nrow=n.off, ncol=nrow(pair1))

for(i in 1:nrow(pair1)){
	which.bands <- rbinom(n.off, 1, 0.5)+1
	bands.off2[,i] <- pair2[i,1:2][which.bands]
	which.allele <- rbinom(n.off, 1, 0.5)+1
	red.off2[,i] <- pair2[i,3:4][which.allele]
	linked.off2[,i] <- pair2[i,5:6][which.allele]
	which.neu <- rbinom(n.off, 1, 0.5)+1
	neutral.off2[,i] <- pair2[i,7:8][which.neu]
}

offspring <- cbind(as.vector(bands.off1), as.vector(bands.off2), as.vector(red.off1),
as.vector(red.off2), as.vector(linked.off1), as.vector(linked.off2), as.vector(neutral.off1),
as.vector(neutral.off2))

return(offspring)

}

# negative frequency dependent selection - inputs are a matrix of phenotypes+genotypes, the base attack rate, and the similarity matrix

NFDS <- function(pgmat, base.attack, similarity){

# this gives the frequencies of each morph in the input matrix
pt <- c(sum(pgmat[,1]==1), sum(pgmat[,1]==2), sum(pgmat[,1]==3), sum(pgmat[,1]==4))
	
	pheno1 <- pgmat[pgmat[,1]==1,]
	pheno2 <- pgmat[pgmat[,1]==2,]
	pheno3 <- pgmat[pgmat[,1]==3,]
	pheno4 <- pgmat[pgmat[,1]==4,]

# this is the rate at which each morph would be eaten in the absence of NFDS
tildeN <- base.attack*pt

# make a matrix of the number of each of of the base attack rates multiplied by the switching similarity matrix
	
switches <- rbind(tildeN*similarity[1,], tildeN*similarity[2,],tildeN*similarity[3,],
tildeN*similarity[4,])	

# total number of individuals taken per morph
totals <- tildeN*rowSums(switches)

# find the proportion of individuals taken in each morph
proportions <- totals/sum(totals)}


###################################
# function ########################
###################################

# parameters - number of breeders
percent.breed <- 0.5
carrying.capacity <- 500
start.pop <- 200
n.gen <- 300

# the function - takes a two element vector of percent migrating and the percentage of the total number of mortalities per generation that are due to NFDS

migLD <- function(vec){

# get the starting genotypes - this needs to be inside the function because
# we will do multiple iterations later - so we need independent starting populations 
# for each run of the simulation

geno1 <- matrix(rbinom(start.pop*6, 1, (2/3)), ncol=6)
colnames(geno1) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

geno2 <- matrix(rbinom(start.pop*6, 1, (1/3)), ncol=6)
colnames(geno2) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

# do the recombination

recombination1 <- rbinom(start.pop, 1, LF)
recombination2 <- rbinom(start.pop, 1, LF)

linked1 <- matrix(NA, nrow=start.pop, ncol=2)

	for(i in 1:start.pop){
	if(recombination1[i]==0){linked1[i,] <- geno1[,3:4][i,]} else
		linked1[i,]<- geno1[,3:4][i,c(2,1)]
	}
	
linked2 <- matrix(NA, nrow=start.pop, ncol=2)

	for(i in 1:start.pop){
	if(recombination2[i]==0){linked2[i,] <- geno2[,3:4][i,]} else
		linked2[i,] <- geno2[,3:4][i,c(2,1)]
	}

geno1 <- cbind(geno1[,1:4], linked1, geno1[,5:6])
g1ph <- phenotype(rowSums(cbind(geno1[,1:2], geno1[,3:4]*3)))
geno1 <- cbind(g1ph, geno1[,1:4], linked1, geno1[,5:6])
colnames(geno1) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")

geno2 <- cbind(geno2[,1:4], linked2, geno2[,5:6])
g2ph <- phenotype(rowSums(cbind(geno2[,1:2], geno2[,3:4]*3)))
geno2 <- cbind(g2ph, geno2[,1:4], linked2, geno2[,5:6])
colnames(geno2) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")

pops <- list()

pops[[1]] <- list(geno1, geno2)


# now we do the for loop to fill the list 

for(z in 1:n.gen){

g1 <- pops[[z]][[1]][,2:9]
g2 <- pops[[z]][[2]][,2:9]

# exchange migrants
n.mig <- round(nrow(g1)*vec[1])

if(n.mig==0){
	geno1m <- g1
	geno2m <- g2
}else{
geno1m <- rbind(g2[1:n.mig,], g1[(n.mig+1):nrow(g1),])

geno2m <- rbind(g1[1:n.mig,], g2[(n.mig+1):nrow(g2),])
}

# make offspring in each populaiton
off1 <- make.off(4, geno1m, nrow(geno1m), percent.breed)
off2 <- make.off(4, geno2m, nrow(geno2m), percent.breed)

# make phenotypes for parental generation (all of them, including the non-breeders), plus the offspring

G1 <- rbind(geno1m, off1)
# get phenotypes
pheno1 <- phenotype(rowSums(cbind(G1[,1:2], G1[,3:4]*3)))
pg1 <- cbind(pheno1, G1)
# group the phenotypes together
order1 <- order(pg1[,1])
pg1 <- pg1[order1,]

# repeat for second population
G2 <- rbind(geno2m, off2)
pheno2 <- phenotype(rowSums(cbind(G2[,1:2], G2[,3:4]*3)))
pg2 <- cbind(pheno2, G2)
order2 <- order(pg2[,1])
pg2 <- pg2[order2,]

# get the number of each morph
pt1 <- c(sum(pg1[,1]==1), sum(pg1[,1]==2), sum(pg1[,1]==3), sum(pg1[,1]==4))
	
phen <- 1:4

# make a list of matrices of genotypes for each morph
ptlist1 <- list()

for(i in 1:4){
	if(pt1[i]==0){ptlist1[[i]]=c()}
	else{ptlist1[[i]]=pg1[pg1[,1]==phen[i],]}
}
	
pt2 <- c(sum(pg2[,1]==1), sum(pg2[,1]==2), sum(pg2[,1]==3), sum(pg2[,1]==4))
	
ptlist2 <- list()

for(i in 1:4){
	if(pt2[i]==0){ptlist2[[i]]=c()}
	else{ptlist2[[i]]=pg2[pg2[,1]==phen[i],]}
}


# NFDS #######################################
# bind all genotpes into an input matrix for the NFDS funciton 
fds1 <- do.call(rbind, ptlist1)
fds2 <- do.call(rbind, ptlist2)

# find the percentage of each morph that will be removed by NFDS - sums to 1
# the predator can only see the population in which they are predating
NF1 <- NFDS(fds1, baseAttack, sim)

NF2 <- NFDS(fds2, baseAttack, sim)

# non-NFDS mortality ##########################

# find the number of individuals will survive - this is set by the logistic mortality function

# population 1
threshold1 <- carrying.capacity*nrow(fds1)*exp(n.off*percent.breed)/(carrying.capacity + nrow(fds1)*(exp(n.off*percent.breed)-1))

# population 2
threshold2 <- carrying.capacity*nrow(fds2)*exp(n.off*percent.breed)/(carrying.capacity + nrow(fds2)*(exp(n.off*percent.breed)-1))

# no predation happens if the threshold is greater than the current population size
# if there is some mortality, this will return the number of each morph that will be removed due to NFDS - this is set by the the total mortality (nrow(fds1)-threshold1) time the percentage of killed individuals exposed to NFDS - vec[2], times a vector that divides those mortalities into a per-morph basis

if(nrow(fds1)>threshold1){
	n.pred1 <- round(((nrow(fds1)-threshold1)*vec[2])*NF1)
}else{n.pred1 <- c(0,0,0,0)}

if(nrow(fds2)>threshold2){
	n.pred2 <- round(((nrow(fds2)-threshold2)*vec[2])*NF2)
}else{n.pred2 <- c(0,0,0,0)}

# no remove the NFDS predated invididuals from the per-morph matrices

predated1 <- list()

for(i in 1:4){
	if(pt1[i] > (n.pred1[i]+2)){predated1[[i]] <- ptlist1[[i]][1:(pt1[i]-n.pred1[i]),]}else{predated1[[i]] <- c()}
}

predated2 <- list()

for(i in 1:4){
	if(pt2[i]>(n.pred2[i]+2)){predated2[[i]] <- ptlist2[[i]][1:(pt2[i]-n.pred2[i]),]}else{
		predated2[[i]] <- c()}
}

# combine the post-predation matrices and randomize

fin1 <- do.call(rbind, predated1)
fin1 <- fin1[sample(nrow(fin1)),]
fin2 <- do.call(rbind, predated2)
fin2 <- fin2[sample(nrow(fin2)),]

# clip to the length of threshold, as long as threshold is less than the current population size

if(nrow(fin1)>threshold1){
	fin1 <- fin1[1:threshold1,]
} else if(nrow(fin1)>1){fin1 <- fin1}else{fin1 <- c()}


if(nrow(fin2)>threshold2){
	fin2 <- fin2[1:threshold2,]
} else if(nrow(fin2)>1){fin2 <- fin2}else{fin2<- c()}


fin <- list(fin1, fin2)

pops[[z+1]] <- fin

}

return(pops)
}

# decide on the ranges of the migration % and recomb frequency we want to test

pm1 <- seq(0, 0.1, by=0.01)
nfds1 <- seq(0, 1, by=0.2)

# now repeat the complete first vector the same number of times as the length of second vector

pm <- rep(pm1, length(nfds1))

# repeat each element of the second vector the same number of times as the length of the first vector
nfds <- rep(nfds1, each=length(pm1))

# now make a matrix of the two vectors bound together - this way each value of migration
# is paired with each value of recomb frequency to test the entire range of parameters

test <- cbind(pm, nfds)

# make each row of the matrix into an element in a list - just makes the apply easier

ltest <- list()

for(i in 1:nrow(test)){
	ltest[[i]] <- test[i,]
}

# now iterate the function and lapply multiple times to get averages of behavior 
# of the model at each paramter value

sonoraX1 <- list()

for(q in 1:length(ltest)){
sonoraX1[[q]] <- migLD(ltest[[q]])	
}

# get the allele frequencies for a single population at a single time point

alleleFreqs <- function(list){
	a1 <- colMeans(list[[1]])
	m1 <- cbind(mean(a1[2], a1[3]), mean(a1[4], a1[5]), mean(a1[6], a1[7]), mean(a1[8], a1[9]))
	a2 <- colMeans(list[[2]])
	m2 <- cbind(mean(a2[2], a2[3]), mean(a2[4], a2[5]), mean(a2[6], a2[7]), mean(a2[8], a2[9]))
return(list(m1, m2))
}

# get the allele frequencies over time
afTime <- function(timelist){
y <- lapply(timelist, alleleFreqs)

# pull out the allele freqs for each population to make a matrix

pop1FreqMat <- matrix(NA, ncol=4, nrow=n.gen+1)
colnames(pop1FreqMat) <- c("black", "red","ul", "linked")

for(i in 1:n.gen+1){
	pop1FreqMat[i,] <- y[[i]][[1]]
}


pop2FreqMat <- matrix(NA, ncol=4, nrow=n.gen+1)
colnames(pop2FreqMat) <- c("black", "red", "ul", "linked")

for(i in 1:n.gen+1){
	pop2FreqMat[i,] <- y[[i]][[2]]
}

return(list(pop1FreqMat, pop2FreqMat))
}

time1 <- afTime(x[[10]])

plot(seq(0,1, by=1/n.gen), type="n")
lines(time1[[1]][,2], col="red")
lines(time1[[1]][,1], lty=3, col="red")
lines(time1[[2]][,2], col="blue")
lines(time1[[2]][,1], lty=3, col="blue")

# get difference in allele frequencies over time between the two pops

freqDiffs <- function(biglist){
d <- lapply(biglist, afTime)

e <- lapply(d, function(list){return(abs(list[[1]]-list[[2]]))})

f <- lapply(e, function(list){return(colMeans(list[10:nrow(list),]))})

g <- matrix(unlist(f), ncol=4, byrow=T)
colnames(g) <- c("black", "red","ul", "linked")

return(g)
	
}

sonorafDs <- freqDiffs(sonoraX1)

redMat <- matrix(sonorafDs[,1], nrow=length(pm1), byrow=F)

blackMat <- matrix(sonorafDs[,2], nrow=length(pm1), byrow=F)



# plots!

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))

persp(pm1, rf1, redMat,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.5))

persp(pm1, rf1, blackMat,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.5))


#########################################################
# predators see both populations ########################
#########################################################

migLDboth <- function(vec){

geno1 <- matrix(rbinom(start.pop*6, 1, (2/3)), ncol=6)
colnames(geno1) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

geno2 <- matrix(rbinom(start.pop*6, 1, (1/3)), ncol=6)
colnames(geno2) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

# do the recombination

recombination1 <- rbinom(start.pop, 1, LF)
recombination2 <- rbinom(start.pop, 1, LF)

linked1 <- matrix(NA, nrow=start.pop, ncol=2)

	for(i in 1:start.pop){
	if(recombination1[i]==0){linked1[i,] <- geno1[,3:4][i,]} else
		linked1[i,]<- geno1[,3:4][i,c(2,1)]
	}
	
linked2 <- matrix(NA, nrow=start.pop, ncol=2)

	for(i in 1:start.pop){
	if(recombination2[i]==0){linked2[i,] <- geno2[,3:4][i,]} else
		linked2[i,] <- geno2[,3:4][i,c(2,1)]
	}

geno1 <- cbind(geno1[,1:4], linked1, geno1[,5:6])
g1ph <- phenotype(rowSums(cbind(geno1[,1:2], geno1[,3:4]*3)))
geno1 <- cbind(g1ph, geno1[,1:4], linked1, geno1[,5:6])
colnames(geno1) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")

geno2 <- cbind(geno2[,1:4], linked2, geno2[,5:6])
g2ph <- phenotype(rowSums(cbind(geno2[,1:2], geno2[,3:4]*3)))
geno2 <- cbind(g2ph, geno2[,1:4], linked2, geno2[,5:6])
colnames(geno2) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")

pops <- list()

pops[[1]] <- list(geno1, geno2)


for(z in 1:n.gen){

g1 <- pops[[z]][[1]][,2:9]
g2 <- pops[[z]][[2]][,2:9]

# exchange migrants
n.mig <- round(nrow(g1)*vec[1])

if(n.mig==0){
	geno1m <- g1
	geno2m <- g2
}else{
geno1m <- rbind(g2[1:n.mig,], g1[(n.mig+1):nrow(g1),])

geno2m <- rbind(g1[1:n.mig,], g2[(n.mig+1):nrow(g2),])
}

off1 <- make.off(4, geno1m, nrow(geno1m), percent.breed)
off2 <- make.off(4, geno2m, nrow(geno2m), percent.breed)

# make phenotypes

G1 <- rbind(geno1m, off1)
pheno1 <- phenotype(rowSums(cbind(G1[,1:2], G1[,3:4]*3)))
pg1 <- cbind(pheno1, G1)
order1 <- order(pg1[,1])
pg1 <- pg1[order1,]

G2 <- rbind(geno2m, off2)
pheno2 <- phenotype(rowSums(cbind(G2[,1:2], G2[,3:4]*3)))
pg2 <- cbind(pheno2, G2)
order2 <- order(pg2[,1])
pg2 <- pg2[order2,]

pt1 <- c(sum(pg1[,1]==1), sum(pg1[,1]==2), sum(pg1[,1]==3), sum(pg1[,1]==4))
	
phen <- 1:4
ptlist1 <- list()

for(i in 1:4){
	if(pt1[i]==0){ptlist1[[i]]=c()}
	else{ptlist1[[i]]=pg1[pg1[,1]==phen[i],]}
}
	
pt2 <- c(sum(pg2[,1]==1), sum(pg2[,1]==2), sum(pg2[,1]==3), sum(pg2[,1]==4))
	
ptlist2 <- list()

for(i in 1:4){
	if(pt2[i]==0){ptlist2[[i]]=c()}
	else{ptlist2[[i]]=pg2[pg2[,1]==phen[i],]}
}

##############################################
# NFDS #######################################
##############################################
fds1 <- do.call(rbind, ptlist1)
fds2 <- do.call(rbind, ptlist2)

fds <- rbind(fds1, fds2)

# the predator only works in a single population at a given generation, but "sees" both populations

NF1 <- NFDS(fds, baseAttack, sim)

NF2 <- NFDS(fds, baseAttack, sim)

threshold1 <- carrying.capacity*nrow(fds1)*exp(n.off*percent.breed)/(carrying.capacity + nrow(fds1)*(exp(n.off*percent.breed)-1))

threshold2 <- carrying.capacity*nrow(fds2)*exp(n.off*percent.breed)/(carrying.capacity + nrow(fds2)*(exp(n.off*percent.breed)-1))

if(nrow(fds1)>threshold1){
	n.pred1 <- round(((nrow(fds1)-threshold1)*vec[2])*NF1)
}else{n.pred1 <- c(0,0,0,0)}

if(nrow(fds2)>threshold2){
	n.pred2 <- round(((nrow(fds2)-threshold2)*vec[2])*NF2)
}else{n.pred2 <- c(0,0,0,0)}

predated1 <- list()

for(i in 1:4){
	if(pt1[i] > (n.pred1[i]+2)){predated1[[i]] <- ptlist1[[i]][1:(pt1[i]-n.pred1[i]),]}else{predated1[[i]] <- c()}
}

predated2 <- list()

for(i in 1:4){
	if(pt2[i]>(n.pred2[i]+2)){predated2[[i]] <- ptlist2[[i]][1:(pt2[i]-n.pred2[i]),]}else{
		predated2[[i]] <- c()}
}

fin1 <- do.call(rbind, predated1)
fin1 <- fin1[sample(nrow(fin1)),]
fin2 <- do.call(rbind, predated2)
fin2 <- fin2[sample(nrow(fin2)),]

if(nrow(fin1)>threshold1){
	fin1 <- fin1[1:threshold1,]
} else if(nrow(fin1)>1){fin1 <- fin1}else{fin1 <- c()}


if(nrow(fin2)>threshold2){
	fin2 <- fin2[1:threshold2,]
} else if(nrow(fin2)>1){fin2 <- fin2}else{fin2<- c()}


fin <- list(fin1, fin2)

pops[[z+1]] <- fin

}


return(pops)
}


sonoraX2 <- list()

for(q in 1:length(ltest)){
sonoraX2[[q]] <- migLDboth(ltest[[q]])	
}

timeSonora1 <- afTime(sonoraX1[[45]])
timeSonora2 <- afTime(sonoraX2[[45]])

par(mar=c(3,3,3,3))
par(mfrow=c(2,1))
plot(seq(0,1, by=1/n.gen), type="n", main="one", xlab="generation")
lines(timeSonora1[[1]][,1], col="red")
lines(timeSonora1[[1]][,2], lty=2, col="red")
lines(timeSonora1[[1]][,4], lty=3, col="red")
lines(timeSonora1[[2]][,1], col="blue")
lines(timeSonora1[[2]][,2], lty=2, col="blue")
lines(timeSonora1[[2]][,4], lty=3, col="blue")
legend(5, 0.3, lty=c(1,2,3), legend=c("bands", "red", "unlinked"), col="red")

plot(seq(0,1, by=1/n.gen), type="n", main="both", xlab="generation")
lines(timeSonora2[[1]][,1], col="red")
lines(timeSonora2[[1]][,2], lty=2, col="red")
lines(timeSonora2[[1]][,4], lty=3, col="red")
lines(timeSonora2[[2]][,1], col="blue")
lines(timeSonora2[[2]][,2], lty=2, col="blue")
lines(timeSonora2[[2]][,4], lty=3, col="blue")
legend(5, 0.3, lty=c(1,2,3), legend=c("bands", "red", "unlinked"), col="red")



sonora1fDs <- freqDiffs(sonoraX1)
sonora2fDs <- freqDiffs(sonoraX2)

redMat1 <- matrix(sonora1fDs[,1], nrow=length(pm1), byrow=F)

ulMat1 <- matrix(sonora1fDs[,4], nrow=length(pm1), byrow=F)

redMat2 <- matrix(sonora2fDs[,1], nrow=length(pm1), byrow=F)

ulMat2 <- matrix(sonora2fDs[,4], nrow=length(pm1), byrow=F)

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))

persp(pm1, nfds1, redMat1,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,1), main="coding, one", xlab="migration", ylab="NFDS", zlab="allele freq difference")

persp(pm1, nfds1, ulMat1 ,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,1),main="unlinked, one", xlab="migration", ylab="NFDS", zlab="allele freq difference")


persp(pm1, nfds1, redMat2,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,1),main="coding, both", xlab="migration", ylab="NFDS", zlab="allele freq difference")

persp(pm1, nfds1, ulMat2,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,1),main="unlinked, both", xlab="migration", ylab="NFDS", zlab="allele freq difference")


############################################################
# code to parse outcomes of multiple runs of same variables#
############################################################
# the function to take in 
parseMorphs <- function(listFreqs){
	# get just the morphs for each population
	intList <- lapply(listFreqs, function(list2){return(list(list2[[1]][,1], list2[[2]][,1]))})
	# see which morphs are present
	freqList <- lapply(intList, function(list3){return(list(table(list3[[1]]), table(list3[[2]])))})
	# get the lengths of the list - will be 4 for case 1 (all 4 morphs present in both populations), 2 for case 2 (2 morphs in pop 1, the other two in pop 2), 1 for case 3 (each population fixed for a different morph)

	lengthList <- lapply(freqList, function(list4){return(list(length(list4[[1]]), length(list4[[2]])))})
	
	# get generations 200-300

	lengthListCut <- lengthList[200:300]
	# separate pop 1 from pop 2	
	length1 <- c()
	length2 <- c()
	
	for(i in 1:length(lengthListCut)){
		length1[i] <- lengthListCut[[i]][[1]]
		length2[i] <- lengthListCut[[i]][[2]]
	}
	# find the average # of morphs for each population across generations 200-500
	mean1 <- mean(length1)
	mean2 <- mean(length2)
		
	# parse the results - if the mean number of morphs across 300 generations is greater than 3, then throughout most of the time, both populations will have all four morphs. If the mean number is less than 1.5, then in most cases both populations were fixed for a single morph (feel free to change the cutoff points as makes sense to you). 
	morphResult <- c()
	if(mean1 > 3 & mean2 > 3){morphResult[1]=1} else if(mean1 < 1.5 & mean2 > 2.5 | mean1 > 2.5 & mean2 < 1.5 ){morphResult[1] = 3}else{morphResult[1]=2}
	
	return(morphResult)
}

# re-running the simulation with the same parameters 10 times
bpSeeBoth <- list()

for(j in 1:100){
	bpSeeBoth[[j]] <- migLDboth(c(0.01,0.4))
}

# getting the outcomes of the simulations
bpBoth <- lapply(bpSeeBoth, parseMorphs)
# see how many simulations lead to each outcome
tabBoth <- table(unlist(bpBoth))
# make a barplot - this can get much fancier!
barplot(tabBoth)

# do the same with the predator only seeing one pop at a time
bpSeeOne <- list()

for(j in 3:10){
	bpSeeOne[[j]] <- migLD(c(0.01,0.95))
}

barPlotOne <- lapply(bpSeeOne, parseMorphs)
tabone <- table(unlist(barPlotOne))

barplot(rbind(c(0,22,8), c(30,0,0)), beside=T)



pm1 <- seq(0, 0.1, by=0.01)
nfds1 <- seq(0, 1, by=0.05)

# now repeat the complete first vector the same number of times as the length of second vector

pm <- rep(pm1, length(nfds1))

# repeat each element of the second vector the same number of times as the length of the first vector
nfds <- rep(nfds1, each=length(pm1))

# now make a matrix of the two vectors bound together - this way each value of migration
# is paired with each value of recomb frequency to test the entire range of parameters

test <- cbind(pm, nfds)

# make each row of the matrix into an element in a list - just makes the apply easier

ltest <- list()

for(i in 1:nrow(test)){
	ltest[[i]] <- test[i,]
}


outcome <- list()

for(b in 1:length(ltest)){
bpSeeBoth <- list()

for(j in 1:10){
	bpSeeBoth[[j]] <- migLD(ltest[[b]])
}

# getting the outcomes of the simulations
bpBoth <- lapply(bpSeeBoth, parseMorphs)
# see how many simulations lead to each outcome
outcome[[b]] <- mean(unlist(bpBoth))
}


outcome1 <- unlist(outcome)

outMat <- matrix(outcome1, nrow=length(pm1))
outMat <- outMat -1

write.table(outMat, "~/Desktop/Sonora/sonora1pop.csv")
matrix(pm, nrow=11)
matrix(nfds, nrow=11)

sonora1pop <- as.matrix(read.csv("~/Desktop/Sonora/sonora1pop.csv"))
colnames(sonora1pop) <- nfds1
rownames(sonora1pop) <- pm1
sonora2pop <- as.matrix(read.csv("~/Desktop/Sonora/sonora2pop.csv"))
colnames(sonora2pop) <- nfds1
rownames(sonora2pop) <- pm1

hel1pop <- as.matrix(read.csv("~/Dropbox/AmNat_2015/heliconius1pop.csv"))+1
hel2pop <- as.matrix(read.csv("~/Dropbox/AmNat_2015/heliconius2pop.csv"))

cep2pop <- as.matrix(read.csv("~/Dropbox/AmNat_2015/Cepaea_Outcomes_csv.files/Cepaea_2_pop.csv"))
cep1pop <- as.matrix(read.csv("~/Dropbox/AmNat_2015/Cepaea_Outcomes_csv.files/Cepaea1pop.csv"))

library(lattice)
library(gridExtra)
library(png)

heatCol=rainbow(200)

par(mar=c(1,1,1,1))


mim <- readPNG(file.choose())
poly <- readPNG(file.choose())
redfix <- readPNG(file.choose())
yelfix <- readPNG(file.choose())
yelmim <- readPNG(file.choose())


sonora1 <- wireframe(sonora1pop, drape=T, zlim=c(0,4), xlab="m", ylab="s", pretty=F, zlab="", col.regions=heatCol, at=seq(from=0, to=4, by=0.2), main="local", default.scales=list(arrows=T))

sonora2 <- wireframe(sonora2pop, drape=T, zlim=c(0,4), xlab="m", ylab="s", pretty=F, zlab="",col.regions=heatCol, at=seq(from=0, to=4, by=0.2), main="regional",default.scales=list(arrows=T))

hel1 <- wireframe(hel1pop, drape=T, zlim=c(0,4), xlab="m", ylab="s", pretty=F, zlab="",col.regions=heatCol, at=seq(from=0, to=4, by=0.2),default.scales=list(arrows=T))

hel2 <- wireframe(hel2pop, drape=T, zlim=c(0,4), xlab="m", ylab="s", pretty=F, zlab="",col.regions=heatCol, at=seq(from=0, to=4, by=0.2),default.scales=list(arrows=T))

cep1 <- wireframe(cep1pop, drape=T, zlim=c(0,4), xlab="m", ylab="s", pretty=F, zlab="",col.regions=heatCol, at=seq(from=0, to=4, by=0.2),default.scales=list(arrows=T))

cep2 <- wireframe(cep2pop, drape=T, zlim=c(0,4), xlab="m", ylab="s", pretty=F, zlab="",col.regions=heatCol, at=seq(from=0, to=4, by=0.2),default.scales=list(arrows=T))

quartz(height=11, width=9.5)
trellis.par.set("axis.line",list(col="black",lty=1,lwd=1))
grid.arrange(sonora1, sonora2, cep1, cep2, hel1, hel2,ncol=2)
grid.text("Sonora", x=0.04, y=0.84, rot=90)
grid.text("Cepaea", x=0.04, y=0.5, rot=90)
grid.text("Heliconius", x=0.04, y=0.22, rot=90)


grid.raster(poly, width=0.02, x=0.07, y=0.815)
grid.raster(poly, width=0.02, x=0.09, y=0.815)

grid.raster(redfix, width=0.02, x=0.07, y=0.835)
grid.raster(yelfix, width=0.02, x=0.09, y=0.835)

grid.raster(mim, width=0.02, x=0.07, y=0.855)
grid.raster(yelfix, width=0.02, x=0.09, y=0.855)

grid.raster(mim, width=0.02, x=0.07, y=0.875)
grid.raster(yelmim, width=0.02, x=0.09, y=0.875)

grid.raster(mim, width=0.02, x=0.07, y=0.895)
grid.raster(mim, width=0.02, x=0.09, y=0.895)

########################################################
grid.raster(poly, width=0.02, x=0.57, y=0.815)
grid.raster(poly, width=0.02, x=0.59, y=0.815)

grid.raster(redfix, width=0.02, x=0.57, y=0.835)
grid.raster(yelfix, width=0.02, x=0.59, y=0.835)

grid.raster(mim, width=0.02, x=0.57, y=0.855)
grid.raster(yelfix, width=0.02, x=0.59, y=0.855)

grid.raster(mim, width=0.02, x=0.57, y=0.875)
grid.raster(yelmim, width=0.02, x=0.59, y=0.875)

grid.raster(mim, width=0.02, x=0.57, y=0.895)
grid.raster(mim, width=0.02, x=0.59, y=0.895)
 
###########################################

grid.raster(poly, width=0.02, x=0.57, y=0.49)
grid.raster(poly, width=0.02, x=0.59, y=0.49)

grid.raster(redfix, width=0.02, x=0.57, y=0.51)
grid.raster(yelfix, width=0.02, x=0.59, y=0.51)

grid.raster(mim, width=0.02, x=0.57, y=0.53)
grid.raster(yelfix, width=0.02, x=0.59, y=0.53)

grid.raster(mim, width=0.02, x=0.57, y=0.55)
grid.raster(yelmim, width=0.02, x=0.59, y=0.55)

grid.raster(mim, width=0.02, x=0.57, y=0.57)
grid.raster(mim, width=0.02, x=0.59, y=0.57)
 
############################################

grid.raster(poly, width=0.02, x=0.07, y=0.49)
grid.raster(poly, width=0.02, x=0.09, y=0.49)

grid.raster(redfix, width=0.02, x=0.07, y=0.51)
grid.raster(yelfix, width=0.02, x=0.09, y=0.51)

grid.raster(mim, width=0.02, x=0.07, y=0.53)
grid.raster(yelfix, width=0.02, x=0.09, y=0.53)

grid.raster(mim, width=0.02, x=0.07, y=0.55)
grid.raster(yelmim, width=0.02, x=0.09, y=0.55)

grid.raster(mim, width=0.02, x=0.07, y=0.57)
grid.raster(mim, width=0.02, x=0.09, y=0.57)

##############################################


grid.raster(poly, width=0.02, x=0.57, y=0.165)
grid.raster(poly, width=0.02, x=0.59, y=0.165)

grid.raster(redfix, width=0.02, x=0.57, y=0.185)
grid.raster(yelfix, width=0.02, x=0.59, y=0.185)

grid.raster(mim, width=0.02, x=0.57, y=0.205)
grid.raster(yelfix, width=0.02, x=0.59, y=0.205)

grid.raster(mim, width=0.02, x=0.57, y=0.225)
grid.raster(yelmim, width=0.02, x=0.59, y=0.225)

grid.raster(mim, width=0.02, x=0.57, y=0.245)
grid.raster(mim, width=0.02, x=0.59, y=0.245)
 
############################################

grid.raster(poly, width=0.02, x=0.07, y=0.165)
grid.raster(poly, width=0.02, x=0.09, y=0.165)

grid.raster(redfix, width=0.02, x=0.07, y=0.185)
grid.raster(yelfix, width=0.02, x=0.09, y=0.185)

grid.raster(mim, width=0.02, x=0.07, y=0.205)
grid.raster(yelfix, width=0.02, x=0.09, y=0.205)

grid.raster(mim, width=0.02, x=0.07, y=0.225)
grid.raster(yelmim, width=0.02, x=0.09, y=0.225)

grid.raster(mim, width=0.02, x=0.07, y=0.245)
grid.raster(mim, width=0.02, x=0.09, y=0.245)



