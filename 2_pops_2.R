# I made most of the chunks into indpendent functions so that we can run them on as many
# populations as we want

start.pop <- 50
LF <- 0.3
percent.breed <- 0.5
carrying.capacity <- 2000
baseAttack <- c(.1, .1, .1, .1)
n.off <- 4

s1=c(1,.1,.1,.1)
s2=c(.1,1,.1,.1)
s3=c(.1,.1,1,.1)
s4=c(.1,.1,.1,1)

sim <- rbind(s1,s2,s3,s4)

T1=c(1,1,1,1)
T2=c(1,1,1,1)
T3=c(1,1,1,1)
T4=c(1,1,1,1)

hand <- rbind(T1,T2,T3,T4)

# make the starting matrices - we'll make the linked allele later

geno1 <- matrix(rbinom(start.pop*6, 1, (1/3)), ncol=6)
colnames(geno1) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

geno2 <- matrix(rbinom(start.pop*6, 1, (1/3)), ncol=6)
colnames(geno2) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")


# set recombination frequency, select which individuals will recombine

recombination=rbinom(start.pop, 1, LF)

# do the recombination

linked1 <- matrix(NA, nrow=start.pop, ncol=2)

	for(i in 1:start.pop){
	if(recombination[i]==0){linked1[i,] <- geno1[,3:4][i,]} else
		linked1[i,]<- geno1[,3:4][i,c(2,1)]
	}
	
linked2 <- matrix(NA, nrow=start.pop, ncol=2)

	for(i in 1:start.pop){
	if(recombination[i]==0){linked2[i,] <- geno2[,3:4][i,]} else
		linked2[i,] <- geno2[,3:4][i,c(2,1)]
	}


# the function for getting phenotypes from genotypes

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
	
geno1 <- cbind(geno1[,1:4], linked1, geno1[,5:6])
g1ph <- phenotype(rowSums(cbind(geno1[,1:2], geno1[,3:4]*3)))
geno1 <- cbind(g1ph, geno1[,1:4], linked1, geno1[,5:6])
colnames(geno1) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")

geno2 <- cbind(geno2[,1:4], linked2, geno2[,5:6])
g2ph <- phenotype(rowSums(cbind(geno2[,1:2], geno2[,3:4]*3)))
geno2 <- cbind(g2ph, geno2[,1:4], linked2, geno2[,5:6])
colnames(geno2) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")



# breeding

make.off <- function(n.off, mat, start.pop, percent.breed){
	
lucky <- sample(start.pop, percent.breed*start.pop)

pairs <- mat[lucky,]

pair1 <- pairs[1:(nrow(pairs)/2),]
pair2 <- pairs[(1+nrow(pairs)/2):nrow(pairs),]

bands.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
red.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
linked.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
neutral.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))


for(i in 1:nrow(pair1)){
	which.bands <- rbinom(n.off, 1, 0.5)+1
	bands.off1[,i] <- pair1[i,1:2][which.bands]
	which.allele <- rbinom(n.off, 1, 0.5)+1
	red.off1[,i] <- pair1[i,3:4][which.allele]
	linked.off1[,i] <- pair1[i,5:6][which.allele]
	which.neu <- rbinom(n.off, 1, 0.5)+1
	neutral.off1[,i] <- pair1[i,7:8][which.neu]
}

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


# negative frequency dependent selection

NFDS <- function(pgmat, base.attack, similarity, handling){

pt <- c(sum(pgmat[,1]==1), sum(pgmat[,1]==2), sum(pgmat[,1]==3), sum(pgmat[,1]==4))
	
	pheno1 <- pgmat[pgmat[,1]==1,]
	pheno2 <- pgmat[pgmat[,1]==2,]
	pheno3 <- pgmat[pgmat[,1]==3,]
	pheno4 <- pgmat[pgmat[,1]==4,]

denom1=matrix(NA, nrow=4, ncol=4)

for(k in 1:4){
	for (j in 1:4){
		denom1[j,k]=base.attack[k]*pt[k]*(1+similarity[k,j]*handling[k,j]*base.attack[j]*pt[j])
	}
}

denom=sum(denom1)

f1=matrix(NA, nrow=4, ncol=4)
for(l in 1:4){
	for(m in 1:4){
	f1[l,m]=base.attack[l]*pt[l]*similarity[l,m]*base.attack[m]*pt[m]
	}
	
}
f=colSums(f1)

surv1=round(pt-pt*(f/denom))
surv=(abs(surv1)+surv1)/2


both=ifelse(pt<surv, pt, surv)

phenolist=list(pheno1, pheno2, pheno3,pheno4)

phenolist2=list()

phenosub=c()

for(q in 1:4){
	if(both[q]>1){phenolist2[[q]] <- phenolist[[q]][1:both[q],]}
	else if(both[q]==1){phenolist2[[q]] <- phenolist[[q]]}
	else if(both[q]==0){phenolist2[[q]] <- phenosub}
	else{phenolist2[[q]] <- phenosub}
}


# now we have a matrix of individuals that survived the morph-specific
# predation

next.gen.2 <- do.call(rbind, phenolist2)

return(next.gen.2)	
}

# normal LV selection

LV <- function(NFmat, carrying.capacity, percent.breed, n.off){
rate.inc <- percent.breed*n.off

nt <- nrow(NFmat)
	
threshold <- nt*exp(rate.inc*(1-nt/carrying.capacity))
	
if(threshold > nrow(NFmat)){
	rand <- sample(nt)
	next.gen.1 <- NFmat[rand,]
	next.gen <- next.gen.1	
	}else{
	rand1 <- sample(nt)
	next.gen.1 <- NFmat[rand1,]
	next.gen <- next.gen.1[1:threshold,]
}

return(next.gen)
}


# make two alleles worth of genotypes - don't differentiate sexes - these are the first elements in a list

# this is for later, to get the average difference in allele frequency between the two populations

freqDiffs <- function(list){
	a1 <- colMeans(list[[1]])
	m1 <- cbind(mean(a1[2], a1[3]), mean(a1[4], a1[5]), mean(a1[6], a1[7]), mean(a1[8], a1[9]))
	a2 <- colMeans(list[[2]])
	m2 <- cbind(mean(a2[2], a2[3]), mean(a2[4], a2[5]), mean(a2[6], a2[7]), mean(a2[8], a2[9]))
	# I include an absolute value because we care about the magnitude of the distance, not the 
	# sign
	diff <- abs(m1-m2)
}



####################################			
# test for loop ####################
####################################

pops <- list()
pops[[1]] <- list(geno1, geno2)


for(i in 1:n.gen){

g1 <- pops[[i]][[1]][,2:9]
g2 <- pops[[i]][[2]][,2:9]

# exchange migrants
n.mig <- round(nrow(g1)*percent.migrate)

geno1m <- rbind(g2[1:n.mig,], g1[(n.mig+1):start.pop,])

geno2m <- rbind(g1[1:n.mig,], g2[(n.mig+1):start.pop,])



off1 <- make.off(4, geno1m, start.pop, percent.breed)
off2 <- make.off(4, geno2m, start.pop, percent.breed)

# make phenotypes

g1 <- rbind(geno1m, off1)
pheno1 <- phenotype(rowSums(cbind(g1[,1:2], g1[,3:4]*3)))
pg1 <- cbind(pheno1, g1)
order1 <- order(pg1[,1])
pg1 <- pg1[order1,]

g2 <- rbind(geno2m, off2)
pheno2 <- phenotype(rowSums(cbind(g2[,1:2], g2[,3:4]*3)))
pg2 <- cbind(pheno2, g2)
order2 <- order(pg2[,1])
pg2 <- pg2[order2,]

pt1 <- c(sum(pg1[,1]==1), sum(pg1[,1]==2), sum(pg1[,1]==3), sum(pg1[,1]==4))
	
	pheno1.1 <- pg1[pg1[,1]==1,]
	pheno1.2 <- pg1[pg1[,1]==2,]
	pheno1.3 <- pg1[pg1[,1]==3,]
	pheno1.4 <- pg2[pg2[,1]==4,]

pt2 <- c(sum(pg2[,1]==1), sum(pg2[,1]==2), sum(pg2[,1]==3), sum(pg2[,1]==4))
	
	pheno2.1 <- pg2[pg2[,1]==1,]
	pheno2.2 <- pg2[pg2[,1]==2,]
	pheno2.3 <- pg2[pg2[,1]==3,]
	pheno2.4 <- pg2[pg2[,1]==4,]

##############################################
# NFDS #######################################
##############################################

# this needs more thought - should frequencies in one population affect what the predator sees?
# we'll probably need two separate functions for that

NF1 <- NFDS(pg1, baseAttack, sim, hand)

NF2 <- NFDS(pg2, baseAttack, sim, hand)


# randomize

NF1 <- NF1[sample(nrow(NF1)),]
NF2 <- NF2[sample(nrow(NF2)),]

# normal selection



fin1 <- LV(NF1, carrying.capacity, percent.breed, n.off)
fin2 <- LV(NF2, carrying.capacity, percent.breed, n.off)

fin <- list(fin1, fin2)
# output this final pop to a list and pull it back to start over

pops[[i+1]] <- fin

}



allele.freq <- function(list){
	a1 <- colSums(list[[1]])/nrow(list[[1]])
	af1 <- c(mean(a1[2], a1[3]), mean(a1[4], a1[5]), mean(a1[6], a1[7]), mean(a1[8], a1[9]))
	a2 <- colSums(list[[2]])/nrow(list[[2]])
	af2 <- c(mean(a2[2], a2[3]), mean(a2[4], a2[5]), mean(a2[6], a2[7]), mean(a2[8], a2[9]))
	return(list(af1, af2))
}


###################################################################
# bands vs. red plot ##############################################
###################################################################


allele.freq.br <- function(list){
	a1 <- colSums(list[[1]])/nrow(list[[1]])
	af1 <- c(mean(a1[2], a1[3]), mean(a1[4], a1[5]))
	a2 <- colSums(list[[2]])/nrow(list[[2]])
	af2 <- c(mean(a2[2], a2[3]), mean(a2[4], a2[5]))
	return(list(af1, af2))
}


br.freq <- lapply(pops, allele.freq.br)

plot(x=seq(0,1,by=0.1), y=seq(0,1,by=0.1), type="n", xlab="bands frequency", ylab="red frequency")

bands.x <- c()
red.y <- c()

for(i in 1:length(pops)){
	bands.x[i] <- br.freq[[i]][[1]][1]
	red.y[i] <- br.freq[[i]][[1]][2]
}

bands.x2 <- c()
red.y2 <- c()

for(i in 1:length(pops)){
	bands.x2[i] <- br.freq[[i]][[2]][1]
	red.y2[i] <- br.freq[[i]][[2]][2]
}

points(bands.x, red.y, col='red')

points(bands.x2, red.y2, pch=15)


###################################
# function ########################
###################################

# set the parameters. The purpose of this function is to compare average
# differences in allele frequencies between the two populations, so n.gen
# should be >50 to get a decent average

percent.breed <- 0.5
carrying.capacity <- 100
start.pop <- 50
n.gen <- 50

# the function - takes a two element vector of percent migrating and recomb. frequency
# everything else is set. This is so we can feed it a wide range of parameter values 
# quickly and easily

 migLD <- function(vec){

# get the starting genotypes - this needs to be inside the function because
# we will do multiple iterations later - so we need independent starting populations 
# for each run of the simulation
#vec=c(0.1,0.1)
geno1 <- matrix(rbinom(start.pop*6, 1, (2/3)), ncol=6)
colnames(geno1) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

geno2 <- matrix(rbinom(start.pop*6, 1, (1/3)), ncol=6)
colnames(geno2) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

# do the recombination

recombination1 <- rbinom(start.pop, 1, vec[2])
recombination2 <- rbinom(start.pop, 1, vec[2])

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
#i=1
for(i in 1:n.gen){

g1 <- pops[[i]][[1]][,2:9]
g2 <- pops[[i]][[2]][,2:9]

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
	
	pheno1.1 <- pg1[pg1[,1]==1,]
	pheno1.2 <- pg1[pg1[,1]==2,]
	pheno1.3 <- pg1[pg1[,1]==3,]
	pheno1.4 <- pg1[pg1[,1]==4,]

pt2 <- c(sum(pg2[,1]==1), sum(pg2[,1]==2), sum(pg2[,1]==3), sum(pg2[,1]==4))
	
	pheno2.1 <- pg2[pg2[,1]==1,]
	pheno2.2 <- pg2[pg2[,1]==2,]
	pheno2.3 <- pg2[pg2[,1]==3,]
	pheno2.4 <- pg2[pg2[,1]==4,]

##############################################
# NFDS #######################################
##############################################

# this needs more thought - should frequencies in one population affect what the predator sees?
# we'll probably need two separate functions for that

NF1 <- NFDS(pg1, baseAttack, sim, hand)

NF2 <- NFDS(pg2, baseAttack, sim, hand)


# randomize

NF1 <- NF1[sample(nrow(NF1)),]
NF2 <- NF2[sample(nrow(NF2)),]

# normal selection

fin1 <- LV(NF1, carrying.capacity, percent.breed, n.off)
fin2 <- LV(NF2, carrying.capacity, percent.breed, n.off)

# make sure they recombine again 
r1 <- rbinom(nrow(fin1), 1, vec[2])
r2 <- rbinom(nrow(fin2), 1, vec[2])

l1 <- matrix(NA, nrow=nrow(fin1), ncol=2)

	for(k in 1:nrow(fin1)){
	if(r1[k]==0){l1[k,] <- fin1[,6:7][k,]} else
		l1[k,]<- fin1[,6:7][k,c(2,1)]
	}
	
l2 <- matrix(NA, nrow=nrow(fin2), ncol=2)

	for(k in 1:nrow(fin2)){
	if(r2[k]==0){l2[k,] <- fin2[,6:7][k,]} else
		l2[k,] <- fin2[,6:7][k,c(2,1)]
	}

FIN1 <- cbind(fin1[,2:5], l1, fin1[,8:9])
FINPH1 <- phenotype(rowSums(cbind(FIN1[,1:2], FIN1[,3:4]*3)))
fin.1 <- cbind(FINPH1, FIN1)
colnames(fin.1) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")

FIN2 <- cbind(fin2[,2:5], l2, fin2[,8:9])
FINPH2 <- phenotype(rowSums(cbind(FIN2[,1:2], FIN2[,3:4]*3)))
fin.2 <- cbind(FINPH2, FIN2)
colnames(fin.2) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")


fin <- list(fin.1, fin.2)
# output this final pop to a list and pull it back to start over

pops[[i+1]] <- fin

}

# once the list is made, we find the difference in allele frequency between the 
# two populations at each generation

diffs <- lapply(pops, freqDiffs)

fMat <- matrix(unlist(diffs), ncol=4, byrow=T)

return(fMat)
}

# decide on the ranges of the migration % and recomb frequency we want to test

pm1 <- seq(0, 0.1, by=0.01)
rf1 <- seq(0, 0.5, by=0.1)

# now repeat the complete first vector the same number of times as the length of second vector

pm <- rep(pm1, length(rf1))

# repeat each element of the second vector the same number of times as the length of the first vector
rf <- rep(rf1, each=length(pm1))

# now make a matrix of the two vectors bound together - this way each value of migration
# is paired with each value of recomb frequency to test the entire range of parameters

test <- cbind(pm, rf)

# make each row of the matrix into an element in a list - just makes the apply easier

ltest <- list()

for(i in 1:nrow(test)){
	ltest[[i]] <- test[i,]
}

# now iterate the function and lapply multiple times to get averages of behavior 
# of the model at each paramter value

repLD <- list()

for(j in 1:10){

# ltest is the same for each iteration, but re-running migLD will get us different
# starting points and progression through the generations
	
af <- lapply(ltest, migLD)

# get colmeans for each run - the columns are the loci, the rows are the
# difference in allele frequencies between population 1 and population 2
# at each generation, so taking colmeans gets you the mean difference between
# populations at that locus across mutliple generations

means <- lapply(af, function(mat){x <- colMeans(mat); return(x)})

# this gets the list of means into a matrix, which is output into a list

repLD[[j]] <- matrix(unlist(means), ncol=4, byrow=T)
	
}

# get the mean of the means across runs - each row is an allele 
# bands, red, linked, unlinked
# each row is a set of parameter values

mean <- Reduce('+', repLD, repLD[[1]])/10

# take the mean values for the "band" locus, make them into a matrix
# with values of pm along the rows and values of rf for the columns

xbandMeans <- matrix(mean[,1], ncol=length(rf1))

xredMeans <- matrix(mean[,2], ncol=length(rf1))

xlMeans <- matrix(mean[,3], ncol=length(rf1))

xulMeans <- matrix(mean[,4], ncol=length(rf1))


# plots!

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))

persp(pm1, rf1, xbandMeans,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.25))

persp(pm1, rf1, xredMeans,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.25))

persp(pm1, rf1, xlMeans, theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.5))

persp(pm1, rf1, xulMeans,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.5))



#########################################################
# predators see both populations ########################
#########################################################

NFDS2 <- function(fullpg, poppg, base.attack, similarity, handling){
#fullpg <- pg
#poppg <- pg1
#base.attack <- baseAttack
#similarity <- sim
#handling <- hand

pt <- c(sum(fullpg[,1]==1), sum(fullpg[,1]==2), sum(fullpg[,1]==3), sum(fullpg[,1]==4))
	
	pheno1 <- fullpg[fullpg[,1]==1,]
	pheno2 <- fullpg[fullpg[,1]==2,]
	pheno3 <- fullpg[fullpg[,1]==3,]
	pheno4 <- fullpg[fullpg[,1]==4,]

denom1 <- matrix(NA, nrow=4, ncol=4)

for(k in 1:4){
	for (j in 1:4){
		denom1[j,k]=base.attack[k]*pt[k]*(1+similarity[k,j]*handling[k,j]*base.attack[j]*pt[j])
	}
}

denom <- sum(denom1)

f1 <- matrix(NA, nrow=4, ncol=4)
for(l in 1:4){
	for(m in 1:4){
	f1[l,m]=base.attack[l]*pt[l]*similarity[l,m]*base.attack[m]*pt[m]
	}
	
}
f <- colSums(f1)

# apply to actual matrices

pt2 <- c(sum(poppg[,1]==1), sum(poppg[,1]==2), sum(poppg[,1]==3), sum(poppg[,1]==4))
	
	pheno1.2 <- poppg[poppg[,1]==1,]
	pheno2.2 <- poppg[poppg[,1]==2,]
	pheno3.2 <- poppg[poppg[,1]==3,]
	pheno4.2 <- poppg[poppg[,1]==4,]

sur2 <- round(pt2-pt2*(f/denom))

surv2 <- (abs(sur2)+sur2)/2

both <- ifelse(pt2<surv2, pt2, surv2)

phenolist=list(pheno1.2, pheno2.2, pheno3.2, pheno4.2)

phenolist2=list()

phenosub=c()

for(q in 1:4){
	if(both[q]>1){phenolist2[[q]] <- phenolist[[q]][1:both[q],]}
	else if(both[q]==1){phenolist2[[q]] <- phenolist[[q]]}
	else if(both[q]==0){phenolist2[[q]] <- phenosub}
	else{phenolist2[[q]] <- phenosub}
}


# now we have a matrix of individuals that survived the morph-specific
# predation

next.gen.2 <- do.call(rbind, phenolist2)

return(next.gen.2)	
}


migLD2 <- function(vec){

# get the starting genotypes - this needs to be inside the function because
# we will do multiple iterations later - so we need independent starting populations 
# for each run of the simulation
#vec=c(0.1,0.1)
geno1 <- matrix(rbinom(start.pop*6, 1, (1/2)), ncol=6)
colnames(geno1) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

geno2 <- matrix(rbinom(start.pop*6, 1, (1/2)), ncol=6)
colnames(geno2) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

# do the recombination

recombination1 <- rbinom(start.pop, 1, vec[2])
recombination2 <- rbinom(start.pop, 1, vec[2])

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
#i=1
for(i in 1:n.gen){

g1 <- pops[[i]][[1]][,2:9]
g2 <- pops[[i]][[2]][,2:9]

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

G <- rbind(G1, G2)
ph <- phenotype(rowSums(cbind(G[,1:2], G[,3:4]*3)))
pg <- cbind(ph, G)
order <- order(pg[,1])
pg <- pg[order,]

pt1 <- c(sum(pg1[,1]==1), sum(pg1[,1]==2), sum(pg1[,1]==3), sum(pg1[,1]==4))
	
	pheno1.1 <- pg1[pg1[,1]==1,]
	pheno1.2 <- pg1[pg1[,1]==2,]
	pheno1.3 <- pg1[pg1[,1]==3,]
	pheno1.4 <- pg1[pg1[,1]==4,]

pt2 <- c(sum(pg2[,1]==1), sum(pg2[,1]==2), sum(pg2[,1]==3), sum(pg2[,1]==4))
	
	pheno2.1 <- pg2[pg2[,1]==1,]
	pheno2.2 <- pg2[pg2[,1]==2,]
	pheno2.3 <- pg2[pg2[,1]==3,]
	pheno2.4 <- pg2[pg2[,1]==4,]

# NFDS #######################################


# this needs more thought - should frequencies in one population affect what the predator sees?
# we'll probably need two separate functions for that

NF1 <- NFDS2(pg,pg1, baseAttack, sim, hand)

NF2 <- NFDS2(pg,pg2, baseAttack, sim, hand)

# randomize

NF1 <- NF1[sample(nrow(NF1)),]
NF2 <- NF2[sample(nrow(NF2)),]

# normal selection

fin1 <- LV(NF1, carrying.capacity, percent.breed, n.off)
fin2 <- LV(NF2, carrying.capacity, percent.breed, n.off)

# make sure they recombine again 
r1 <- rbinom(nrow(fin1), 1, vec[2])
r2 <- rbinom(nrow(fin2), 1, vec[2])

l1 <- matrix(NA, nrow=nrow(fin1), ncol=2)

	for(k in 1:nrow(fin1)){
	if(r1[k]==0){l1[k,] <- fin1[,6:7][k,]} else
		l1[k,]<- fin1[,6:7][k,c(2,1)]
	}
	
l2 <- matrix(NA, nrow=nrow(fin2), ncol=2)

	for(k in 1:nrow(fin2)){
	if(r2[k]==0){l2[k,] <- fin2[,6:7][k,]} else
		l2[k,] <- fin2[,6:7][k,c(2,1)]
	}

FIN1 <- cbind(fin1[,2:5], l1, fin1[,8:9])
FINPH1 <- phenotype(rowSums(cbind(FIN1[,1:2], FIN1[,3:4]*3)))
fin.1 <- cbind(FINPH1, FIN1)
colnames(fin.1) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")

FIN2 <- cbind(fin2[,2:5], l2, fin2[,8:9])
FINPH2 <- phenotype(rowSums(cbind(FIN2[,1:2], FIN2[,3:4]*3)))
fin.2 <- cbind(FINPH2, FIN2)
colnames(fin.2) <- c("phenotype","bands1", "bands2", "red1", "red2", "linked1", "linked2","neutral1", "neutral2")


fin <- list(fin.1, fin.2)
# output this final pop to a list and pull it back to start over

pops[[i+1]] <- fin

}

# once the list is made, we find the difference in allele frequency between the 
# two populations at each generation

diffs <- lapply(pops, freqDiffs)

fMat <- matrix(unlist(diffs), ncol=4, byrow=T)

return(fMat)
}

repLD2 <- list()

for(j in 1:5){

# ltest is the same for each iteration, but re-running migLD will get us different
# starting points and progression through the generations
	
af <- lapply(ltest, migLD2)

# get colmeans for each run - the columns are the loci, the rows are the
# difference in allele frequencies between population 1 and population 2
# at each generation, so taking colmeans gets you the mean difference between
# populations at that locus across mutliple generations

means <- lapply(af, function(mat){x <- colMeans(mat); return(x)})

# this gets the list of means into a matrix, which is output into a list

repLD2[[j]] <- matrix(unlist(means), ncol=4, byrow=T)
	
}

# get the mean of the means across runs - each row is an allele 
# bands, red, linked, unlinked
# each row is a set of parameter values

mean <- Reduce('+', repLD2, repLD2[[1]])/10

# take the mean values for the "band" locus, make them into a matrix
# with values of pm along the rows and values of rf for the columns

xbandMeans <- matrix(mean[,1], ncol=length(rf1))

xredMeans <- matrix(mean[,2], ncol=length(rf1))

xlMeans <- matrix(mean[,3], ncol=length(rf1))

xulMeans <- matrix(mean[,4], ncol=length(rf1))


# plots!

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))

persp(pm1, rf1, xbandMeans,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.25))

persp(pm1, rf1, xredMeans,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.25))

persp(pm1, rf1, xlMeans, theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.5))

persp(pm1, rf1, xulMeans,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.5))



