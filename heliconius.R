# getting rid of any residual named objects 
rm(list=ls())

ba1 <- c(0.9, 0.01, 0.9, 0.9)
ba2 <- c(0.9, 0.0, 0.01, 0.9)

LF <- 0.3

n.mig <- 0.1

matchSurv <- 0.9
nomatchSurv <- 0.1

rep.dist <- 0.8
rep.nondist <- 0.5
 
s1=c(0,.9,.9,.9)
s2=c(.9,0,.9,.9)
s3=c(.9,.9,0,.9)
s4=c(.9,.9,.9,0)

sim <- rbind(s1,s2,s3,s4)

# the offspring function - this is specific to Heliconius
# the "Mat" is the genotypes of all of the individuals in the population
# for Cepea you will need to apply this twice, to each species
# make sure the matrices are randomly ordered going into this function

make.off <- function(n.off, mat){

# split the population into "males" and "females"

pair1 <- mat[1:(nrow(mat)/2),]
pair2 <- mat[(1+nrow(mat)/2):nrow(mat),]

# make the matrices for the gametes from the "males"

bar.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
spot.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
sprec.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
mphrec.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
ul.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
linked.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))

for(i in 1:nrow(pair1)){
# choose the bar alleles
	which.bar <- rbinom(n.off, 1, 0.5)+1
# put them in the matrix
	bar.off1[,i] <- pair1[i,2:3][which.bar]
# since bar and spot are linked, choose the same alleles are bar
	spot.off1[,i] <- pair1[i,4:5][which.bar]
# this is the species recognition allele - probably not relevant to Cepea
	which.sprec <- rbinom(n.off, 1, 0.5)+1
	sprec.off1[,i] <- pair1[i,6:7][which.sprec]
# moprh recognition allele
	which.mphrec <- rbinom(n.off, 1, 0.5)+1
	mphrec.off1[,i] <- pair1[i,8:9][which.mphrec]
	which.ul <- rbinom(n.off, 1, 0.5)+1
	ul.off1[,i] <- pair1[i,10:11][which.ul]		
	linked.off1[,i] <- pair1[i,12:13][which.bar]
}

# make them into vectors

bar.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
spot.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
sprec.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
mphrec.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
ul.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
linked.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))

# same for the "females"

for(i in 1:nrow(pair2)){
	which.bar <- rbinom(n.off, 1, 0.5)+1
	bar.off2[,i] <- pair2[i,2:3][which.bar]
	spot.off2[,i] <- pair2[i,4:5][which.bar]
	which.sprec <- rbinom(n.off, 1, 0.5)+1
	sprec.off2[,i] <- pair2[i,6:7][which.sprec]
	which.mphrec <- rbinom(n.off, 1, 0.5)+1
	mphrec.off2[,i] <- pair2[i,8:9][which.mphrec]
	which.ul <- rbinom(n.off, 1, 0.5)+1
	ul.off2[,i] <- pair2[i,10:11][which.ul]		
	linked.off2[,i] <- pair2[i,12:13][which.bar]
}

# make the offspring

offspring <- cbind(as.vector(bar.off1), as.vector(bar.off2), as.vector(spot.off1),
as.vector(spot.off2), as.vector(sprec.off1), as.vector(sprec.off2), as.vector(mphrec.off1),
as.vector(mphrec.off2), as.vector(ul.off1), as.vector(ul.off2))

# make the recombination 
recombination <- rbinom(nrow(offspring), 1, LF)

linked1 <- matrix(NA, nrow=nrow(offspring), ncol=2)

for(i in 1:nrow(offspring)){
	if(recombination[i]==0){linked1[i,] <- offspring[,4:5][i,]} else
		linked1[i,]<- offspring[,4:5][i,c(2,1)]
}

recombSpot <- rbinom(nrow(offspring), 1, LF)

# do the recombination between bar and spot
spotPre <- cbind(as.vector(spot.off1), as.vector(spot.off2))

spot1 <- matrix(NA, nrow=nrow(spotPre), ncol=2)

for(i in 1:nrow(spotPre)){
	if(recombSpot[i]==0){spot1[i,] <- spotPre[,1:2][i,]} else
	spot1[i,]<- spotPre[,1:2][i,c(2,1)]
}

# make the offspring a second time with the recombined spot and linked alleles

offspring <- cbind(offspring[,1:2], spot1, offspring[,5:10], linked1)
colnames(offspring) <- c("bar1", "bar2", "spot1", "spot2", "sprec1", "sprec2", "mphrec1", "mphrec2","ul1", "ul2", "linked1", "linked2")

return(offspring)
}

phenotype <- function(vector, mat){
genos <- mat[,vector]

offspring.phenotype <- rowSums(cbind(genos[,1:2], genos[,3:4]*3))

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


# this does the parsing out of the potential parents based on whether they can 
# recognize morphs and species - not relevant to Cepea

#mat <- geno1m
get.heliconius.offspring <- function(mat, n.off){

# four phenotypes from the two color alleles

pheno <- phenotype(1:4,mat)

# breeding
# split species distinguishing from not
pg <- cbind(pheno, mat)

# first divide by species recognition

sprec <- rowSums(cbind(mat[,5], mat[,6]))
sprec[sprec==2]=1

# do the if/else thing to account for times that there
# are none of one genotype or another

rec <- list()

if(sum(sprec==1)/length(sprec)==0){rec[[1]]=pg[1:round(length(sprec)*rep.nondist),]
	} else if(sum(sprec==1)/length(sprec)==1){rec[[1]]=pg[1:round(length(sprec)*rep.dist),]
	} else if(sum(sprec==1)*rep.dist < 1 & sum(sprec==0)*rep.nondist >1){rec[[1]] <- pg[1:round(length(sprec)*rep.nondist),]} else if(sum(sprec==1)*rep.dist > 1 & sum(sprec==0)*rep.nondist <1){rec[[1]]=pg[1:round(length(sprec)*rep.dist),]}else{
	r <- pg[sprec==1,]
	rec[[1]] <- r[1:round(nrow(r)*rep.dist),]
	nr <- pg[sprec==0,]
	rec[[2]] <- nr[1:round(nrow(nr)*rep.nondist),]
}

# now do the same for the morph recognizers

recLucky <- do.call(rbind, rec)

# break up by morph rec

rr <- rowSums(cbind(recLucky[,8], recLucky[,9]))
rr[rr==2]=1

morphs <- unique(recLucky[,1])

morphBred <- list()

if(sum(rr==1)/length(rr)==0){
	morphBred[[1]] <- recLucky	
} else if(sum(rr==1)/length(rr)==1){
	morphBred[[1]] <- recLucky
} else{
	morphBred[[1]]=recLucky[rr==0,]
	
	recMorph <- recLucky[rr==1,]
	for(i in 1:length(morphs)){
		morphBred[[i+1]] = recMorph[recMorph[,1]==morphs[i],]
	}
}

percent.breed <- sum(unlist(lapply(morphBred, nrow)))/nrow(mat)

off <- list()	

for(i in 1:length(morphBred)){
	if(class(morphBred[[i]])=="numeric"){
		off[[i]] <- morphBred[[i]][2:13]
		} else if(nrow(morphBred[[i]]) < 6){
			off[[i]] <- morphBred[[i]][,2:13]
		} else {
			off[[i]] <- make.off(n.off, morphBred[[i]])
		}
}


offspring <- do.call(rbind, off)
offspring <- offspring[sample(nrow(offspring)),]

return(list(offspring, percent.breed))

}

# positive (or negative, depending on the "sim" matrix) frequency dependent selection function
# takes in a matrix of genotypes, base attack vector, and similarity matrix

PFDS <- function(pgmat, base.attack, similarity){

pt <- c(sum(pgmat[,1]==1), sum(pgmat[,1]==2), sum(pgmat[,1]==3), sum(pgmat[,1]==4))

# parses matrix into each phenotype	

	pheno1 <- pgmat[pgmat[,1]==1,]
	pheno2 <- pgmat[pgmat[,1]==2,]
	pheno3 <- pgmat[pgmat[,1]==3,]
	pheno4 <- pgmat[pgmat[,1]==4,]
	
# find how many would be eaten without frequency dependence

tildeN <- base.attack*pt

# make a matrix of the number of each of of the base attack rates multiplied by the switching similarity matrix
	
switches <- rbind(tildeN*similarity[1,], tildeN*similarity[2,],tildeN*similarity[3,],
tildeN*similarity[4,])	

# total number of individuals taken per morph
totals <- tildeN*rowSums(switches)

# find the proportion of individuals taken in each morph
proportions <- totals/sum(totals)
}

####################################
# I made two functions - one where the predator sees both populations, one where the predator sees one
####################################
# predator sees single pop #########
####################################


start.pop <- 1000
p1start <- 0.5
p2start <- 0.5
n.off <- 3
carrying.capacity <- 1000
ngen <- 100


#vec <- c(0.1, 0.1)
#vec <- ltest[[12]]

migHel1pop <- function(vec){

# get the genotypes started	
geno1 <- matrix(rbinom(start.pop*8, 1, p1start), ncol=8)

recombSpot <- rbinom(nrow(geno1), 1, LF)

# do the recombination
spot1 <- matrix(NA, nrow=nrow(geno1), ncol=2)

for(i in 1:start.pop){
	if(recombSpot[i]==0){spot1[i,] <- geno1[,1:2][i,]} else
	spot1[i,]<- geno1[,1:2][i,c(2,1)]
}

recombination <- rbinom(nrow(geno1), 1, LF)

linked1 <- matrix(NA, nrow=nrow(geno1), ncol=2)

for(i in 1:start.pop){
	if(recombination[i]==0){linked1[i,] <- geno1[,1:2][i,]} else
		linked1[i,]<- geno1[,1:2][i,c(2,1)]
}

pheno1 <- phenotype(1:4, geno1)

geno1 <- cbind(pheno1, geno1[,1:2], spot1, geno1[,3:8], linked1)
colnames(geno1) <- c("phenotype","bar1", "bar2", "spot1", "spot2", "sprec1", "sprec2", "mphrec1", "mphrec2","ul1", "ul2", "linked1", "linked2")

# second pop

geno2 <- matrix(rbinom(start.pop*8, 1, p2start), ncol=8)

recombSpot2 <- rbinom(nrow(geno2), 1, LF)

# do the recombination
spot1.2 <- matrix(NA, nrow=nrow(geno2), ncol=2)

for(i in 1:start.pop){
	if(recombSpot2[i]==0){spot1.2[i,] <- geno2[,1:2][i,]} else
	spot1.2[i,]<- geno2[,1:2][i,c(2,1)]
}

recombination2 <- rbinom(nrow(geno2), 1, LF)

linked1.2 <- matrix(NA, nrow=nrow(geno2), ncol=2)

for(i in 1:start.pop){
	if(recombination2[i]==0){linked1.2[i,] <- geno2[,1:2][i,]} else
		linked1.2[i,]<- geno2[,1:2][i,c(2,1)]
}

pheno2 <- phenotype(1:4, geno2)

geno2 <- cbind(pheno2, geno2[,1:2], spot1.2, geno2[,3:8], linked1.2)
colnames(geno2) <- c("phenotype","bar1", "bar2", "spot1", "spot2", "sprec1", "sprec2", "mphrec1", "mphrec2","ul1", "ul2", "linked1", "linked2")

# start for loop ###############################

pops <- list()

pops[[1]] <-list(geno1, geno2)

#z=9
for(z in 1:ngen){
		# make offspring
	g1 <- pops[[z]][[1]][,2:13]	
	g2 <- pops[[z]][[2]][,2:13]

	# migrate
	n.mig <- round(nrow(g1)*vec[1])

if(n.mig==0){
	geno1m <- g1
	geno2m <- g2
	}else{
	geno1m <- rbind(g2[1:n.mig,], g1[(n.mig+1):nrow(g1),])

	geno2m <- rbind(g1[1:n.mig,], g2[(n.mig+1):nrow(g2),])
}
	
	# make offspring - this accounts for who gets to reproduce according to morph/species recognition,
	# and the resulting offspring
	offspring1 <- get.heliconius.offspring(geno1m,n.off)
	offspring2 <- get.heliconius.offspring(geno2m,n.off)
	off1 <- offspring1[[1]]
	off2 <- offspring2[[1]]
	pb1 <- offspring1[[2]]
	pb2 <- offspring2[[2]]
	
	# make phenotypes - break them in lists and use the if/else to account for morphs that aren't in the population
	
	pheno1 <- phenotype(1:4, off1)
	pheno2 <- phenotype(1:4, off2)
	
	pg1 <- cbind(pheno1, off1)
	order1 <- order(pg1[,1])
	pg1 <- pg1[order1,]
	
pt1 <- c(sum(pg1[,1]==1), sum(pg1[,1]==2), sum(pg1[,1]==3), sum(pg1[,1]==4))

phen <- 1:4
ptlist1 <- list()

for(i in 1:4){
	if(pt1[i]==0){ptlist1[[i]]=c()}
	else{ptlist1[[i]]=pg1[pg1[,1]==phen[i],]}
}
	
	pg2 <- cbind(pheno2, off2)
	order2 <- order(pg2[,1])
	pg2 <- pg2[order2,]
	
pt2 <- c(sum(pg2[,1]==1), sum(pg2[,1]==2), sum(pg2[,1]==3), sum(pg2[,1]==4))

ptlist2 <- list()

for(i in 1:4){
	if(pt2[i]==0){ptlist2[[i]]=c()}
	else{ptlist2[[i]]=pg2[pg2[,1]==phen[i],]}
}

### frequency dependent selection ####################
fds1 <- do.call(rbind, ptlist1)
fds2 <- do.call(rbind, ptlist2)

# do the frequency dependent selection - you can put any matrix into the function

FDS1 <- PFDS(fds1, ba1, sim)

FDS2 <- PFDS(fds2, ba2, sim)
######################################################

# use the LV equation to decide how many will die

threshold1 <- round(nrow(fds1)*exp(n.off*pb1*(1-nrow(fds1)/carrying.capacity)))
threshold2 <- round(nrow(fds2)*exp(n.off*pb2*(1-nrow(fds2)/carrying.capacity)))

# if/else to deal with non-existent morphs - vec[2] decides what percentage of the total 
# mortality will be morph-dependent. You take that percentage and divide it up between 
# morphs to get the final number that will be removed

if(nrow(fds1)>threshold1){
	n.pred1 <- round(((nrow(fds1)-threshold1)*vec[2])*FDS1)
}else{n.pred1 <- c(0,0,0,0)}

if(nrow(fds2)>threshold2){
	n.pred2 <- round(((nrow(fds2)-threshold2)*vec[2])*FDS2)
}else{n.pred2 <- c(0,0,0,0)}

# remove the predated ones

pt1[is.na(pt1)] <- 0
n.pred1[is.na(n.pred1)] <- 0

pt2[is.na(pt2)] <- 0
n.pred2[is.na(n.pred2)] <- 0

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

# now get rid of the rest that will die, taking random individuals

LV1 <- round(threshold1-(nrow(fds1)-threshold1)*(vec[2]))
LV2 <- round(threshold2-(nrow(fds2)-threshold2)*(vec[2]))

if(nrow(fin1)>threshold1){
	fin1 <- fin1[1:LV1,]
} else if(nrow(fin1)>1){fin1 <- fin1}else{fin1 <- c()}


if(nrow(fin2)>threshold2){
	fin2 <- fin2[1:LV2,]
} else if(nrow(fin2)>1){fin2 <- fin2}else{fin2<- c()}


fin <- list(fin1, fin2)

pops[[z+1]] <- fin

}

# get the frequency differences between populations using this code
#diffs <- lapply(pops, freqDiffs)

#fMat <- matrix(unlist(diffs), ncol=6, byrow=T)

return(pops)

}

# this is all to make figures

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


r <- migHel1pop(ltest[[13]])

x <- list()


for(q in 1:length(ltest)){
x[[q]] <- migHel1pop(ltest[[q]])	
}

# we have a list of lists of lists - each element of the 36 main list is a list of n.gen matrices, each of which
# has two lists for the two populations. We want the change in allele frequency of each population over the generations

# make a function to pull allele frequencies from the 2-pop sub list for one generation

alleleFreqs <- function(list){
	a1 <- colMeans(list[[1]])
	m1 <- cbind(mean(a1[2], a1[3]), mean(a1[4], a1[5]), mean(a1[6], a1[7]), mean(a1[8], a1[9]), mean(a1[10], a1[11]), mean(a1[12], a1[13]))
	a2 <- colMeans(list[[2]])
	m2 <- cbind(mean(a2[2], a2[3]), mean(a2[4], a2[5]), mean(a2[6], a2[7]), mean(a2[8], a2[9]),mean(a2[10], a2[11]), mean(a2[12], a2[13]))
return(list(m1, m2))
}

# apply the allele freqs function to the full list for all generations for one set of input parameters


afTime <- function(timelist){
y <- lapply(timelist, alleleFreqs)

# pull out the allele freqs for each population to make a matrix

pop1FreqMat <- matrix(NA, ncol=6, nrow=ngen+1)
colnames(pop1FreqMat) <- c("bar", "spot", "sprec", "mphrec", "ul", "linked")

for(i in 1:ngen+1){
	pop1FreqMat[i,] <- y[[i]][[1]]
}


pop2FreqMat <- matrix(NA, ncol=6, nrow=ngen+1)
colnames(pop2FreqMat) <- c("bar", "spot", "sprec", "mphrec", "ul", "linked")

for(i in 1:ngen+1){
	pop2FreqMat[i,] <- y[[i]][[2]]
}

return(list(pop1FreqMat, pop2FreqMat))
}


time1 <- afTime(x[[62]])

plot(seq(0,1, by=1/ngen), type="n")
lines(time1[[1]][,5], col="red")
lines(time1[[1]][,1], lty=3, col="red")
lines(time1[[2]][,5], col="blue")
lines(time1[[2]][,1], lty=3, col="blue")

# could we do something with mean/sd across generations?

# 3-d plot

# get frequency differences in each generation, return average for each allele at each parameter value


freqDiffs <- function(biglist){
d <- lapply(biglist, afTime)

e <- lapply(d, function(list){return(abs(list[[1]]-list[[2]]))})

f <- lapply(e, function(list){return(colMeans(list[10:nrow(list),]))})


g <- matrix(unlist(f), ncol=6, byrow=T)
colnames(g) <- c("bar", "spot", "sprec", "mphrec", "ul", "linked")

return(g)
	
}

fDs <- freqDiffs(x)

# make a matrix of percent migrate vs. percent nfds for each allele
# the rows will be percent migrate (6 rows), the columns percent nfds (also 6)

barMat <- matrix(fDs[,1], nrow=length(pm1), byrow=F)

ulMat <- matrix(fDs[,5], nrow=length(pm1), byrow=F)



# plots!

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))

persp(pm1, nfds1, barMat,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.3))

persp(pm1, nfds1, ulMat,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.3))




#############################################
# two pops ##################################
#############################################


start.pop <- 1000
p1start <- 0.5
p2start <- 0.5
n.off <- 3
carrying.capacity <- 1000
ngen <- 100


#vec <- c(0.1, 0.1)
#vec <- ltest[[12]]

migHel2pop <- function(vec){

# get the genotypes started	
geno1 <- matrix(rbinom(start.pop*8, 1, p1start), ncol=8)

recombSpot <- rbinom(nrow(geno1), 1, LF)

# do the recombination
spot1 <- matrix(NA, nrow=nrow(geno1), ncol=2)

for(i in 1:start.pop){
	if(recombSpot[i]==0){spot1[i,] <- geno1[,1:2][i,]} else
	spot1[i,]<- geno1[,1:2][i,c(2,1)]
}

recombination <- rbinom(nrow(geno1), 1, LF)

linked1 <- matrix(NA, nrow=nrow(geno1), ncol=2)

for(i in 1:start.pop){
	if(recombination[i]==0){linked1[i,] <- geno1[,1:2][i,]} else
		linked1[i,]<- geno1[,1:2][i,c(2,1)]
}

pheno1 <- phenotype(1:4, geno1)

geno1 <- cbind(pheno1, geno1[,1:2], spot1, geno1[,3:8], linked1)
colnames(geno1) <- c("phenotype","bar1", "bar2", "spot1", "spot2", "sprec1", "sprec2", "mphrec1", "mphrec2","ul1", "ul2", "linked1", "linked2")

# second pop

geno2 <- matrix(rbinom(start.pop*8, 1, p2start), ncol=8)

recombSpot2 <- rbinom(nrow(geno2), 1, LF)

# do the recombination
spot1.2 <- matrix(NA, nrow=nrow(geno2), ncol=2)

for(i in 1:start.pop){
	if(recombSpot2[i]==0){spot1.2[i,] <- geno2[,1:2][i,]} else
	spot1.2[i,]<- geno2[,1:2][i,c(2,1)]
}

recombination2 <- rbinom(nrow(geno2), 1, LF)

linked1.2 <- matrix(NA, nrow=nrow(geno2), ncol=2)

for(i in 1:start.pop){
	if(recombination2[i]==0){linked1.2[i,] <- geno2[,1:2][i,]} else
		linked1.2[i,]<- geno2[,1:2][i,c(2,1)]
}

pheno2 <- phenotype(1:4, geno2)

geno2 <- cbind(pheno2, geno2[,1:2], spot1.2, geno2[,3:8], linked1.2)
colnames(geno2) <- c("phenotype","bar1", "bar2", "spot1", "spot2", "sprec1", "sprec2", "mphrec1", "mphrec2","ul1", "ul2", "linked1", "linked2")

# start for loop ###############################

pops <- list()

pops[[1]] <-list(geno1, geno2)

#z=9
for(z in 1:ngen){
		# make offspring
	g1 <- pops[[z]][[1]][,2:13]	
	g2 <- pops[[z]][[2]][,2:13]

	# migrate
	n.mig <- round(nrow(g1)*vec[1])

if(n.mig==0){
	geno1m <- g1
	geno2m <- g2
	}else{
	geno1m <- rbind(g2[1:n.mig,], g1[(n.mig+1):nrow(g1),])

	geno2m <- rbind(g1[1:n.mig,], g2[(n.mig+1):nrow(g2),])
}
	
	# make offspring - this accounts for who gets to reproduce according to morph/species recognition,
	# and the resulting offspring
	offspring1 <- get.heliconius.offspring(geno1m,n.off)
	offspring2 <- get.heliconius.offspring(geno2m,n.off)
	off1 <- offspring1[[1]]
	off2 <- offspring2[[1]]
	pb1 <- offspring1[[2]]
	pb2 <- offspring2[[2]]
	
	# make phenotypes - break them in lists and use the if/else to account for morphs that aren't in the population
	
	pheno1 <- phenotype(1:4, off1)
	pheno2 <- phenotype(1:4, off2)
	
	pg1 <- cbind(pheno1, off1)
	order1 <- order(pg1[,1])
	pg1 <- pg1[order1,]
	
pt1 <- c(sum(pg1[,1]==1), sum(pg1[,1]==2), sum(pg1[,1]==3), sum(pg1[,1]==4))

phen <- 1:4
ptlist1 <- list()

for(i in 1:4){
	if(pt1[i]==0){ptlist1[[i]]=c()}
	else{ptlist1[[i]]=pg1[pg1[,1]==phen[i],]}
}
	
	pg2 <- cbind(pheno2, off2)
	order2 <- order(pg2[,1])
	pg2 <- pg2[order2,]
	
pt2 <- c(sum(pg2[,1]==1), sum(pg2[,1]==2), sum(pg2[,1]==3), sum(pg2[,1]==4))

ptlist2 <- list()

for(i in 1:4){
	if(pt2[i]==0){ptlist2[[i]]=c()}
	else{ptlist2[[i]]=pg2[pg2[,1]==phen[i],]}
}

### frequency dependent selection ####################
fds1 <- do.call(rbind, ptlist1)
fds2 <- do.call(rbind, ptlist2)

fds <- rbind(fds1, fds2)

# do the frequency dependent selection - you can put any matrix into the function

FDS1 <- PFDS(fds, ba1, sim)

FDS2 <- PFDS(fds, ba2, sim)
######################################################

# use the LV equation to decide how many will die

threshold1 <- round(nrow(fds1)*exp(n.off*pb1*(1-nrow(fds1)/carrying.capacity)))
threshold2 <- round(nrow(fds2)*exp(n.off*pb2*(1-nrow(fds2)/carrying.capacity)))

# if/else to deal with non-existent morphs - vec[2] decides what percentage of the total 
# mortality will be morph-dependent. You take that percentage and divide it up between 
# morphs to get the final number that will be removed

if(nrow(fds1)>threshold1){
	n.pred1 <- round(((nrow(fds1)-threshold1)*vec[2])*FDS1)
}else{n.pred1 <- c(0,0,0,0)}

if(nrow(fds2)>threshold2){
	n.pred2 <- round(((nrow(fds2)-threshold2)*vec[2])*FDS2)
}else{n.pred2 <- c(0,0,0,0)}

# remove the predated ones

pt1[is.na(pt1)] <- 0
n.pred1[is.na(n.pred1)] <- 0

pt2[is.na(pt2)] <- 0
n.pred2[is.na(n.pred2)] <- 0

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

# now get rid of the rest that will die, taking random individuals

LV1 <- round(threshold1-(nrow(fds1)-threshold1)*(vec[2]))
LV2 <- round(threshold2-(nrow(fds2)-threshold2)*(vec[2]))

if(nrow(fin1)>threshold1){
	fin1 <- fin1[1:LV1,]
} else if(nrow(fin1)>1){fin1 <- fin1}else{fin1 <- c()}


if(nrow(fin2)>threshold2){
	fin2 <- fin2[1:LV2,]
} else if(nrow(fin2)>1){fin2 <- fin2}else{fin2<- c()}


fin <- list(fin1, fin2)

pops[[z+1]] <- fin

}

# get the frequency differences between populations using this code
#diffs <- lapply(pops, freqDiffs)

#fMat <- matrix(unlist(diffs), ncol=6, byrow=T)

return(pops)

}

# this is all to make figures

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



x2 <- list()


for(q in 1:length(ltest)){
x2[[q]] <- migHel2pop(ltest[[q]])	
}

# we have a list of lists of lists - each element of the 36 main list is a list of n.gen matrices, each of which
# has two lists for the two populations. We want the change in allele frequency of each population over the generations

# make a function to pull allele frequencies from the 2-pop sub list for one generation


time2 <- afTime(x2[[59]])
time1 <- afTime(x[[59]])

par(mfrow=c(2,1))
par(mar=c(2,2,2,2))

plot(seq(0,1, by=1/ngen), type="n", main="Predator sees one")
lines(time2[[1]][,3], lty=3, col="red")
lines(time2[[1]][,1], col="red")
lines(time2[[2]][,3], lty=3,  col="blue")
lines(time2[[2]][,1], col="blue")


plot(seq(0,1, by=1/ngen), type="n", main="Predator sees both")
lines(time1[[1]][,3], lty=3, col="red")
lines(time1[[1]][,1], col="red")
lines(time1[[2]][,3],lty=3,  col="blue")
lines(time1[[2]][,1], col="blue")

# could we do something with mean/sd across generations?

# 3-d plot

# get frequency differences in each generation, return average for each allele at each parameter value


fDs2 <- freqDiffs(x2)

# make a matrix of percent migrate vs. percent nfds for each allele
# the rows will be percent migrate (6 rows), the columns percent nfds (also 6)

barMat2 <- matrix(fDs2[,1], nrow=length(pm1), byrow=F)

ulMat2 <- matrix(fDs2[,5], nrow=length(pm1), byrow=F)



# plots!

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))

persp(pm1, nfds1, barMat2,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.3))

persp(pm1, nfds1, ulMat2,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(0,0.3))

