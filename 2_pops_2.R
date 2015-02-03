# I made most of the chunks into indpendent functions so that we can run them on as many
# populations as we want

start.pop <- 500
LF <- 0.3
percent.migrate <- 0.5
percent.breed <- 0.5
carrying.capacity <- 2000
baseAttack <- c(.01, .1, .1, .3)
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

n.gen <- 10

geno1 <- matrix(rbinom(start.pop*6, 1, (1/3)), ncol=6)
colnames(geno1) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

geno2 <- matrix(rbinom(start.pop*6, 1, (1/3)), ncol=6)
colnames(geno2) <- c("bands1", "bands2", "red1", "red2", "neutral1", "neutral2")

recombination=rbinom(start.pop, 1, LF)

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

#make.off(4, geno1, start.pop, percent.breed)

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

next.gen.2=do.call(rbind, phenolist2)

return(next.gen.2)	
}

# normal LV selection

LV <- function(NFmat, carrying.capacity, percent.breed, n.off){
rate.inc <- percent.breed*n.off
	
threshold <- abs(nrow(NFmat)+(nrow(NFmat)*rate.inc*(1-(nrow(NFmat)/carrying.capacity))))
	
if(threshold > nrow(NFmat)){
	rand=sample(nrow(NFmat))
	next.gen.1=NFmat[rand,]
	next.gen=next.gen.1	
	}else{
	rand=sample(nrow(NFmat))
	next.gen.1=next.gen.2[rand,]
	next.gen=next.gen.1[1:threshold,]
}

return(next.gen)
}

#LV(NF1, carrying.capacity, percent.breed, n.off)
# parameters

# make two alleles worth of genotypes - don't differentiate sexes - these are the first elements in a list



# make the linked alleles

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

NF1 <- NFDS(pg1, ba, sim, hand)

NF2 <- NFDS(pg2, ba, sim, hand)


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


pops0 <- lapply(pops, allele.freq)

diff.red.0 <- c()
diff.linked.0 <- c()
diff.unlinked.0 <- c()

for(i in 1:n.gen){
diff.red.0[i] <- abs(pops0[[i]][[1]][2]-pops0[[i]][[2]][2])
diff.linked.0[i] <- abs(pops0[[i]][[1]][3]-pops0[[i]][[2]][3])
diff.unlinked.0[i] <- abs(pops0[[i]][[1]][4]-pops0[[i]][[2]][4])
}

pops0.1 <- lapply(pops, allele.freq)


diff.red.01 <- c()
diff.linked.01 <- c()
diff.unlinked.01 <- c()

for(i in 1:n.gen){
diff.red.01[i] <- abs(pops0.1[[i]][[1]][2]-pops0.1[[i]][[2]][2])
diff.linked.01[i] <- abs(pops0.1[[i]][[1]][3]-pops0.1[[i]][[2]][3])
diff.unlinked.01[i] <- abs(pops0.1[[i]][[1]][4]-pops0.1[[i]][[2]][4])
}

pops0.2 <- lapply(pops, allele.freq)

diff.red.02 <- c()
diff.linked.02 <- c()
diff.unlinked.02 <- c()

for(i in 1:n.gen){
diff.red.02[i] <- abs(pops0.2[[i]][[1]][2]-pops0.2[[i]][[2]][2])
diff.linked.02[i] <- abs(pops0.2[[i]][[1]][3]-pops0.2[[i]][[2]][3])
diff.unlinked.02[i] <- abs(pops0.2[[i]][[1]][4]-pops0.2[[i]][[2]][4])
}

pops0.3 <- lapply(pops, allele.freq)

diff.red.03 <- c()
diff.linked.03 <- c()
diff.unlinked.03 <- c()

for(i in 1:n.gen){
diff.red.03[i] <- abs(pops0.3[[i]][[1]][2]-pops0.3[[i]][[2]][2])
diff.linked.03[i] <- abs(pops0.3[[i]][[1]][3]-pops0.3[[i]][[2]][3])
diff.unlinked.03[i] <- abs(pops0.3[[i]][[1]][4]-pops0.3[[i]][[2]][4])
}


pops0.4 <- lapply(pops, allele.freq)


diff.red.04 <- c()
diff.linked.04 <- c()
diff.unlinked.04 <- c()

for(i in 1:n.gen){
diff.red.04[i] <- abs(pops0.4[[i]][[1]][2]-pops0.4[[i]][[2]][2])
diff.linked.04[i] <- abs(pops0.4[[i]][[1]][3]-pops0.4[[i]][[2]][3])
diff.unlinked.04[i] <- abs(pops0.4[[i]][[1]][4]-pops0.4[[i]][[2]][4])
}


pops0.5 <- lapply(pops, allele.freq)

diff.red.05 <- c()
diff.linked.05 <- c()
diff.unlinked.05 <- c()

for(i in 1:n.gen){
diff.red.05[i] <- abs(pops0.5[[i]][[1]][2]-pops0.5[[i]][[2]][2])
diff.linked.05[i] <- abs(pops0.5[[i]][[1]][3]-pops0.5[[i]][[2]][3])
diff.unlinked.05[i] <- abs(pops0.5[[i]][[1]][4]-pops0.5[[i]][[2]][4])
}

plot(seq(0,0.5,1/20), seq(0,1,1/10), type="n", xlab="migration", ylab="difference in allele freq")

points(x=rep(0, length=n.gen), y=diff.red.0, col="red")
points(x=rep(0.005, length=n.gen), y=diff.unlinked.0, col="blue")
points(x=rep(0.01, length=n.gen), y=diff.linked.0, col="yellow")

points(x=rep(0.1, length=n.gen), y=diff.red.01, col="red")
points(x=rep(0.105, length=n.gen), y=diff.unlinked.01, col="blue")
points(x=rep(0.11, length=n.gen), y=diff.linked.01, col="yellow")

points(x=rep(0.2, length=n.gen), y=diff.red.02, col="red")
points(x=rep(0.205, length=n.gen), y=diff.unlinked.02, col="blue")
points(x=rep(0.21, length=n.gen), y=diff.linked.02, col="yellow")

points(x=rep(0.3, length=n.gen), y=diff.red.03, col="red")
points(x=rep(0.305, length=n.gen), y=diff.unlinked.03, col="blue")
points(x=rep(0.31, length=n.gen), y=diff.linked.03, col="yellow")

points(x=rep(0.4, length=n.gen), y=diff.red.04, col="red")
points(x=rep(0.405, length=n.gen), y=diff.unlinked.04, col="blue")
points(x=rep(0.41, length=n.gen), y=diff.linked.04, col="yellow")


points(x=rep(0.5, length=n.gen), y=diff.red.05, col="red")
points(x=rep(0.505, length=n.gen), y=diff.unlinked.05, col="blue")
points(x=rep(0.51, length=n.gen), y=diff.linked.05, col="yellow")

###################################################################
#
###################################################################

allele.freq.br <- function(list){
	a1 <- colSums(list[[1]])/nrow(list[[1]])
	af1 <- c(mean(a1[2], a1[3]), mean(a1[4], a1[5]))
	a2 <- colSums(list[[2]])/nrow(list[[2]])
	af2 <- c(mean(a2[2], a2[3]), mean(a2[4], a2[5]))
	return(list(af1, af2))
}


br.freq <- lapply(pops, allele.freq.br)

plot(x=seq(0,1,by=0.1), y=seq(0,1,by=0.1), type="n")

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









