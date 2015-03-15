# NFDS
start.pop = 50
LF = 0.3
percent.breed = 0.5
n.gen = 20
thresh = 200

make.off <- function(n.off, mat, percent.breed){
	
lucky <- sample(nrow(mat), percent.breed*nrow(mat))

pairs <- mat[lucky,]

pair1 <- pairs[1:(nrow(pairs)/2),]
pair2 <- pairs[(1+nrow(pairs)/2):nrow(pairs),]

coding.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
linked.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
neutral.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))

for(i in 1:nrow(pair1)){
	which.coding <- rbinom(n.off, 1, 0.5)+1
	coding.off1[,i] <- pair1[i,1:2][which.coding]
	linked.off1[,i] <- pair1[i,5:6][which.coding]
	which.neu <- rbinom(n.off, 1, 0.5)+1
	neutral.off1[,i] <- pair1[i,3:4][which.neu]
}

coding.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
linked.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
neutral.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))

for(i in 1:nrow(pair2)){
	which.coding <- rbinom(n.off, 1, 0.5)+1
	coding.off2[,i] <- pair1[i,1:2][which.coding]
	linked.off2[,i] <- pair1[i,5:6][which.coding]
	which.neu <- rbinom(n.off, 1, 0.5)+1
	neutral.off2[,i] <- pair1[i,3:4][which.neu]
}


offspring <- cbind(as.vector(coding.off1), as.vector(coding.off2),as.vector(neutral.off1), as.vector(neutral.off2), as.vector(linked.off1), as.vector(linked.off2))

rec <- rbinom(nrow(offspring), 1, LF)

linked1 <- matrix(NA, nrow=nrow(offspring), ncol=2)

for(i in 1:nrow(offspring)){
	if(rec[i]==0){linked1[i,] <- offspring[,5:6][i,]} else
		linked1[i,]<- offspring[,5:6][i,c(2,1)]
}

offspring <- cbind(offspring[,1:4], linked1)

return(offspring)

}

geno1 <- matrix(rbinom(start.pop*4, 1, (1/3)), ncol=4)
colnames(geno1) <- c("coding1", "coding2", "neutral1", "neutral2")

geno2 <- matrix(rbinom(start.pop*4, 1, (1/3)), ncol=4)
colnames(geno2) <- c("coding1", "coding2", "neutral1", "neutral2")

recombination=rbinom(start.pop, 1, LF)

# do the recombination

linked1 <- matrix(NA, nrow=start.pop, ncol=2)

for(i in 1:start.pop){
	if(recombination[i]==0){linked1[i,] <- geno1[,1:2][i,]} else
		linked1[i,]<- geno1[,1:2][i,c(2,1)]
}
	
linked2 <- matrix(NA, nrow=start.pop, ncol=2)

for(i in 1:start.pop){
	if(recombination[i]==0){linked2[i,] <- geno2[,1:2][i,]} else
		linked2[i,] <- geno2[,1:2][i,c(2,1)]
}

geno1 <- cbind(geno1, linked1)
colnames(geno1) <- c("coding1", "coding2", "neutral1", "neutral2", "linked1", "linked2")
geno2 <- cbind(geno2, linked2)
colnames(geno2) <- c("coding1", "coding2", "neutral1", "neutral2", "linked1", "linked2")


#########################
# start the function ####
#########################


negative <- function(vec){

geno1 <- matrix(rbinom(start.pop*4, 1, (1/3)), ncol=4)
colnames(geno1) <- c("coding1", "coding2", "neutral1", "neutral2")

geno2 <- matrix(rbinom(start.pop*4, 1, (1/3)), ncol=4)
colnames(geno2) <- c("coding1", "coding2", "neutral1", "neutral2")

recombination1=rbinom(start.pop, 1, vec[2])

recombination2=rbinom(start.pop, 1, vec[2])

# do the recombination

linked1 <- matrix(NA, nrow=start.pop, ncol=2)

for(i in 1:start.pop){
	if(recombination1[i]==0){linked1[i,] <- geno1[,1:2][i,]} else
		linked1[i,]<- geno1[,1:2][i,c(2,1)]
}
	
linked2 <- matrix(NA, nrow=start.pop, ncol=2)

for(i in 1:start.pop){
	if(recombination2[i]==0){linked2[i,] <- geno2[,1:2][i,]} else
		linked2[i,] <- geno2[,1:2][i,c(2,1)]
}

pheno1 <- rowSums(geno1[,1:2])
pheno1[pheno1==2]=1

pheno2 <- rowSums(geno2[,1:2])
pheno2[pheno2==2]=1

geno1 <- cbind(pheno1, geno1, linked1)
colnames(geno1) <- c("pheno","coding1", "coding2", "neutral1", "neutral2", "linked1", "linked2")
geno2 <- cbind(pheno2, geno2, linked2)
colnames(geno2) <- c("pheno","coding1", "coding2", "neutral1", "neutral2", "linked1", "linked2")

pops <- list()

pops[[1]] <- list(geno1, geno2)

#i=1
for(i in 1:20){

g1 <- pops[[i]][[1]][,2:7]
g2 <- pops[[i]][[2]][,2:7]

n.mig <- round(nrow(g1)*vec[1])

if(n.mig==0){
	geno1m <- g1
	geno2m <- g2
}else{
geno1m <- rbind(g2[1:n.mig,], g1[(n.mig+1):nrow(g1),])

geno2m <- rbind(g1[1:n.mig,], g2[(n.mig+1):nrow(g2),])
}

off1 <- make.off(10, geno1m, 0.8)
off2 <- make.off(10, geno2m, 0.8)
# negative frequency dependent selection

pheno1 <- rowSums(off1[,1:2])
pheno1[pheno1==2]=1

percent1 <- sum(pheno1)/length(pheno1)

pg1 <- cbind(pheno1, off1)
ones1 <- pg1[pg1[,1]==1,]
zeros1 <- pg1[pg1[,1]==0,]

ones1 <- ones1[1:round(nrow(ones1)*(1-percent1)),]
zeros1 <- zeros1[1:round(nrow(zeros1)*(percent1)),]

ng1 <- rbind(ones1, zeros1)
ng1 <- ng1[sample(nrow(ng1)),]

if(nrow(ng1)<thresh){ng1 <- ng1}else(ng1 <- ng1[1:thresh,])

######################

pheno2 <- rowSums(off2[,1:2])
pheno2[pheno2==2]=1

percent2 <- sum(pheno2)/length(pheno2)

pg2 <- cbind(pheno2, off2)
ones2 <- pg2[pg2[,1]==1,]
zeros2 <- pg2[pg2[,1]==0,]

ones2 <- ones2[1:round(nrow(ones2)*(1-percent2)),]
zeros2 <- zeros2[1:round(nrow(ones2)*(percent2)),]

ng2 <- rbind(ones2, zeros2)
ng2 <- ng2[sample(nrow(ng2)),]
if(nrow(ng2)<thresh){ng2 <- ng2}else(ng2 <- ng2[1:thresh,])

pops[[i+1]] <- list(ng1, ng2)
}

return(pops)
}

###########################################
pm1 <- seq(0, 0.1, by=0.05)
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


R <- list()

for(m in 1:10){

# ltest is the same for each iteration, but re-running migLD will get us different
# starting points and progression through the generations
	
af <- lapply(ltest, negative)

diff <- list()

for(k in 1:18){
diff[[k]] = matrix(nrow=21, ncol=7)
for(i in 1:21){
diff[[k]][i,] <- abs(colMeans(af[[k]][[i]][[1]]) - colMeans(af[[k]][[i]][[2]]))
}
}


means <- lapply(diff, function(p){x <- colMeans(p); return(x)})

R[[m]]	<- matrix(unlist(means), ncol=7, byrow=T)
}

# get the mean of the means across runs - each row is an allele 
# bands, red, linked, unlinked
# each row is a set of parameter values

mean <- Reduce('+', R, R[[1]])/10

mean2 <- cbind(rowMeans(cbind(mean[,2], mean[,3])), rowMeans(cbind(mean[,4], mean[,5])),
rowMeans(cbind(mean[,6], mean[,7])))

# take the mean values for the "band" locus, make them into a matrix
# with values of pm along the rows and values of rf for the columns

xcodeMeans <- matrix(mean2[,1], ncol=length(rf1))

xulMeans <- matrix(mean2[,2], ncol=length(rf1))

xlinkedMeans <- matrix(mean2[,3], ncol=length(rf1))

NFDScodeCorrected <- xcodeMeans-xulMeans
NFDSlinkedCorrected <- xlinkedMeans-xulMeans

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))

persp(pm1, rf1, NFDScodeCorrected,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(-.15,0.15))

persp(pm1, rf1, NFDSlinkedCorrected, theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(-.15,0.3))


##########################################
# positive frequency dependence ##########
##########################################

#vec = c(0,0)

positive <- function(vec){

geno1 <- matrix(rbinom(start.pop*4, 1, (1/3)), ncol=4)
colnames(geno1) <- c("coding1", "coding2", "neutral1", "neutral2")

geno2 <- matrix(rbinom(start.pop*4, 1, (1/3)), ncol=4)
colnames(geno2) <- c("coding1", "coding2", "neutral1", "neutral2")

recombination1=rbinom(start.pop, 1, vec[2])

recombination2=rbinom(start.pop, 1, vec[2])

# do the recombination

linked1 <- matrix(NA, nrow=start.pop, ncol=2)

for(j in 1:start.pop){
	if(recombination1[j]==0){linked1[j,] <- geno1[,1:2][j,]} else
		linked1[j,]<- geno1[,1:2][j,c(2,1)]
}
	
linked2 <- matrix(NA, nrow=start.pop, ncol=2)

for(n in 1:start.pop){
	if(recombination2[n]==0){linked2[n,] <- geno2[,1:2][n,]} else
		linked2[n,] <- geno2[,1:2][n,c(2,1)]
}

pheno1 <- rowSums(geno1[,1:2])
pheno1[pheno1==2]=1

pheno2 <- rowSums(geno2[,1:2])
pheno2[pheno2==2]=1

geno1 <- cbind(pheno1, geno1, linked1)
colnames(geno1) <- c("pheno","coding1", "coding2", "neutral1", "neutral2", "linked1", "linked2")
geno2 <- cbind(pheno2, geno2, linked2)
colnames(geno2) <- c("pheno","coding1", "coding2", "neutral1", "neutral2", "linked1", "linked2")

pops <- list()

pops[[1]] <- list(geno1, geno2)

#i=1
for(i in 1:20){

g1 <- pops[[i]][[1]][,2:7]
g2 <- pops[[i]][[2]][,2:7]

n.mig <- round(nrow(g1)*vec[1])

if(n.mig==0){
	geno1m <- g1
	geno2m <- g2
}else{
geno1m <- rbind(g2[1:n.mig,], g1[(n.mig+1):nrow(g1),])

geno2m <- rbind(g1[1:n.mig,], g2[(n.mig+1):nrow(g2),])
}

off1 <- make.off(10, geno1m, 0.8)
off2 <- make.off(10, geno2m, 0.8)
# negative frequency dependent selection

pheno1 <- rowSums(off1[,1:2])
pheno1[pheno1==2]=1

percent1 <- sum(pheno1)/length(pheno1)

pg1 <- cbind(pheno1, off1)
ones1 <- pg1[pg1[,1]==1,]
zeros1 <- pg1[pg1[,1]==0,]

genlength <- c(nrow(ones1), nrow(zeros1))
per <- c(percent1, (1-percent1))

genolist <- list(ones1, zeros1)
genolist2 <- list()
sub <- c()

for(q in 1:2){
	if(genlength[q]>1){genolist2[[q]] <- genolist[[q]][1:(nrow(genolist[[q]])*per[q]),]}
	else if(genlength[q]==1){genolist2[[q]] <- genolist[[q]]}
	else if(genlength[q]==0){genolist2[[q]] <- sub}
	else{genolist2[[q]] <- sub}
}

ng1 <- do.call(rbind, genolist2)

ng1 <- ng1[sample(nrow(ng1)),]

if(nrow(ng1)<thresh){ng1 <- ng1}else(ng1 <- ng1[1:thresh,])

######################

pheno2 <- rowSums(off2[,1:2])
pheno2[pheno2==2]=1

percent2 <- sum(pheno2)/length(pheno2)

pg2 <- cbind(pheno2, off2)
ones2 <- pg2[pg2[,1]==1,]
zeros2 <- pg2[pg2[,1]==0,]

genlength2 <- c(nrow(ones2), nrow(zeros2))
per2 <- c(percent2, (1-percent2))

genolist1.2 <- list(ones2, zeros2)
genolist2.2 <- list()
sub2 <- c()

for(b in 1:2){
	if(genlength2[b]>1){genolist2.2[[b]] <- genolist1.2[[b]][1:(nrow(genolist1.2[[b]])*per2[b]),]}
	else if(genlength2[b]==1){genolist2.2[[b]] <- genolist1.2[[b]]}
	else if(genlength2[b]==0){genolist2.2[[b]] <- sub2}
	else{genolist2.2[[b]] <- sub2}
}

ng2 <- do.call(rbind, genolist2.2)

ng2 <- ng2[sample(nrow(ng2)),]

if(nrow(ng2)<thresh){ng2 <- ng2}else(ng2 <- ng2[1:thresh,])

pops[[i+1]] <- list(ng1, ng2)
}

return(pops)
}

POS <- list()

for(m in 1:10){

# ltest is the same for each iteration, but re-running migLD will get us different
# starting points and progression through the generations
	
af2 <- lapply(ltest, positive)

diffPos <- list()

for(w in 1:18){
diffPos[[w]] = matrix(nrow=21, ncol=7)
for(z in 1:21){
diffPos[[w]][z,] <- abs(colMeans(af2[[w]][[z]][[1]]) - colMeans(af2[[w]][[z]][[2]]))
}
}

meansPos <- lapply(diffPos, function(p){x <- colMeans(p); return(x)})

POS[[5]] <- matrix(unlist(meansPos), ncol=7, byrow=T)
}

# get the mean of the means across runs - each row is an allele 
# bands, red, linked, unlinked
# each row is a set of parameter values

meanPos <- Reduce('+', POS, POS[[1]])/5

meanPos2 <- cbind(rowMeans(cbind(meanPos[,2], meanPos[,3])), rowMeans(cbind(meanPos[,4], meanPos[,5])),
rowMeans(cbind(meanPos[,6], meanPos[,7])))

# take the mean values for the "band" locus, make them into a matrix
# with values of pm along the rows and values of rf for the columns

xcodeMeansPos <- matrix(meanPos2[,1], ncol=length(rf1))

xulMeansPos <- matrix(meanPos2[,2], ncol=length(rf1))

xlinkedMeansPos <- matrix(meanPos2[,3], ncol=length(rf1))


PFDScodeCorrected <- xcodeMeansPos-xulMeansPos
PFDSlinkedCorrected <- xlinkedMeansPos-xulMeansPos

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))

persp(pm1, rf1, NFDScodeCorrected,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(-.5,0.5))

persp(pm1, rf1, NFDSlinkedCorrected, theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(-.5,0.5))

persp(pm1, rf1, PFDScodeCorrected,theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(-.5,0.5))

persp(pm1, rf1, PFDSlinkedCorrected, theta=30, phi=30, col="lightblue", shade=0.4,
ticktype="detailed", zlim=c(-.5,0.5))