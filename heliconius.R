# make the initial population
rm(list=ls())

start.pop <- 50
p1start <- 1/2
p2start <- 1/2

rep.dist <- 0.8
rep.nondist <- 0.2

LF <- 0.3

matchSurv <- 0.9
nomatchSurv <- 0.1

n.mig <- 

make.off <- function(n.off, mat){

pair1 <- mat[1:(nrow(mat)/2),]
pair2 <- mat[(1+nrow(mat)/2):nrow(mat),]

bar.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
spot.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
sprec.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
mphrec.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
ul.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))
linked.off1 <- matrix(nrow=n.off, ncol=nrow(pair1))

for(i in 1:nrow(pair1)){
	which.bar <- rbinom(n.off, 1, 0.5)+1
	bar.off1[,i] <- pair1[i,2:3][which.bar]
	spot.off1[,i] <- pair1[i,4:5][which.bar]
	which.sprec <- rbinom(n.off, 1, 0.5)+1
	sprec.off1[,i] <- pair1[i,6:7][which.sprec]
	which.mphrec <- rbinom(n.off, 1, 0.5)+1
	mphrec.off1[,i] <- pair1[i,8:9][which.mphrec]
	which.ul <- rbinom(n.off, 1, 0.5)+1
	ul.off1[,i] <- pair1[i,10:11][which.ul]		
	linked.off1[,i] <- pair1[i,12:13][which.bar]
}

bar.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
spot.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
sprec.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
mphrec.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
ul.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))
linked.off2 <- matrix(nrow=n.off, ncol=nrow(pair2))

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

offspring <- cbind(as.vector(bar.off1), as.vector(bar.off2), as.vector(spot.off1),
as.vector(spot.off2), as.vector(sprec.off1), as.vector(sprec.off2), as.vector(mphrec.off1),
as.vector(mphrec.off2), as.vector(ul.off1), as.vector(ul.off2))

recombination <- rbinom(nrow(offspring), 1, LF)

linked1 <- matrix(NA, nrow=nrow(offspring), ncol=2)

for(i in 1:nrow(offspring)){
	if(recombination[i]==0){linked1[i,] <- offspring[,4:5][i,]} else
		linked1[i,]<- offspring[,4:5][i,c(2,1)]
}

recombSpot <- rbinom(nrow(offspring), 1, LF)

# do the recombination
spotPre <- cbind(as.vector(spot.off1), as.vector(spot.off2))

spot1 <- matrix(NA, nrow=nrow(spotPre), ncol=2)

for(i in 1:nrow(spotPre)){
	if(recombSpot[i]==0){spot1[i,] <- spotPre[,1:2][i,]} else
	spot1[i,]<- spotPre[,1:2][i,c(2,1)]
}


offspring <- cbind(offspring[,1:2], spot1, offspring[,5:10], linked1)
colnames(offspring) <- c("bar1", "bar2", "spot1", "spot2", "sprec1", "sprec2", "mphrec1", "mphrec2","ul1", "ul2", "linked1", "linked2")

return(offspring)
}

phenotype <- function(vec, mat){
genos <- mat[,vec]

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

# "bar" allele
# "spot" allele
# species recognition allele
# morph recognition allele
# linked allele
# unlinked allele

get.heliconius.offspring <- function(mat){
# phenotype
# four phenotypes from the two color alleles

pheno <- phenotype(1:4,mat)

# breeding
# split species distinguishing from not
pg <- cbind(pheno, mat)

# first divide by species recognition
# phenotype species recognition

sprec <- rowSums(cbind(mat[,5], mat[,6]))
sprec[sprec==2]=1

rec <- list()

if(sum(sprec==1)/length(sprec)==0){rec[[1]]=pg[1:round(nrow(rec)*rep.nondist),]
	} else if(sum(sprec==1)/length(sprec)==1){rec[[1]]=pg[1:round(nrow(nr)*rep.dist),]
	} else {
	r <- pg[sprec==1,]
	rec[[1]] <- r[1:round(nrow(r)*rep.dist),]
	nr <- pg[sprec==0,]
	rec[[2]] <- nr[1:round(nrow(nr)*rep.nondist),]
}

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

off <- list()	

for(i in 1:length(morphBred)){
	if(class(morphBred[[i]])=="numeric"){
		off[[i]] <- morphBred[[i]][2:13]
		} else if(nrow(morphBred[[i]]) < 6){
			off[[i]] <- morphBred[[i]][,2:13]
		} else {
			off[[i]] <- make.off(3, morphBred[[i]])
		}
}


offspring <- do.call(rbind, off)
offspring <- offspring[sample(nrow(offspring)),]

return(offspring)

}
####################################

####################################
# predator sees single pop #########
####################################

vec=c(0.1, 0.1)

migHel1pop <- function(vec){
	
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

geno2 <- matrix(rbinom(start.pop*8, 1, p1start), ncol=8)

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

pops[[1]] <- list(geno1, geno2)

i=1

for(i in 1:ngen){
		# make offspring
	g2 <- pops[[i]][[2]][,2:13]	
	g1 <- pops[[i]][[1]][,2:13]

	# migrate
	n.mig <- round(nrow(g1)*vec[1])

	if(n.mig==0){
	geno1m <- g1
	geno2m <- g2
	}else{
	geno1m <- rbind(g2[1:n.mig,], g1[(n.mig+1):nrow(g1),])

	geno2m <- rbind(g1[1:n.mig,], g2[(n.mig+1):nrow(g2),])
	}
	
	# make offspring
	off1 <- get.heliconius.offspring(geno1m)
	off2 <- get.heliconius.offspring(geno2m)
	
	# make phenotypes
	
	pheno1 <- phenotype(1:4, off1)
	pheno2 <- phenotype(1:4, off2)
	
	pg1 <- cbind(pheno1, off1)
	order1 <- order(pg1[,1])
	pg1 <- pg1[order1,]
	
	pg2 <- cbind(pheno2, off2)
	order2 <- order(pg2[,1])
	pg2 <- pg2[order2,]
}

}

# selection
phenoSel <- phenotype(1:4,offspring)

offSel <- cbind(phenoSel, offspring)

oneBS <- offSel[offSel[,1]==1,]
twoBS <- offSel[offSel[,1]==2,]
threeBS <- offSel[offSel[,1]==3,]
fourBS <- offSel[offSel[,1]==4,]

oneS <- oneBS[1:round(nrow(oneBS)*nomatchSurv),]
twoS <- twoBS[1:round(nrow(twoBS)*matchSurv),]
threeS <- threeBS[1:round(nrow(threeBS)*nomatchSurv),]
fourS <- fourBS[1:round(nrow(fourBS)*nomatchSurv),]


selNames <- c("oneS", "twoS", "threeS", "fourS")

selList <- list()
selSub <- c()

for(z in 1:4){
	if(exists(selNames[z])){selList[[z]] <- get(selNames[z])}
	else(selList[[z]] <- selSub)
}

survivors <- do.call(rbind,selList)

survivors <- survivors[sample(nrow(survivors)),]



# migration




















