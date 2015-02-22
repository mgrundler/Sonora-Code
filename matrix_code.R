
start.pop <- 20
ncols <- 2
nrows <- 2
npops <- ncols*nrows
npred <- 0.5
percent.breed <- 0.8
n.off <- 2
ngen <- 10

# generate populations

pophist <- list()

for(k in 1:ngen){
pops <- list()

for(i in 1:npops){
	p <- matrix(rbinom(start.pop*2, 1, (1/3)), ncol=2)
	phen <- rowSums(p)
	phen[phen>1]=1
	pops[[i]]<- cbind(phen, p)
}

# migration
# pick which populations will migrate to the other
# should we restrict to one individual per generation?
# select a source and sink population

source <- c(sample(ncols,1), sample(nrows,1))
sink <- c(sample(ncols,1), sample(nrows,1))

# take one from source to sink

mat <- matrix(1:npops, ncol=ncols, nrow=nrows, byrow=T)

pops[[mat[sink[1],sink[2]]]] <- rbind(pops[[mat[source[1],source[2]]]][1,], pops[[mat[sink[1],sink[2]]]])

# frequency dependent selection

fds <- function(mat){
	mat <- mat[order(mat[,1]),]
	p1 <- sum(mat[,1])/nrow(mat)
	amtpred0 <- round(npred*nrow(mat)*(1-p1))
	amtpred1 <- round(npred*nrow(mat))-amtpred0
	mat0 <- mat[mat[,1]==0,]
	mat0out <- mat0[1:(nrow(mat0)-amtpred0),]
	mat1 <- mat[mat[,1]==1,]
	mat1out <- mat1[1:(nrow(mat1)-amtpred1),]
	matOut <- rbind(mat0out, mat1out)
	return(matOut[sample(nrow(matOut)),])
}

popsFds <- lapply(pops, fds)

# reproduction

make.off <- function(mat){
	
lucky <- sample(nrow(mat), percent.breed*nrow(mat))

pairs <- mat[lucky,]

pair1 <- pairs[1:(nrow(pairs)/2),]
pair2 <- pairs[(1+nrow(pairs)/2):nrow(pairs),]

off1 <- matrix(nrow=n.off, ncol=nrow(pair1))

for(i in 1:nrow(pair1)){
	which <- rbinom(n.off, 1, 0.5)+1
	off1[,i] <- pair1[i,1:2][which]
}

off2 <- matrix(nrow=n.off, ncol=nrow(pair1))

for(i in 1:nrow(pair1)){
	which <- rbinom(n.off, 1, 0.5)+1
	off2[,i] <- pair2[i,1:2][which]
}

offspring <- cbind(as.vector(off1), as.vector(off2))
phenoff <- rowSums(offspring)
phenoff[phenoff>1]=1
offspring1 <- cbind(phenoff, offspring)

newpop <- rbind(mat, offspring1)
return(newpop)
}

pophist[[k]] <- lapply(popsFds, make.off)
}






