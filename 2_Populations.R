###########
#functions#
###########

# this is a test comment to see what github does

#This function takes in a vector of 8 possible genotypes, and organizes them
#by simple dominance into 4 phenotypes (morphs)

#ADD 4 genotypes to represent locus 3 neutral allele, each gives one of four existing morphs to maintain original ratio

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

#Second population
# we should need only one of these - it's a function, so we can add any input we want
phenotype2=function(offspring.phenotype2){
offspring.phenotype1.2=ifelse(offspring.phenotype2==0, 1, offspring.phenotype2)

offspring.phenotype2.2=ifelse(offspring.phenotype2==1, 2, offspring.phenotype1.2)

offspring.phenotype3.2=ifelse(offspring.phenotype2==2, 2, offspring.phenotype2.2)

offspring.phenotype4.2=ifelse(offspring.phenotype2==3, 3, offspring.phenotype3.2)

offspring.phenotype5.2=ifelse(offspring.phenotype2==4, 4, offspring.phenotype4.2)

offspring.phenotype6.2=ifelse(offspring.phenotype2==5, 4, offspring.phenotype5.2)

offspring.phenotype7.2=ifelse(offspring.phenotype2==6, 3, offspring.phenotype6.2)

offspring.phenotype8.2=ifelse(offspring.phenotype2==7, 4, offspring.phenotype7.2)

offspring.phenotype9.2=ifelse(offspring.phenotype2==8, 4, offspring.phenotype8.2)

offspring.phenotype10.2=ifelse(offspring.phenotype2==9, 1, offspring.phenotype9.2)

offspring.phenotype11.2=ifelse(offspring.phenotype2==10, 2, offspring.phenotype10.2)

offspring.phenotype12.2=ifelse(offspring.phenotype2==11, 3, offspring.phenotype11.2)

offspring.phenotype13.2=ifelse(offspring.phenotype2==12, 4, offspring.phenotype12.2)

return(offspring.phenotype13.2)
}


#The function produces "n.off" offspring genotypes from 2 parental 
#genotypes
#n.off = 4
get.offspring=function(v){
which.allele=rbinom(n.off, 1, 0.5)+1
offspring.set=c()
for(i in 1:n.off){
offspring.set[i]=v[which.allele[i]]
}
return(offspring.set)
}

##########
#FUNCTION#
##########

# We put the starting genotype percentages outside the function. We
# calculate the percentage of each of 2 alleles at each of 2 loci, 
# separated into male and female groups

#Repeat for population 2

geno.percent=c()

geno.percent.2=c()

# morph percentages are calculated from the genotype percentages
# the 5th entry is the population size at each generation

morph.percent=c()

morph.percent.2=c()

#starting population size
start.pop=40

#fill the genotype matrix: male locus 1, chr. 1, male locus 1, chr. 2, 
# female locus 1, chr.1, female locus 1, chr.2, and the same for locus 2 #AND NEUTRAL LOCUS 3: male locus 3, chr. 1, male locus 3, chr. 2; female locus 3, chr. 1, female locus 3, chr. 2

one=rbinom(start.pop/2, 1, (1/3))
two=rbinom(start.pop/2, 1, (1/3))
three=rbinom(start.pop/2, 1, (1/3))
four=rbinom(start.pop/2, 1, (1/3))
five=rbinom(start.pop/2, 1, (1/3))
six=rbinom(start.pop/2, 1, (1/3))
seven=rbinom(start.pop/2, 1, (1/3))
eight=rbinom(start.pop/2, 1, (1/3))
nine=rbinom(start.pop/2, 1, (1/3))
ten=rbinom(start.pop/2, 1, (1/3))
eleven=rbinom(start.pop/2, 1, (1/3))
twelve=rbinom(start.pop/2, 1, (1/3))

one.2=rbinom(start.pop/2, 1, (1/3))
two.2=rbinom(start.pop/2, 1, (1/3))
three.2=rbinom(start.pop/2, 1, (1/3))
four.2=rbinom(start.pop/2, 1, (1/3))
five.2=rbinom(start.pop/2, 1, (1/3))
six.2=rbinom(start.pop/2, 1, (1/3))
seven.2=rbinom(start.pop/2, 1, (1/3))
eight.2=rbinom(start.pop/2, 1, (1/3))
nine.2=rbinom(start.pop/2, 1, (1/3))
ten.2=rbinom(start.pop/2, 1, (1/3))
eleven.2=rbinom(start.pop/2, 1, (1/3))
twelve.2=rbinom(start.pop/2, 1, (1/3))


# multiply the locus 2 by 3, so that each potential 2-locus diploid genotype
# has a unique number associated with it ----> And multiply locus 3 by 5
# males
geno0=rowSums(cbind(one, two, five*3, six*3, nine*5, ten*5))

geno0.2=rowSums(cbind(one.2, two.2, five.2*3, six.2*3, nine.2*5, ten.2*5))
#females
geno1=rowSums(cbind(three, four, seven*3, eight*3, eleven*5, twelve*5))

geno1.2=rowSums(cbind(three.2, four.2, seven.2*3, eight.2*3, eleven.2*5, twelve.2*5))

geno3=append(geno0, geno1)

geno3.2=append(geno0.2, geno1.2)
# change genotype codes to phenotypes
geno2=phenotype(geno3)

geno2.2=phenotype(geno3.2)

# input the percentages of each genotype into the geno.percent vector
geno.percent[1]=sum(one)/(start.pop/2)
geno.percent[2]=sum(two)/(start.pop/2)
geno.percent[3]=sum(three)/(start.pop/2)
geno.percent[4]=sum(four)/(start.pop/2)
geno.percent[5]=sum(five)/(start.pop/2)
geno.percent[6]=sum(six)/(start.pop/2)
geno.percent[7]=sum(seven)/(start.pop/2)
geno.percent[8]=sum(eight)/(start.pop/2)

geno.percent.2[1]=sum(one.2)/(start.pop/2)
geno.percent.2[2]=sum(two.2)/(start.pop/2)
geno.percent.2[3]=sum(three.2)/(start.pop/2)
geno.percent.2[4]=sum(four.2)/(start.pop/2)
geno.percent.2[5]=sum(five.2)/(start.pop/2)
geno.percent.2[6]=sum(six.2)/(start.pop/2)
geno.percent.2[7]=sum(seven.2)/(start.pop/2)
geno.percent.2[8]=sum(eight.2)/(start.pop/2)



# make the morph.percent vector
morph.percent[1]=length(geno2[geno2==1])/start.pop
morph.percent[2]=length(geno2[geno2==2])/start.pop
morph.percent[3]=length(geno2[geno2==3])/start.pop
morph.percent[4]=length(geno2[geno2==4])/start.pop
morph.percent[5]=start.pop

morph.percent.2[1]=length(geno2.2[geno2.2==1])/start.pop
morph.percent.2[2]=length(geno2.2[geno2.2==2])/start.pop
morph.percent.2[3]=length(geno2.2[geno2.2==3])/start.pop
morph.percent.2[4]=length(geno2.2[geno2.2==4])/start.pop
morph.percent.2[5]=start.pop

#########EXCHANGE MIGRANTS PRIOR TO BREEDING##########
geno.matrix=matrix(geno.percent)
geno.matrix2=matrix(geno.percent.2)



migrants<-geno.matrix[1:3,]
migrants2<-geno.matrix2[1:3,]
mig=matrix(migrants)
mig2=matrix(migrants2)

#bind migrants to new population

geno.matrix.bind=rbind(geno.matrix, mig2)
geno.matrix.2bind=rbind(geno.matrix2, mig)


#remove migrants from initial populations

geno.percent.new<-geno.matrix.new[-(1:3),]
geno.percent.2new<-geno.matrix.2new[-(1:3),]

#generate new phenotype matrix:

# multiply the locus 2 by 3, so that each potential 2-locus diploid genotype
# has a unique number associated with it ----> And multiply locus 3 by 5
# males
new.geno0=rowSums(cbind(one, two, five*3, six*3, nine*5, ten*5))

geno0.2=rowSums(cbind(one.2, two.2, five.2*3, six.2*3, nine.2*5, ten.2*5))
#females
geno1=rowSums(cbind(three, four, seven*3, eight*3, eleven*5, twelve*5))

geno1.2=rowSums(cbind(three.2, four.2, seven.2*3, eight.2*3, eleven.2*5, twelve.2*5))

geno3=append(geno0, geno1)

geno3.2=append(geno0.2, geno1.2)
# change genotype codes to phenotypes
geno2=phenotype(geno3)

geno2.2=phenotype(geno3.2)






# function to 1) put one generation through a round of breeding, 2) morph-specific selection, 
# 3) random selection around a population threshold

# inputs are:

# general parameters for the population in the absence of predation

# Carrying capacity (one number) used to set step 3

# n.off sets the number of offspring per parental pair - not yet variable, but could be

# percent.breed = percentage of incoming breeding population that does breed. currently
# the same for males and females, but could be made variable.

# parameters for the markov chain

# base.attack = a four element vector that described the tendency of the predators
# to attack a given morph, in the absence of any learning

# similarity is a 4x4 matrix that describes the similarity between morphs as perceived
# by the predator. Similarity is bounded between 0 and 1. Diagonal elements are the 
# similarity of a morph to itself, and so should usually be 1. This matrix needs to 
# be symmetrical.

# handling is a 4x4 matrix describing the inverse of the handling time for each morph given 
# that the predator previously handled a morph 1 - 4

one.gen=function(carrying.capacity, n.off, percent.breed, base.attack, similarity, handling, 
geno.percent, morph.percent, geno.percent.2, morph.percent.2 LF){

n.off = 4
percent.breed = 0.5
ngen = 100
carrying.capacity = 200

base.attack = c(0.1, 0.1, 0.1, 0.1)

s1=c(1,.1,.1,.1)
s2=c(.1,1,.1,.1)
s3=c(.1,.1,1,.1)
s4=c(.1,.1,.1,1)

similarity=rbind(s1,s2,s3,s4)

T1=c(1,1,1,1)
T2=c(1,1,1,1)
T3=c(1,1,1,1)
T4=c(1,1,1,1)

handling=rbind(T1,T2,T3,T4)



#make output vectors
morph.percent.out=c()
geno.percent.out=c()

# size describes the number of male and female individuals - currently 
# a 50:50 sex ratio

	size=morph.percent[5]/2
	
# generate diploid genotypes for males and females separately at each 
# locus

#ADD LOCUS 3
	
	males1=cbind(rbinom(size, 1, geno.percent.new[1]),rbinom(size, 1, geno.percent.new[2]))
	LF=0.3
	recombination=rbinom(nrow(males1), 1, LF)
	malesLinked=matrix(NA,ncol=2,nrow=nrow(males1))
	
	#for second population
	
	males1.2=cbind(rbinom(size, 1, geno.percent.2new[1]),rbinom(size, 1, geno.percent.2new[2]))
	recombination.2=rbinom(nrow(males1.2), 1, LF)
	malesLinked.2=matrix(NA,ncol=2,nrow=nrow(males1.2))

	
	for(i in 1:nrow(males1)){
	if(recombination[i]==0){malesLinked[i,]<-males1[i,]} else
		malesLinked[i,]<- males1[i,c(2,1)]
	}
	
	for(i in 1:nrow(males1.2)){
	if(recombination.2[i]==0){malesLinked.2[i,]<-males1.2[i,]} else
		malesLinked.2[i,]<- males1.2[i,c(2,1)]
	}
	
	males2=cbind(rbinom(size, 1, geno.percent.new[5]),rbinom(size, 1, geno.percent.new[6]))
	males2.2=cbind(rbinom(size, 1, geno.percent.2new[5]),rbinom(size, 1, geno.percent.2new[6]))
	
	females1=cbind(rbinom(size, 1, geno.percent.new[3]),rbinom(size, 1, geno.percent.new[4]))
	females1.2=cbind(rbinom(size, 1, geno.percent.2new[3]),rbinom(size, 1, geno.percent.2new[4]))
	
		recombination2=rbinom(nrow(females1), 1, LF)
		recombination2.2=rbinom(nrow(females1.2), 1, LF)
	
	femalesLinked=matrix(NA,ncol=2,nrow=nrow(females1))
		femalesLinked.2=matrix(NA,ncol=2,nrow=nrow(females1.2))
	
	for(i in 1:nrow(females1)){
	if(recombination2[i]==0){femalesLinked[i,]<-females1[i,]} else
		femalesLinked[i,]<- females1[i,c(2,1)]
	}
	
	for(i in 1:nrow(females1.2)){
	if(recombination2.2[i]==0){femalesLinked.2[i,]<-females1.2[i,]} else
		femalesLinked.2[i,]<- females1.2[i,c(2,1)]
	}
	
	females2=cbind(rbinom(size, 1, geno.percent.new[7]),rbinom(size, 1, geno.percent.new[8]))
	females2.2=cbind(rbinom(size, 1, geno.percent.2new[7]),rbinom(size, 1, geno.percent.2new[8]))
	
# identify the individuals that get to breed. Currently no reproductive skew
	
	lucky.males1=sample(nrow(males1), (percent.breed)*(nrow(males1)))
	lucky.males1.2=sample(nrow(males1.2), (percent.breed)*(nrow(males1.2)))

	lucky.males2=sample(nrow(males2), (percent.breed)*(nrow(males1)))
	lucky.males2.2=sample(nrow(males2.2), (percent.breed)*(nrow(males1.2)))

	lucky.females1=sample(nrow(females1), (percent.breed)*(nrow(males1)))
	lucky.females1.2=sample(nrow(females1.2), (percent.breed)*(nrow(males1.2)))
	
	lucky.females2=sample(nrow(females2), (percent.breed)*(nrow(males1)))
	lucky.females2.2=sample(nrow(females2.2), (percent.breed)*(nrow(males1.2)))

	
# reduce the matrices to the lucky ones - breeding is completely random in
# this population

####EXPAND male and female matrices (locus 1) to reflect linkage with neutral locus

	male.pairs1=males1[lucky.males1,]
	male.pairs1.2=males1.2[lucky.males1.2,]
	maleLinked.pairs <- malesLinked[lucky.males1,]
	maleLinked.pairs.2 <- malesLinked.2[lucky.males1.2,]

	male.pairs2=males2[lucky.males2,]
	male.pairs2.2=males2.2[lucky.males2.2,]

	female.pairs1=females1[lucky.females1,]
	female.pairs1.2=females1.2[lucky.females1.2,]
	
	femaleLinked.pairs <- femalesLinked[lucky.females1,]
	femaleLinked.pairs.2 <- femalesLinked.2[lucky.females1.2,]

	female.pairs2=females2[lucky.females2,]
	female.pairs2.2=females2.2[lucky.females2.2,]

# generate male gametes (4 per male) for the two loci
#Alter for linkage between locus one and neutral locas 3 (retain get.offspring for unaffected locus)	

#(Locus 1 + linked):

loc1 <- matrix(nrow=4, ncol=10)
linked <- matrix(nrow=4, ncol=10)

loc1.2 <- matrix(nrow=4, ncol=10)
linked.2 <- matrix(nrow=4, ncol=10)


for(i in 1:10){
	which.allele <- rbinom(4, 1, 0.5)+1
	loc1[,i] <- male.pairs1[i,][which.allele]
	linked[,i] <- maleLinked.pairs[i,][which.allele]
}

for(i in 1:10){
	which.allele <- rbinom(4, 1, 0.5)+1
	loc1.2[,i] <- male.pairs1.2[i,][which.allele]
	linked.2[,i] <- maleLinked.pairs.2[i,][which.allele]
}


#(Locus 2):
	offspring1.2=matrix(NA, ncol=nrow(male.pairs1), nrow=n.off)
	offspring1.2.2=matrix(NA, ncol=nrow(male.pairs1.2), nrow=n.off)


for(m in 1:nrow(male.pairs1)){
	which.allele <- rbinom(4, 1, 0.5)+1
	offspring1.2[,m]=male.pairs2[m,][which.allele]	
}

for(m in 1:nrow(male.pairs1.2)){
	which.allele <- rbinom(4, 1, 0.5)+1
	offspring1.2.2[,m]=male.pairs2.2[m,][which.allele]	
}


male.offspring <- cbind(as.vector(loc1), as.vector(linked), as.vector(offspring1.2))
male.offspring.2 <- cbind(as.vector(loc1.2), as.vector(linked.2), as.vector(offspring1.2.2))

	


# generate female gametes

loc1.female <- matrix(nrow=4, ncol=10)
linked.female <- matrix(nrow=4, ncol=10)


loc1.female.2 <- matrix(nrow=4, ncol=10)
linked.female.2 <- matrix(nrow=4, ncol=10)


for(i in 1:10){
	which.allele <- rbinom(4, 1, 0.5)+1
	loc1.female[,i] <- female.pairs1[i,][which.allele]
	linked.female[,i] <- femaleLinked.pairs[i,][which.allele]
}

for(i in 1:10){
	which.allele <- rbinom(4, 1, 0.5)+1
	loc1.female.2[,i] <- female.pairs1.2[i,][which.allele]
	linked.female.2[,i] <- femaleLinked.pairs.2[i,][which.allele]
}


#Locus 2 female:	
	offspring2.2=matrix(NA, ncol=nrow(female.pairs1), nrow=n.off)
	offspring2.2.2=matrix(NA, ncol=nrow(female.pairs1.2), nrow=n.off)



for(o in 1:nrow(female.pairs1)){
	offspring2.2[,o]=get.offspring(female.pairs2[o,])	
}

for(o in 1:nrow(female.pairs1.2)){
	offspring2.2.2[,o]=get.offspring(female.pairs2.2[o,])	
}
	offspring.2.2=as.vector(offspring2.2)
	offspring.2.2.2=as.vector(offspring2.2.2)
	
	female.offspring <- cbind(as.vector(loc1.female), as.vector(linked.female), offspring.2.2)
	female.offspring.2 <- cbind(as.vector(loc1.female.2), as.vector(linked.female.2), offspring.2.2.2)

	
# make offspring diploid genotypes from the gametes
#alter to include linked gametes


	offspring=cbind(male.offspring[,1], female.offspring[,1], male.offspring[,2], female.offspring[,2], male.offspring[,3], female.offspring[,3])
	
	offspring.2=cbind(male.offspring.2[,1], female.offspring.2[,1], male.offspring.2[,2], female.offspring.2[,2], male.offspring.2[,3], female.offspring.2[,3])

	
	colnames(offspring) <- c("loc1.1", "loc1.2", "linked.1", "linked.2", "loc2.1", "loc2.2")
    colnames(offspring.2) <- c("loc1.1", "loc1.2", "linked.1", "linked.2", "loc2.1", "loc2.2")

	
	offspring.phenotype=rowSums(cbind(offspring[,1:2], offspring[,5]*3, offspring[,6]*3))
	offspring.phenotype.2=rowSums(cbind(offspring.2[,1:2], offspring.2[,5]*3, offspring.2[,6]*3))


# get the offspring phenotypes

	off.pheno=phenotype(offspring.phenotype)
	off.pheno.2=phenotype(offspring.phenotype.2)


# get the parent phenotypes - we assume that the non-parents have already died - 
# maybe we should revisit that? It would only take replacing the  male.pairs
# and female.pairs vectors with males1, males2, females1 and females2

	parent.phenotype=rowSums(cbind(rbind(males1, females1), rbind(3*males2, 3*females2)))
	parent.phenotype.2=rowSums(cbind(rbind(males1.2, females1.2), rbind(3*males2.2, 3*females2.2)))
	
	parent.pheno=phenotype(parent.phenotype)
	parent.pheno.2=phenotype(parent.phenotype.2)

# get the phenotypes and genotypes for the entire population

	phenotypes=append(off.pheno, parent.pheno)
	phenotypes.2=append(off.pheno.2, parent.pheno.2)
	

	genotypes=rbind(offspring, cbind(males1, malesLinked, males2), cbind(females1, femalesLinked, females2))
	genotypes.2=rbind(offspring.2, cbind(males1.2, malesLinked.2, males2.2), cbind(females1.2, femalesLinked.2, females2.2))

	pheno.geno=cbind(phenotypes, genotypes)
	pheno.geno.2=cbind(phenotypes.2, genotypes.2)

	pg1=order(pheno.geno[,1])
	pg=pheno.geno[pg1,]
	
	pg1.2=order(pheno.geno.2[,1])
	pg.2=pheno.geno.2[pg1.2,]


# split into four matrices, one for each phenotype
	
	pt=c(sum(pg[,1]==1), sum(pg[,1]==2), sum(pg[,1]==3), sum(pg[,1]==4))
	
	pheno_1=pg[pg[,1]==1,]
	pheno_2=pg[pg[,1]==2,]
	pheno_3=pg[pg[,1]==3,]
	pheno_4=pg[pg[,1]==4,]

pt.2=c(sum(pg.2[,1]==1), sum(pg.2[,1]==2), sum(pg.2[,1]==3), sum(pg.2[,1]==4))
	
	pheno_1.2=pg.2[pg.2[,1]==1,]
	pheno_2.2=pg.2[pg.2[,1]==2,]
	pheno_3.2=pg.2[pg.2[,1]==3,]
	pheno_4.2=pg.2[pg.2[,1]==4,]




##################
# this is the negative frequency dependence
#treat populations 1 and 2 as composite population for purposes of predation
N=c(pt,pt.2)
# first we find the total number of predators, in terms of 
# base attack rate, similarity, handling time, and number of prey

denom1=matrix(NA, nrow=4, ncol=4)

for(k in 1:4){
	for (j in 1:4){
		denom1[j,k]=base.attack[k]*N[k]*(1+similarity[k,j]*handling[k,j]*base.attack[j]*N[j])
	}
}

denom=sum(denom1)

# then we find the relative number of each morph that is eaten by the predators

f1=matrix(NA, nrow=4, ncol=4)
for(l in 1:4){
	for(m in 1:4){
	f1[l,m]=base.attack[l]*N[l]*similarity[l,m]*base.attack[m]*N[m]
	}
	
}
f=colSums(f1)

# these induce negative frequency dependence because a predator is more
# likely to eat a morph when it's last meal was also that morph. Common
# morphs are more frequently encountered, and so are more likely to be 
# that last meal - therefore they experience more mortality than rare morphs.s

# now we find the absolute number of each morph eaten by the predators,
# and clip the matrices by this number. The rest of this code is so that
# the matrices don't give us an error message when we index a matrix
# with no rows, or index from 0:1 or 1:1

surv1=round(N-N*(f/denom))
surv=(abs(surv1)+surv1)/2


both=ifelse(pt<surv, pt, surv)
phenolist=list(pheno_1, pheno_2,pheno_3,pheno_4)
phenolist2=list()
phenosub=c()

for(q in 1:4){
	if(both[q]>1){phenolist2[[q]]=phenolist[[q]][1:both[q],]}
	else if(both[q]==1){phenolist2[[q]]=phenolist[[q]]}
	else if(both[q]==0){phenolist2[[q]]=phenosub}
	else{phenolist2[[q]]=phenosub}
}


both2=ifelse(pt.2<surv, pt.2, surv)
pop2.phenolist=list(pheno_1.2, pheno_2.2, pheno_3.2, pheno_4.2)
pop2.phenolist2=list()
phenosub2=c()


for(r in 1:4){
	if(both2[r]>1){pop2.phenolist2[[r]]=pop2.phenolist[[2]][1:both2[r],]}
	else if(both2[r]==1){pop2.phenolist2[[r]]=pop2.phenolist[[r]]}
	else if(both2[r]==0){pop2.phenolist2[[r]]=phenosub2}
	else{pop2.phenolist2[[r]]=phenosub2}
}





# now we have a matrix of individuals that survived the morph-specific
# predation

next.gen.2=do.call(rbind, phenolist2)
pop2.next.gen.2=do.call(rbind, pop2.phenolist2)


##########################

# now we induce random, lotka-volterra type mortality on a randomized list

	last.gen=morph.percent[5]
	rate.inc=(percent.breed)*n.off
	threshold=abs(last.gen+(last.gen*rate.inc*(1-(last.gen/carrying.capacity))))
	if(threshold > nrow(next.gen.2)){
	rand=sample(nrow(next.gen.2))
	next.gen.1=next.gen.2[rand,]
	next.gen=next.gen.1	
	}else{
	rand=sample(nrow(next.gen.2))
	next.gen.1=next.gen.2[rand,]
	next.gen=next.gen.1[1:threshold,]
	}

#and fill in the output vectors

	#morph.percent.out[5]=length(next.gen[,1])
	
	#morph.percent.out[1]=(table(next.gen[,1])["1"])/length(next.gen[,1])
	#morph.percent.out[2]=(table(next.gen[,1])["2"])/length(next.gen[,1])
	#morph.percent.out[3]=(table(next.gen[,1])["3"])/length(next.gen[,1])
	#morph.percent.out[4]=(table(next.gen[,1])["4"])/length(next.gen[,1])
	
	#geno.percent.out[1]=nrow(next.gen[next.gen[,2]==1,])/nrow(next.gen)
	#geno.percent.out[2]=nrow(next.gen[next.gen[,3]==1,])/nrow(next.gen)
	#geno.percent.out[3]=nrow(next.gen[next.gen[,4]==1,])/nrow(next.gen)
	#geno.percent.out[4]=nrow(next.gen[next.gen[,5]==1,])/nrow(next.gen)
	#geno.percent.out[5]=nrow(next.gen[next.gen[,2]==1,])/nrow(next.gen)
	#geno.percent.out[6]=nrow(next.gen[next.gen[,3]==1,])/nrow(next.gen)
	#geno.percent.out[7]=nrow(next.gen[next.gen[,4]==1,])/nrow(next.gen)
	#geno.percent.out[8]=nrow(next.gen[next.gen[,5]==1,])/nrow(next.gen)
	#geno.percent.out[9]=nrow(next.gen[next.gen[,2]==1,])/nrow(next.gen)
	#geno.percent.out[10]=nrow(next.gen[next.gen[,3]==1,])/nrow(next.gen)
	#geno.percent.out[11]=nrow(next.gen[next.gen[,4]==1,])/nrow(next.gen)
	#geno.percent.out[12]=nrow(next.gen[next.gen[,5]==1,])/nrow(next.gen)

return(list(next.gen))
}


C=c(.01, .1, .1, .3)

s1=c(1,.1,.1,.1)
s2=c(.1,1,.1,.1)
s3=c(.1,.1,1,.1)
s4=c(.1,.1,.1,1)

s=rbind(s1,s2,s3,s4)

T1=c(1,1,1,1)
T2=c(1,1,1,1)
T3=c(1,1,1,1)
T4=c(1,1,1,1)

t=rbind(T1,T2,T3,T4)

ngen=10000

geno.percent=matrix(0, nrow=ngen, ncol=12)
morph.percent=matrix(0, nrow=ngen, ncol=5)
geno.percent.2=matrix(0, nrow=ngen, ncol=12)
morph.percent.2=matrix(0, nrow=ngen, ncol=5)

one=rbinom(start.pop/2, 1, (1/3))
two=rbinom(start.pop/2, 1, (1/3))
three=rbinom(start.pop/2, 1, (1/3))
four=rbinom(start.pop/2, 1, (1/3))
five=rbinom(start.pop/2, 1, (1/3))
six=rbinom(start.pop/2, 1, (1/3))
seven=rbinom(start.pop/2, 1, (1/3))
eight=rbinom(start.pop/2, 1, (1/3))
nine=rbinom(start.pop/2, 1, (1/3))
ten=rbinom(start.pop/2, 1, (1/3))
eleven=rbinom(start.pop/2, 1, (1/3))
twelve=rbinom(start.pop/2, 1, (1/3))

one.2=rbinom(start.pop/2, 1, (1/3))
two.2=rbinom(start.pop/2, 1, (1/3))
three.2=rbinom(start.pop/2, 1, (1/3))
four.2=rbinom(start.pop/2, 1, (1/3))
five.2=rbinom(start.pop/2, 1, (1/3))
six.2=rbinom(start.pop/2, 1, (1/3))
seven.2=rbinom(start.pop/2, 1, (1/3))
eight.2=rbinom(start.pop/2, 1, (1/3))
nine.2=rbinom(start.pop/2, 1, (1/3))
ten.2=rbinom(start.pop/2, 1, (1/3))
eleven.2=rbinom(start.pop/2, 1, (1/3))
twelve.2=rbinom(start.pop/2, 1, (1/3))


geno0=rowSums(cbind(one, two,five*3, six*3, nine*5, ten*5))
geno1=rowSums(cbind(three, four, seven*3, eight*3, eleven*5, twelve*5))
geno3=append(geno0, geno1)
geno2=phenotype(geno3)

geno0.2=rowSums(cbind(one.2, two.2, five.2*3, six.2*3, nine.2*5, ten.2*5))
geno1.2=rowSums(cbind(three.2, four.2, seven.2*3, eight.2*3, eleven.2*5, twelve.2*5))
geno3.2=append(geno0.2, geno1.2)
geno2.2=phenotype(geno3.2)


geno.percent[1]=sum(one)/(start.pop/2)
geno.percent[2]=sum(two)/(start.pop/2)
geno.percent[3]=sum(three)/(start.pop/2)
geno.percent[4]=sum(four)/(start.pop/2)
geno.percent[5]=sum(five)/(start.pop/2)
geno.percent[6]=sum(six)/(start.pop/2)
geno.percent[7]=sum(seven)/(start.pop/2)
geno.percent[8]=sum(eight)/(start.pop/2)
geno.percent[9]=sum(nine)/(start.pop/2)
geno.percent[10]=sum(ten)/(start.pop/2)
geno.percent[11]=sum(eleven)/(start.pop/2)
geno.percent[12]=sum(twelve)/(start.pop/2)


geno.percent.2[1]=sum(one.2)/(start.pop/2)
geno.percent.2[2]=sum(two.2)/(start.pop/2)
geno.percent.2[3]=sum(three.2)/(start.pop/2)
geno.percent.2[4]=sum(four.2)/(start.pop/2)
geno.percent.2[5]=sum(five.2)/(start.pop/2)
geno.percent.2[6]=sum(six.2)/(start.pop/2)
geno.percent.2[7]=sum(seven.2)/(start.pop/2)
geno.percent.2[8]=sum(eight.2)/(start.pop/2)
geno.percent.2[9]=sum(nine.2)/(start.pop/2)
geno.percent.2[10]=sum(ten.2)/(start.pop/2)
geno.percent.2[11]=sum(eleven.2)/(start.pop/2)
geno.percent.2[12]=sum(twelve.2)/(start.pop/2)


morph.percent[1]=length(geno2[geno2==1])/start.pop
morph.percent[2]=length(geno2[geno2==2])/start.pop
morph.percent[3]=length(geno2[geno2==3])/start.pop
morph.percent[4]=length(geno2[geno2==4])/start.pop
morph.percent[5]=start.pop

morph.percent.2[1]=length(geno2.2[geno2.2==1])/start.pop
morph.percent.2[2]=length(geno2.2[geno2.2==2])/start.pop
morph.percent.2[3]=length(geno2.2[geno2.2==3])/start.pop
morph.percent.2[4]=length(geno2.2[geno2.2==4])/start.pop
morph.percent.2[5]=start.pop


#one.gen=function(carrying.capacity, n.off, start.pop, percent.breed, base.attack, similarity, handling, 
#geno.percent, morph.percent)


one.gen(carrying.capacity, n.off, percent.breed, C, s, t, geno.percent, morph.percent, geno.percent.2, morph.percent.2)


##########
#FOR LOOP#
##########

# now we can run the function over multiple generations


start.pop=40
n.gen=10

# this is my base attack vector
C=c(.01, .1, .1, .3)

# similarity matrix
s1=c(1,0,0,0)
s2=c(0,1,0,0)
s3=c(0,0,1,0)
s4=c(0,0,0,1)

s=rbind(s1,s2,s3,s4)

# handling time matrix
T1=c(1,1,1,1)
T2=c(1,1,1,1)
T3=c(1,1,1,1)
T4=c(1,1,1,1)

t=rbind(T1,T2,T3,T4)

# initiate the population

one=rbinom(start.pop/2, 1, (1/3))
two=rbinom(start.pop/2, 1, (1/3))
three=rbinom(start.pop/2, 1, (1/3))
four=rbinom(start.pop/2, 1, (1/3))
five=rbinom(start.pop/2, 1, (1/3))
six=rbinom(start.pop/2, 1, (1/3))
seven=rbinom(start.pop/2, 1, (1/3))
eight=rbinom(start.pop/2, 1, (1/3))
nine=rbinom(start.pop/2, 1, (1/3))
ten=rbinom(start.pop/2, 1, (1/3))
eleven=rbinom(start.pop/2, 1, (1/3))
twelve=rbinom(start.pop/2, 1, (1/3))

one.2=rbinom(start.pop/2, 1, (1/3))
two.2=rbinom(start.pop/2, 1, (1/3))
three.2=rbinom(start.pop/2, 1, (1/3))
four.2=rbinom(start.pop/2, 1, (1/3))
five.2=rbinom(start.pop/2, 1, (1/3))
six.2=rbinom(start.pop/2, 1, (1/3))
seven.2=rbinom(start.pop/2, 1, (1/3))
eight.2=rbinom(start.pop/2, 1, (1/3))
nine.2=rbinom(start.pop/2, 1, (1/3))
ten.2=rbinom(start.pop/2, 1, (1/3))
eleven.2=rbinom(start.pop/2, 1, (1/3))
twelve.2=rbinom(start.pop/2, 1, (1/3))


#EDIT HERE FOR LOCUS 3

geno0=rowSums(cbind(one, two,five*3, six*3, nine*5, ten*5))
geno1=rowSums(cbind(three, four, seven*3, eight*3, eleven*5, twelve*5))
geno3=append(geno0, geno1)
geno2=phenotype(geno3)

geno0.2=rowSums(cbind(one.2, two.2, five.2*3, six.2*3, nine.2*5, ten.2*5))
geno1.2=rowSums(cbind(three.2, four.2, seven.2*3, eight.2*3, eleven.2*5, twelve.2*5))
geno3.2=append(geno0.2, geno1.2)
geno2.2=phenotype(geno3.2)


GP=matrix(0, nrow=n.gen, ncol=12)
GP2=matrix(0, nrow=n.gen, ncol=12)
MP=matrix(0, nrow=n.gen, ncol=5)
MP2=matrix(0, nrow=n.gen, ncol=5)

GP[1,1]=sum(one)/(start.pop/2)
GP[1,2]=sum(two)/(start.pop/2)
GP[1,3]=sum(three)/(start.pop/2)
GP[1,4]=sum(four)/(start.pop/2)
GP[1,5]=sum(five)/(start.pop/2)
GP[1,6]=sum(six)/(start.pop/2)
GP[1,7]=sum(seven)/(start.pop/2)
GP[1,8]=sum(eight)/(start.pop/2)
GP[1,9]=sum(nine)/(start.pop/2)
GP[1,10]=sum(ten)/(start.pop/2)
GP[1,11]=sum(eleven)/(start.pop/2)
GP[1,12]=sum(twelve)/(start.pop/2)

GP2[1,1]=sum(one.2)/(start.pop/2)
GP2[1,2]=sum(two.2)/(start.pop/2)
GP2[1,3]=sum(three.2)/(start.pop/2)
GP2[1,4]=sum(four.2)/(start.pop/2)
GP2[1,5]=sum(five.2)/(start.pop/2)
GP2[1,6]=sum(six.2)/(start.pop/2)
GP2[1,7]=sum(seven.2)/(start.pop/2)
GP2[1,8]=sum(eight.2)/(start.pop/2)
GP2[1,9]=sum(nine.2)/(start.pop/2)
GP2[1,10]=sum(ten.2)/(start.pop/2)
GP2[1,11]=sum(eleven.2)/(start.pop/2)
GP2[1,12]=sum(twelve.2)/(start.pop/2)


MP[1,1]=length(geno2[geno2==1])/start.pop
MP[1,2]=length(geno2[geno2==2])/start.pop
MP[1,3]=length(geno2[geno2==3])/start.pop
MP[1,4]=length(geno2[geno2==4])/start.pop
MP[1,5]=start.pop

MP2[1,1]=length(geno2.2[geno2.2==1])/start.pop
MP2[1,2]=length(geno2.2[geno2.2==2])/start.pop
MP2[1,3]=length(geno2.2[geno2.2==3])/start.pop
MP2[1,4]=length(geno2.2[geno2.2==4])/start.pop
MP2[1,5]=start.pop




# For loop to add generations
# one.gen(carrying.capacity, n.off, percent.breed, C, s, t, geno.percent, morph.percent)

i=3
for(i in 2:10){
	GP[i,]=one.gen(200,4, 0.5, C, s, t, GP[i-1,], MP[i-1,])[[2]]
	MP[i,]=one.gen(200,4, 0.5, C, s, t, GP[i-1,], MP[i-1,])[[1]]
}

j=3
for(j in 2:10){
	GP2[j,]=one.gen(200,4, 0.5, C, s, t, GP2[i-1,], MP2[i-1,])[[2]]
	MP2[j,]=one.gen(200,4, 0.5, C, s, t, GP2[i-1,], MP2[i-1,])[[1]]
}

# plotting

plot(x=c(1:10000), y=seq(from=0, to=1, by=1/9999), type="n")
lines(MP[,1])
lines(MP[,2], col="red")
lines(MP[,3], col="blue")
lines(MP[,4], col="green")
lines(MP2[,1], col="darkgray")
lines(MP2[,2], col="orange")
lines(MP2[,3], col="purple")
lines(MP2[,4], col="yellow")



loc1=(GP[,1]+GP[,2]+GP[,5]+GP[,6])/4
loc2=(GP[,3]+GP[,4]+GP[,7]+GP[,8])/4
loc3=(GP[,9]+GP[,10]+GP[,11]+GP[,12])/4
loc1.2=(GP[,1]+GP[,2]+GP[,5]+GP[,6])/4
loc2.2=(GP[,3]+GP[,4]+GP[,7]+GP[,8])/4
loc3.2=(GP2[,9]+GP2[,10]+GP2[,11]+GP2[,12])/4
lines(loc1, col="darkgreen")
lines(loc2, col="darkblue")
lines(loc3, col="deeppink")


# concerns: at the moment, individual genotypes are not passed from one generation to the next
# instead new genotypes produced from the previous year's genotype ratios. This is probably
# fine in a large population, but would it bias a small population? The benefit is that doing
# it this way saves a lot of memory. We also currently have no age structure - mortality from
# predation and lotka volterra effects apply equally to young of the year and the parents.
# Thoughts?



