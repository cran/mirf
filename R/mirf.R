#*****************************************************************************************
#Author:Yimin Wu
#Function for Multiple Imputation and Random Forests(MIRF)
#***************************************************************************************** 
#Input 
#geno:a data frame or matrix of alleles, such that each locus has a pair of adjacent 
#         columns of  alleles, and the order of columns corresponds to the order of loci 
#         on a chromosome. If  there are n loci, then ncol(geno) = 2*n. Rows represent 
#         the alleles for each subject. 

#y: a vector of responses for the subjects, length(y) = nrow(geno). Here we accept both 
#   categorical and continuous y. If y is categorical, "residualize" is not allowed.

#gene.column: a vector represents how many columns in the geno matrix belonging to 
#             each gene  separately. If there are k genes and n loci, then 
#             length(gene.column) = k, and sum (gene.column) = ncol(geno) = 2*n.
#             If user doesn't input this variable, then we assume all columns in geno 
#		  matrix belong to one gene. In another word, we set the default value to
#             c(2*n)

#gene.names: vector of names for genes

#SNPnames: vector of names for SNPs.

#M: number of haplotype-imputations within each simulation

#hwe.group : a vector of racial/ethnic group information about the subjects. The default
#            value is set to  NULL. If given racial/ethnic group information, haplotype 
#            frequency estimates are abtained  within each racial/ethnic group separately 
#            due to potential violations of Hardy Weinberg  equilibrium(HWE). 

#covariates: a matrix of covariates, e.g. sex,age. Defaultly  set to NULL.

#residualize: a vector indicates which columns in covariates matrix are used in residualization. 
#             Defaultly set to NULL. If ncol(covariates) = 3, residualize = c(1,3), then we 
#             use the first  and the third columns in covariates matrix to residualize, and 
#             directly include left  columns into randomForest algorithm. If residualize = NULL, 
#             then we include all columns in  covariates matrix into randomForest algorithm.

###############################
#parameters used in haplo.em
################################
#miss.val:vector of values that represent missing alleles in geno. It's a optional parameter in haplo.em
#weight:weights for observations (rows of geno matrix).It's a optional parameter in haplo.em
#control:list of control parameters. 
#        The default is constructed by the function haplo.em.control. 
#        The default behavior of this function results in 
#        the following parameter settings: 
#        loci.insert.order=1:n.loci, insert.batch.size=min(4,n.loci), min.posterior= 0.0001, 
#        tol=0.00001, max.iter=500, random.start=0 (no random start), 
#        iseed=NULL (no saved seed to start random start), verbose=0 (no printout during EM iterations). 
#        See haplo.em.control for more details.

###############################
#parameters used in randomForest
################################
#...:    optional parameters to be passed to the low level function randomForest.

mirf <- function (geno, y, gene.column=NULL, gene.names=NULL, SNPnames=NULL, M=10, hwe.group=NULL, covariates=NULL, residualize=NULL, 
		 miss.val=c(0, NA, "NA"), weight=NULL, control=haplo.em.control(),
             ...)
{
	if(M<=1)
	{
		return("M should be larger than 1! Choose another M and try again!")
	}
	conflictPckg <- "genetics"
	loadedPckg <- (.packages())
	for (i in 1:length(loadedPckg))
	{
		if (conflictPckg == loadedPckg[i])
		{
			return("There is a conflict between the genetics package and the haplo.em package, remove the genetics package and try again!")
		}
	} 

	if(!is.null(covariates))
	{
		covariates <- as.matrix(covariates)
	}
	
	#If residualize is not NULL, then we fit a linear model for the residualized covariates
	#Note: if y is categorical, then we don't allow residualize
	if(!(is.null(residualize)))
	{
		#test whether y is continuous
		if(!is.numeric(y))
		{
			return("When y is categorical, residualize is not allowed!")
		}
		options(contrasts = c("contr.treatment", "contr.treatment")) 
		lmformula <- paste("covariates[,",residualize,sep="")
		lmformula <- paste(lmformula,"]",sep="")
		lmformula <- paste("y~",paste(lmformula,collapse="+"),sep="")
		model1<-lm(as.formula(lmformula))
		if(length(model1$resid)!=length(y))
		{ 
			return("There might have some NA values in the residualized covariate, remove these rows and try again!")
		}
		y<-model1$resid
		
	}	
	N <- dim(geno)[[1]];  #number of individual samples
	if(is.null(gene.column))
	{
		gene.column = c(dim(geno)[[2]])
	}
	G <- length(gene.column); #number of genes
	cum.gene.column <- cumsum(gene.column);
	loci <- list();  # a list of loci names for all genes
	for (i in 1:G)
	{
		if(is.null(SNPnames)){
		#loci[i] is a list of loci names for the ith gene 
			loci[i] <- list(paste("gene_", paste(i, 1:(gene.column[i]/2), sep=""), sep=""));
		}else{
			if(i==1)
			{
				loci[i] <- list(SNPnames[1:(cum.gene.column[1]/2)])
			}else
			{
				loci[i] <- list(SNPnames[(cum.gene.column[i-1]/2+1):(cum.gene.column[i]/2)])
			}
		}
	}

	library(haplo.stats)
	print("haplo.em is processing......")
	groupnum <- 1; # number of hwe groups, default is 1
	if(is.null(hwe.group)) #if hwe.group is NULL, that means we assume all individuals are in one group.
	{
		hwe.group <- rep(1,times=N)
	}
	groupnum <- length(unique(hwe.group)) # total number of groups
	ingroup.num <- array(0,c(1,groupnum)) # store number of samples in each group
	haploem.groups <- list();	
	for (i in 1:groupnum)
	{
		ingroup.num[i] <- length(hwe.group[][hwe.group==unique(hwe.group)[i]])
		haploem <- list() # a list of all haplo.em results for all genes of group i
		for (j in 1:G)
		{
			if(j == 1)
				haploem[j] <- list(haplo.em(geno=geno[,1:gene.column[1]][hwe.group==unique(hwe.group)[i],],locus.label=loci[[j]], miss.val, weight, control))		
			else
				haploem[j] <- list(haplo.em(geno=geno[,(cum.gene.column[j-1]+1):cum.gene.column[j]][hwe.group==unique(hwe.group)[i],],locus.label=loci[[j]],miss.val, weight, control))	
		}
		haploem.groups[[i]] <- haploem	
	}
	cum.ingroup.num <- cumsum(c(0,ingroup.num))

	haploname.freq.gene <- c()	
	haplonames <- c()
	for(i in 1:G)
	{
		haploname.freq <- c()
		for (j in 1:groupnum)
		{
			#tempGeneName <- strsplit(dimnames(haploem.groups[[j]][[i]]$haplotype)[[2]][1],split="_")[[1]][1]	
			#tempHaploFre <- c(tempHaploFre,haploem.groups[[j]][[i]]$hap.prob)
			if(is.null(gene.names)){
				tempGeneName <- paste("gene",i,sep="");
			}else{
				tempGeneName <- gene.names[i];
			}
			for(k in 1:dim(haploem.groups[[j]][[i]]$haplotype)[1])
			{
		 		tempHaploType <- paste(haploem.groups[[j]][[i]]$haplotype[k,],collapse='')
				existHaploType <- 0
				if(length(haploname.freq)>0){
				for(tempIndex in 1:length(haploname.freq))	
				{
					if(haploname.freq[tempIndex][1]==tempHaploType)
					{
						haploname.freq[tempIndex,1+j]=haploem.groups[[j]][[i]]$hap.prob[k]
						existHaploType <- 1
						break
					}
				}
				}
				if(existHaploType == 0)
				{
					haploname.freq <- rbind(haploname.freq,c(tempHaploType,c(rep("NA",j-1),haploem.groups[[j]][[i]]$hap.prob[k],rep("NA",groupnum-j)),tempGeneName))
				}
			}
		}
		haploname.freq.gene <- rbind(haploname.freq.gene,haploname.freq)
		haplonames <- c(haplonames,haploname.freq[,1])
	}	
	haplonum <- length(haplonames) # store the total number of haplotypes
	names <- haplonames
	



	if(!(is.null(covariates)))
	{
		haplonum <- haplonum + dim(covariates)[[2]] - length(residualize)
		#if user didn't input the name of covariates, we use C1,C2,C3.... as the name of covariates
		if(is.null(dimnames(covariates)[[2]]))
		{
			dimnames(covariates)[[2]] <- as.list(paste("C",1:dim(covariates)[[2]],sep=""))
		}
		if(!(is.null(residualize)))
			names <- c(haplonames,dimnames(covariates)[[2]][-residualize])	
		else
			names <- c(haplonames,dimnames(covariates)[[2]])
	}

	#variables for storing importance statistics
	theta.bar.D<-matrix(data=NA,nrow=1,ncol=haplonum)
	W.bar.D<-matrix(data=NA,nrow=1,ncol=haplonum)
	B.D<-matrix(data=NA,nrow=1,ncol=haplonum)
	T.D<-matrix(data=NA,nrow=1,ncol=haplonum)
	teststats<-matrix(data=NA,nrow=1,ncol=haplonum)
	dimnames(teststats)[[2]] <- names	
	df<-matrix(data=NA,nrow=1,ncol=haplonum)
	#start simulations
	#***************************************************************************************************
	#  Iterative sampling of haplotype pairs and fitting a random forest, iterating M times 
	#  each of theses produces standardized and unstandardized score as well as corresponding standrad error
	#  from RF.       
	#  Note that haplo.em objects comprise the 1)indx.subj 2)nreps 3)post(posterior prob) 
	#***************************************************************************************************		
	
	#matrix to keep standardized varimportance of the  haplotypes at each iteration
	stdimportance<-matrix(data=0,nrow=M,ncol=haplonum)
	
	#matrix to keep raw varimportance of the haplotypes at each iteration
	varimportance<-matrix(data=0,nrow=M,ncol=haplonum)
 	dimnames(varimportance)[[2]]<- names
 	
	#matrix of within-imputation standard errors   
	Wd<-matrix(data=0,nrow=M,ncol=haplonum)

	##begin resampling and fitting a random forest each time		
	for (k in 1:M)
	{
		print(paste("iteration=",k,sep=""))
		hap.indic <- array(0,c(N, length(haplonames))) #binary matrix for indicating whether a haplotype exists in samples
		dimnames(hap.indic)[[2]]<-haplonames
		#iterate different groups and different genes
		for (i in 1:groupnum)
		{
			for(j in 1:G)
			{
				attach(haploem.groups[[i]][[j]])
				hap.pairs<-matrix(data=0,ingroup.num[i],ncol=2)
				
				for (s in 1:length(nreps))
				{
					if(nreps[s]>1)
   					{
    						hap.pairs[s,]<-cbind(hap1code,hap2code)[indx.subj==s,][sample(seq(1:nreps[s]),size=1,prob=post[indx.subj==s]),]	
	     				}
	     				#if there is only one haplotype pair for a subject take that one only
	      			if(nreps[s]==1)
    	  				{
 						hap.pairs[s,]<-cbind(hap1code,hap2code)[indx.subj==s,]
          				}	
     				} 
				for (s in (cum.ingroup.num[i]+1):cum.ingroup.num[i+1])
				{
					hap.indic[s,][dimnames(hap.indic)[[2]]==paste(haplotype[hap.pairs[s-cum.ingroup.num[i],1],],collapse="")|dimnames(hap.indic)[[2]]==paste(haplotype[hap.pairs[s-cum.ingroup.num[i],2],],collapse="") ] <- 1
				}		     		    		     		
				detach(haploem.groups[[i]][[j]])
			}
		}
		haplodata<-as.data.frame(hap.indic)  
		#include covariates
		if(!is.null(covariates))
		{
			if(!is.null(residualize))
			{
				haplodata <- cbind(haplodata,covariates[,-residualize])
			}
			else
			{
				haplodata <- cbind(haplodata,covariates)
			}	
			dimnames(haplodata)[[2]]<-names		
		}
		#####Now "haplodata" will be entered into the random forest	
		library(randomForest)  #fit a random forest to haplodata
		tempdata <- cbind(haplodata,y)
		forest <- randomForest(y~., data=tempdata, na.action=na.omit,importance=TRUE,proximity=TRUE,...) #in order to handle missing data
		
		
		stdimportance[k,]<-importance(forest)[,1] #standardized importance scores 
		varimportance[k,]<-forest$importance[,1] #raw importance scores in model with 1 dominant haplotype
		Wd[k,]<-forest$importanceSD #standard deviations of raw scores
	}#end of iterative resampling / MI
		
	### storing the estimates for inference for the current (s'th) simulation
	### this is according to the end of algorithm 2 of the MIRF paper
	for (j in 1:haplonum)
	{
    		theta.bar.D[1,j]<-(1/M)*sum(varimportance[,j]) #average imp. score over m imputations
      	W.bar.D[1,j]<-(1/M)*sum(Wd[,j]^2) #Within imputation variance
      	B.D[1,j]<-(1/(M-1))*sum((varimportance[,j]-theta.bar.D[1,j])^2) #between imputation variance
      	T.D[1,j]<-W.bar.D[1,j]+((M+1)/M)*B.D[1,j]# total variance V(jM)
      	teststats[1,j]<-theta.bar.D[1,j]/sqrt(T.D[1,j]) #Final importance statistic over all M imputations
      	df[1,j]<-(M-1)*(1+(1/(M+1))*W.bar.D[1,j]/B.D[1,j])^2
  	} 

	output <- t(rbind(dimnames(teststats)[[2]],teststats))
      dimnames(output)[[1]]<-c(1:nrow(output))
	dimnames(output)[[2]]<-c("Predictors","ImportanceScore")
	covariates.num <- length(output[,1]) - length(haploname.freq.gene[,3])
	group.haplo.freq <- matrix(as.numeric(haploname.freq.gene[,2:(1+groupnum)]),ncol=groupnum)
	if(groupnum>1)
	{
		dimnames(group.haplo.freq)[[2]]<-unique(hwe.group)
	}
	#dimnames(group.haplo.freq)[[2]]<-paste(paste("haplo.freq(",unique(hwe.group),sep=""),")",sep="")
	if(groupnum>1)
	{
		print(data.frame(sourceGene=c(haploname.freq.gene[,2+groupnum],rep("NA",covariates.num)),
			haplotype=output[,1],
			importanceScore=round(as.numeric(output[,2]),digits=3),
			haplo.freq=rbind(round(group.haplo.freq,digits=3),matrix(rep("<NA>",covariates.num*groupnum),nrow=covariates.num))))
		result <- list(score=data.frame(sourceGene=c(haploname.freq.gene[,2+groupnum],rep("NA",covariates.num)),
		haplotype=output[,1],
		importanceScore=as.numeric(output[,2]),
		haplo.freq=rbind(group.haplo.freq,matrix(rep("<NA>",covariates.num*groupnum),nrow=covariates.num))))
	}else
	{
		print(data.frame(sourceGene=c(haploname.freq.gene[,2+groupnum],rep("NA",covariates.num)),
			haplotype=output[,1],
			importanceScore=round(as.numeric(output[,2]),digits=3),
			haplo.freq=rbind(round(group.haplo.freq,digits=3),rep("<NA>",covariates.num*groupnum))))
		result <- list(score=data.frame(sourceGene=c(haploname.freq.gene[,2+groupnum],rep("NA",covariates.num)),
		haplotype=output[,1],
		importanceScore=as.numeric(output[,2]),
		haplo.freq=rbind(group.haplo.freq,rep("<NA>",covariates.num*groupnum))))

	}


	
}
	

