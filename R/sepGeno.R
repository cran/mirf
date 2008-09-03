#***************************************************************************************************
#Author:Yimin Wu
#Function that takes a genotype matrix that has two letters in each 
#column (with no separation or one symbol as separation such as "/") and makes one column per letter.   
#
#**************************************************************************************************** 
#Input 
#inputGeno: a genotype matrix that has two letters in each column (with no separation or one symbol 
#as separation such as "/") 
#
#*****************************************************************************************
#Output 
#geno: a genotype matrix that has one letter in each column
#SNPnames: SNP names from the input genotype matrix
#
#*****************************************************************************************	
sepGeno <- function(inputGeno)
{
	if(is.null(dimnames(inputGeno)[[2]]))
	{
		SNPnames = NULL
	}else{
		SNPnames = dimnames(inputGeno)[[2]]
	}
	geno = c()
	genoNames = c()
	for(i in 1:dim(inputGeno)[[2]])
	{
		genoNames = cbind(genoNames, paste(dimnames(inputGeno)[[2]][i],".a1",sep=""),paste(dimnames(inputGeno)[[2]][i],".a2",sep=""))
		if(nchar(inputGeno[1,i])==3)
		{
			geno = cbind(geno, substr(inputGeno[,i],1,1),substr(inputGeno[,i],3,3)) 
		}else
		{
			geno = cbind(geno, substr(inputGeno[,i],1,1),substr(inputGeno[,i],2,2)) 
		}
	}
	dimnames(geno)[[2]] <- genoNames;
	dimnames(geno)[[1]] <- dimnames(inputGeno)[[1]]
	result <- list(geno=geno,SNPnames=SNPnames)
	return(result)
}

