\name{mirf}
\alias{mirf}
\title{
MULTIPLE IMPUTATION AND RANDOM FORESTS FOR UNOBSERVABLE PHASE, HIGH-DIMENSIONAL DATA
}
\description{
mirf combines multiple imputation and random forests to characterize haplotype-trait associations.  This approach involves estimation of allelic phase, imputing data according to these estimates and then combining the results of random forests across multiply imputed datasets.  The method depends on two existing R packages: haplo.em and randomForest.  Any parameters in these two packages can also be set in this function.
}

\usage{
mirf(geno, y, gene.column=NULL, gene.names=NULL, SNPnames=NULL,
M=10, hwe.group=NULL, covariates=NULL, residualize=NULL, 
miss.val=c(0, NA, "NA"), weight=NULL, control=haplo.em.control(), ...)
}

\arguments{
\item{geno}{
a data frame or matrix of alleles, such that each locus has a pair of adjacent columns of alleles, and the order of columns corresponds to the order of loci on a chromosome. If there are C loci, then ncol(geno) = 2*C.  Rows represent the alleles for each subject. 
}

\item{y}{
a vector of responses for the subjects, length(y) = nrow(geno).  Both categorical and continuous y are acceptable. (NOTE: missing values are not allowed. If y is categorical, the argument "residualize" is not permitted)
}

\item{gene.column}{
a vector with kth element equal to the number of columns in the geno matrix for gene k. If there are K genes and C loci, then length(gene.column) = K, and sum(gene.column) = ncol(geno) = 2*C. If user doesn't input this variable, then we assume all columns in the geno matrix belong to one gene. In another word, we set the default value to c(2*C).
}

\item{gene.names}{
vector of names for genes. The length of gene.names is equal to the length of gene.column
}

\item{SNPnames}{
vector of names for SNP.
}

\item{M}{
number of haplotype-imputations. The default value is set to 10. This value may be set higher for more accuracy, and is recommended for small sample sizes. (NOTE: increasing this value will increase computation time) 
}

\item{hwe.group}{
a vector containing a categorical variable which defines grouping according to the assumption of Hardy Weinberg equilibrium (HWE).  HWE is assumed within each level of hwe.group.  The default value is set to NULL.  If specified, haplotype frequency estimates are obtained within each group separately.
}

\item{covariates}{
a matrix of covariates that are predictors of the response.  Default set to NULL.  There are two options for dealing with covariates: the first is to include them in the random forest analysis and the second one is to residualize the response (y), see residualize option below. (NOTE: the missing data in the covariates are omitted by default.)
}

\item{residualize}{
a vector indicating which columns in the covariates matrix are to be used in residualizing the response. Default set to NULL. Example: if ncol(covariates) = 3 and residualize = c(1,3), then the first and the third columns in the 
covariates matrix are used to residualize y, and the second column is included in randomForest analysis together with haplotype predictors. If residualize = NULL, then all covriates are included in the randomForest analysis.
}

\item{miss.val}{This is a paramter in "haplo.em". Please refer to the documentation of "haplo.em".
}
\item{weight}{This is a paramter in "haplo.em". Please refer to the documentation of "haplo.em".
}
\item{control}{This is a paramter in "haplo.em". Please refer to the documentation of "haplo.em".
}
\item{...}{You can pass all possible parameters defined in randomForest method. Please refer to the documentation of "randomForest".
}
}
\note{You can also set all remaining parameters accepted by haplo.em and randomForest. Please check their documentation for more information.
}
\value{ 
An object of class mirf which has the following component:
\item{score}{
it should read as a data frame with four(or more) columns. The first column is gene name to which the haplotype predictor corresponds. The second column is imputed haplotypes and covariates included as predictors in randomForest analysis. The third column contains importance scores for the predictors, averaged over imputations taking into account the within- and between-imputation variance. Importance score=NAN implies that the sample size was too small to estimate the average importance score for that predictor.  The remaining columns are the estimated haplotype frequencies from haplo.em for different hwe groups.
}}

\author{
Yimin Wu, B. Aletta S. Nonyane and Andrea S. Foulkes
}

\references{
B.A.S. Nonyane and A.S. Foulkes (2007) Multiple imputation and random forests (MIRF) for unobservable high-dimensional data. 
The International Journal of Biostatistics 3(1): Article 12
}
         
\seealso{
\code{\link[haplo.stats]{haplo.em}}, \code{\link[randomForest]{randomForest}}
}

\examples{\dontrun{
library(mirf)
data(FMSmirfpckg)

# ACTN3 and RESISTIN are genes with 4 and 6 SNPs, respectively. 
genotype <- FMSmirfpckg[,c(2:11)]
gene.names <- c("actn3", "resistin")

# Make a genotype matrix with one letter per column. 
# See help(sepGeno) for more detail.
genoSingleColsResult <- sepGeno(genotype)
genoSingleCols <- genoSingleColsResult$geno
SNPnames <- genoSingleColsResult$SNPnames

# Now ACTN3 gene has 8 columns and RESISTIN has 12 columns
gene.column <- c(8,12) 
trait <- FMSmirfpckg$"NDRM.CH" 

# Assuming HWE for entire cohort and no covariates
# Note: this takes several minutes to run
mirf(geno=genoSingleCols, y=trait, gene.column=gene.column, 
gene.names = gene.names, SNPnames=SNPnames, M=4) 

# Specifying groups within which HWE is expected to hold 
# and specifying covariates to include as predictors 
hwe.group <- FMSmirfpckg$"Race"

# HWE assumed within groups defined by RACE
covariates<-cbind(FMSmirfpckg$"Center",FMSmirfpckg$"Gender",FMSmirfpckg$"Age")
dimnames(covariates)[[2]] <- c("Center","Gender","Age")

mirf(geno=genoSingleCols, y=trait, gene.column=gene.column, 
gene.names = gene.names, SNPnames=SNPnames, M=4, hwe.group=hwe.group, covariates=covariates)

# Residualizing by Center
residualize <-c(1)
mirf(geno=genoSingleCols, y=trait, gene.column=gene.column, 
gene.names = gene.names, SNPnames=SNPnames, M=4, hwe.group=hwe.group, 
covariates=covariates, residualize=residualize)
}}

\keyword{models}          
