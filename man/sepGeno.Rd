\name{sepGeno}
\alias{sepGeno}
\title{
CHANGE A GENOTYPE MATRIX THAT HAS TWO LETTERS IN EACH 
COLUMN  TO A NEW GENOTYPE MATRIX WITH ONE COLUMN PER LETTER
}
\description{
This function is a preprocess function that takes a geno matrix that has two letters in each column (with no separation or any one character as separation) and return a new geno matrix with one column per letter. This is useful because mirf and haplo.em only take a geno matrix with one column per letter. 
}

\usage{
sepGeno(inputGeno)
}

\arguments{
\item{inputGeno}{
a data frame or matrix of alleles, each colunm has two letters (with no separation or one symbol as separation such as "/"). Rows represent the alleles for each subject. 
}}

\value{ 
An object of class sepGeno which has the following component:
\item{geno}{a new geno matrix with one column per letter, which can be used in mirf or haplo.em}

\item{SNPnames}{
a vector of names for SNPs. If it is NULL, that means the inputGeno matrix doesn't have column names. This vector can be input as SNPnames in mirf or as locus.label in haplo.em}
}

\author{
Yimin Wu, B. Aletta S. Nonyane and Andrea S. Foulkes
}

\references{
B.A.S. Nonyane and A.S. Foulkes (2007) Multiple imputation and random forests (MIRF) for unobservable high-dimensional data. 
The International Journal of Biostatistics 3(1): Article 12
}
         
\seealso{
\code{\link{mirf}},\code{\link[haplo.stats]{haplo.em}}, \code{\link[randomForest]{randomForest}}
}

\examples{
library(mirf)
data(FMSmirfpckg)

inputGeno <- FMSmirfpckg[1:4,c(2:5)]

genoSingleCols <- sepGeno(inputGeno)
genoSingleCols$geno
genoSingleCols$SNPnames  
}
\keyword{models}          
