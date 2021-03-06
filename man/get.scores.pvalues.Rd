\name{get.scores.pvalues}
\alias{get.scores.pvalues}
\title{
function to compute scores and p-values
}
\description{
For each tested SNP and each parameter in the model, computes scores by summing family scores over all families and computes the corresponding p-values.  P-values of global and joint tests are also computed.
}
\usage{
get.scores.pvalues(test, joint.tests)
}
\arguments{
  \item{test}{
object returned by \code{fat2Lpoly.withinR}.
}
  \item{joint.tests}{
list of vectors of numbers between 1 and the total number of parameters in the \code{design} function. Each vector gives parameter indices to test the corresponding parameters jointly.
}
}
\value{
data frame of p-values for all the tested SNPs, for the global test (all parameters tested jointly), the individual tests and other joint tests specified by the argument \code{joint.tests}.  The p-values are obtained from scores summed over all families.  These scores of individual tests are also included in this data frame.
}
\author{
Alexandre Bureau and Jordie Croteau
}
\seealso{
\code{\link{fat2Lpoly}, \link{fat2Lpoly.withinR}}
}
\examples{
data(fat2Lpoly.allSNPs)
			   
joint.tests=list(c(2,5),c(3,4))

get.scores.pvalues(fat2Lpoly.allSNPs, joint.tests)	
			   
#    snp.cond  snp.test global_p params.joint_2-5_p params.joint_3-4_p param_1_score 
# 1 snp4.loc1 snp2.loc2 5.80e-03           7.12e-01           0.000954         0.449 
# 2 snp4.loc1 snp4.loc2 2.14e-07           1.24e-05           0.000954         0.449 
# 3 snp4.loc1 snp5.loc2 1.14e-03           1.44e-01           0.000954         0.449 
# 4 snp4.loc1 snp6.loc2 5.59e-04           3.84e-02           0.000954         0.449 
# 5 snp4.loc1 snp8.loc2 1.15e-03           1.55e-01           0.000954         0.449 
# param_2_score param_3_score param_4_score param_5_score param_1_p param_2_p
#         0.333        -1.427         3.638         0.733     0.653     0.739
#         0.890        -1.427         3.638         4.612     0.653     0.373		
#         0.776        -1.427         3.638         1.785     0.653     0.438
#        -0.082        -1.427         3.638         2.553     0.653     0.934
#         0.869        -1.427         3.638         1.695     0.653     0.385		
#   param_3_p param_4_p param_5_p
# 1     0.154  0.000275  0.464000
# 2     0.154  0.000275  0.000004
# 3     0.154  0.000275  0.074200
# 4     0.154  0.000275  0.010700
# 5     0.154  0.000275  0.090100	
}
