### bayes factor caculator 1.0 ###

M=100000
nh=2
nf=2

bayes_factor_caculator <- function(old, new) {

lkd.model1=function(y,n,lambda){ return(exp(-n*lambda+y*log(lambda))) }
lkd.model2=function(y1,n1,y2,n2,lambda1,lambda2){ return(exp(-n1*lambda1+y1*log(lambda1)-n2*lambda2+y2*log(lambda2))) }

BF_MC=function(a,b,y1,n1,y2,n2,M){
    lambda1=rgamma(M,a,b)
  m1=cumsum(lkd.model1(y1+y2,n1+n2,lambda1))/(1:M)
  lambda2.1=rgamma(M,a,b)
  lambda2.2=rgamma(M,a,b)
    m2=cumsum(lkd.model2(y1,n1,y2,n2,lambda2.1,lambda2.2))/(1:M)
    return(m2/m1)
  }

con <- file(old,"rt")
line <- readLines(con,n=1)

while (length(line) > 0){
  line <- unlist(strsplit(line,split = "\t"))
  qh <- as.numeric(line[5])
  qf <- as.numeric(line[4])
  set.seed(1)
  bayes <- BF_MC(10,1,qh,nh,qf,nf,M)[M] 
  newline=t(c(line[1],line[2],line[3],round(bayes,5)))
  write.table(newline, new, col.names = F, row.names = F, sep = '\t', quote=F, append =T)
  line <- readLines(con, n = 1)
}

close(con)

# BF_MC(10,1,qh,nh,qf,nf,M)[M]
# qh: total number of mapped reads for heavy-digestion
# nh: number of replicates for heavy-digestion
# qf: total number of mapped reads for light-digestion
# nf: number of replicates for light-digestion

}

args=commandArgs(T)
bayes_factor_caculator(args[1], paste0(args[1], '.bayes'))
