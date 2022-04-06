library(glmnet)

n<-c(100,200,300,400) #number of samples
p<-c(100,400,900,1600) #number of variables
p_0<-c(3,10,20) #number of true variables

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  (matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL)) 
}

data_gen_4<-function(n,p,p_0)
{
  ma<-matrix(data=0.3,nrow=p,ncol = p)
  diag(ma)<-rep(1,p)
  X<-rmvn(n,rep(0,p),ma)
  beta<-c(rep(3,p_0),rep(0,p-p_0))
  ep<-rcauchy(n,0,1)
  Y= exp(4+ 0.5*(X%*%beta))  + ep
  return(list(X,beta,Y))
}



lambda_rl<-function(p,n){
  return(0.3*sqrt(log(p)/n))
}

weights_arl<-function(lambda,beta)
{
  w<- sapply(1:length(beta), function(x) ifelse(abs(beta[x]) > 0.1*lambda,0.1*lambda/abs(beta[x]),1/abs(beta[x])))
  return(w)
}

lambda_lad<-function(p,n){
  return(1.5*sqrt(log(p)/n))
}

fp_rl_1<-rep(0,200)
fn_rl_1<-rep(0,200)
tp_rl_1<-rep(0,200)
nmp_1<-matrix(0,nrow = 3,ncol = 4)
fdr_1<-matrix(0,nrow = 3,ncol = 4)
power_1<-matrix(0,nrow = 3,ncol = 4)

fp_arl_1<-rep(0,200)
fn_arl_1<-rep(0,200)
tp_arl_1<-rep(0,200)
nmp_arl_1<-matrix(0,nrow = 3,ncol = 4)
fdr_arl_1<-matrix(0,nrow = 3,ncol = 4)
power_arl_1<-matrix(0,nrow = 3,ncol = 4)


fp_thrl_1<-rep(0,200)
fn_thrl_1<-rep(0,200)
tp_thrl_1<-rep(0,200)
nmp_thrl_1<-matrix(0,nrow = 3,ncol = 4)
fdr_thrl_1<-matrix(0,nrow = 3,ncol = 4)
power_thrl_1<-matrix(0,nrow = 3,ncol = 4)

# fp_lad_1<-rep(0,200)
# fn_lad_1<-rep(0,200)
# tp_lad_1<-rep(0,200)
# nmp_lad_1<-matrix(0,nrow = 3,ncol = 4)
# fdr_lad_1<-matrix(0,nrow = 3,ncol = 4)
# power_lad_1<-matrix(0,nrow = 3,ncol = 4)

fp_lasso_1<-rep(0,200)
fn_lasso_1<-rep(0,200)
tp_lasso_1<-rep(0,200)

nmp_lasso_1<-matrix(0,nrow = 3,ncol = 4)
fdr_lasso_1<-matrix(0,nrow = 3,ncol = 4)
power_lasso_1<-matrix(0,nrow = 3,ncol = 4)


for (i in 1:3) 
{
  for(k in 1:4)
  {
    print(k)
    for (j in 1:200)
    {
      X=matrix(data=NA,n[k],n[k])
      beta=rep(NA,p[k])
      Y=rep(NA,n[k])
      c(X,beta,Y):=data_gen_4(n[k],p[k],p_0[i])
      
      model_rl<-glmnet(X,(rank(Y)/n[k]-0.5),family = "gaussian",lambda = lambda_rl(p[k],n[k]),intercept = FALSE)
      beta_est_rl<-model_rl$beta
      
      fp_rl_1[j]<-sum(beta_est_rl[(p_0[i]+1):p[k]]!=beta[(p_0[i]+1):p[k]])
      fn_rl_1[j]<-sum(beta_est_rl[(1:p_0[i])]==0)
      tp_rl_1[j]<-sum(beta_est_rl[(1:p_0[i])]!=0)
      
      model_arl<-glmnet(X,(rank(Y)/n[k]-0.5),family = "gaussian",lambda = 2*lambda_rl(p[k],n[k]),intercept = FALSE,penalty.factor =   weights_arl(lambda=lambda_rl(p[k],n[k]),as.array(beta_est_rl)))
      beta_est_arl<-model_arl$beta
      
      fp_arl_1[j]<-sum(beta_est_arl[(p_0[i]+1):p[k]]!=beta[(p_0[i]+1):p[k]])
      fn_arl_1[j]<-sum(beta_est_arl[(1:p_0[i])]==0)
      tp_arl_1[j]<-sum(beta_est_arl[(1:p_0[i])]!=0)
      
      model_thrl<-cv.glmnet(X,(rank(Y)/n[k]-0.5),family = "gaussian")
      beta_est_thrl<-model_thrl$glmnet.fit$beta[,model_thrl$index[1]]
      ordered_beta_est_thrl<-beta_est_thrl[order(beta_est_thrl,decreasing = TRUE)]
      beta_est_thrl[!(order(beta_est_thrl,decreasing = TRUE) <= sum(beta_est_arl!=0))]<-0
      
      fp_thrl_1[j]<-sum(beta_est_thrl[(p_0[i]+1):p[k]]!=beta[(p_0[i]+1):p[k]])
      fn_thrl_1[j]<-sum(beta_est_thrl[(1:p_0[i])]==0)
      tp_thrl_1[j]<-sum(beta_est_thrl[(1:p_0[i])]!=0)
      
      
      # model_lad<-LADlasso(X,Y,beta.ini=LAD(X, Y),lambda = lambda_lad(p[k],n[k]))
      # beta_est_lad<-model_lad$beta
      # fp_lad_1[j]<-sum(beta_est_lad[(p_0[i]+1):p[k]]!=beta[(p_0[i]+1):p[k]])
      # fn_lad_1[j]<-sum(beta_est_lad[(1:p_0[i])]==0)
      # tp_lad_1[j]<-sum(beta_est_lad[(1:p_0[i])]!=0)
      
      model_lasso<-cv.glmnet(X,Y,family="gaussian")
      beta_est_lasso<-model_lasso$glmnet.fit$beta[,model_lasso$index[1]]
      fp_lasso_1[j]<-sum(beta_est_lasso[(p_0[i]+1):p[k]]!=beta[(p_0[i]+1):p[k]])
      fn_lasso_1[j]<-sum(beta_est_lasso[(1:p_0[i])]==0)
      tp_lasso_1[j]<-sum(beta_est_lasso[(1:p_0[i])]!=0)
      
    }
    nmp_1[i,k]<-mean(fp_rl_1+fn_rl_1)
    fdr_1[i,k]<-mean(fp_rl_1/max((fp_rl_1+tp_rl_1),1))
    power_1[i,k]<-mean(tp_rl_1/p_0[i])
    
    nmp_arl_1[i,k]<-mean(fp_arl_1+fn_arl_1)
    fdr_arl_1[i,k]<-mean(fp_arl_1/max((fp_arl_1+tp_arl_1),1))
    power_arl_1[i,k]<-mean(tp_arl_1/p_0[i])
    
    nmp_thrl_1[i,k]<-mean(fp_thrl_1+fn_thrl_1)
    fdr_thrl_1[i,k]<-mean(fp_thrl_1/max((fp_thrl_1+tp_thrl_1),1))
    power_thrl_1[i,k]<-mean(tp_thrl_1/p_0[i])
    
    # nmp_lad_1[i,k]<-mean(fp_lad_1+fn_lad_1)
    # fdr_lad_1[i,k]<-mean(fp_lad_1/max((fp_lad_1+tp_lad_1),1))
    # power_lad_1[i,k]<-mean(tp_lad_1/p_0[i])
    
    nmp_lasso_1[i,k]<-mean(fp_lasso_1+fn_lasso_1)
    fdr_lasso_1[i,k]<-mean(fp_lasso_1/max((fp_lasso_1+tp_lasso_1),1))
    power_lasso_1[i,k]<-mean(tp_lasso_1/p_0[i])
  }
}

sink(file = "output-4.txt")
nmp_1
fdr_1
power_1
nmp_arl_1
fdr_arl_1
power_arl_1
nmp_thrl_1
fdr_thrl_1
power_thrl_1
# nmp_lad_1
# fdr_lad_1
# power_lad_1
nmp_lasso_1
fdr_lasso_1
power_lasso_1
sink(file = NULL)

save(nmp_1,fdr_1,power_1,nmp_arl_1,fdr_arl_1,power_arl_1,nmp_thrl_1,fdr_thrl_1,power_thrl_1,nmp_lasso_1,fdr_lasso_1,power_lasso_1,file="output-4.RData")
