\name{GORMC}
\alias{GORMC}
\title{
Fitting the GORMC model with interval censored data.
}
\description{
The Generalized Odds Rate Mixture Cure model is fitted for interval censored survival data. The EM algorithm facilitated 
by a gamma-poisson data augmentation is applied for estimating the coefficients in both the cure rate part and the survival
part. The covariance matrix has closed forms based on the Louis method. 
}
\usage{
GORMC(survfun = formula(data), curefun = formula(data), data = parent.frame(), 
r = 0, n.int = 5, order = 3, max.iter = 1000, cov.rate = 0.001)
}
\arguments{
  \item{survfun}{
A formula for the survival part in the GORMC model, defined using the Surv function and type=``interval2".
}
  \item{curefun}{
The formula of predictors of the cure rate part in the GORMC model.
}
  \item{data}{
The interval censored sorvival data, including the left and right end points of the time intervals and the covariates for
 the cure rate part and the survival part. If a subject is left(right) censored, the left(right) end point of the subject 
 should be defined as ``NA", see example. 
}
  \item{r}{
The transformation parameter in the GORMC model, should be greater than or equal to 0. r=0 refers to the PHMC model and
  r=1 refers to the POMC model. The default is 0.
}
  \item{n.int}{
Number of interior knots of the splines. Default is 5.
}
  \item{order}{
Order of the spline basis functions. Default is 3, i.e. the cubic splines.
}
  \item{max.iter}{
The maximum number of interations for the EM algorithm. Default is 1000.
}
  \item{cov.rate}{
The bound for convergence of the algorithm, which defined as the difference between the log-likelihood values of two 
consecutive iterations smaller than this value. Default is 0.001.
}
}
\details{
The formula defined for ``survfun" is based on the Surv() function, where the left and right end points of the time interval
 are included and the type is equal to ``interval2". The left(right) end points of left(right) censored individuals should
 be defined as ``NA" in the data frame before running the function.\\ 
The transformation parameter r is a nonnegative number corresponding to a specific model in the GORMC family of models.
Special cases include the PHMC model(r=0) and the POMC model(r=1). Other positive numbers can also be specified. The 
grid search method is suggested to find the best model in practice. That is, try a sequence of r values and choose the
one with the greatest log-likelihood value. 
}
\value{
\item{ParEst}{A list includes the estimated coefficients (Eta,Beta,gl), the whole hessian matrix (Hessian), AIC, and the 
log-likelihood value(loglik).}
\item{ParVcov}{The estimated covariance matrix of the coefficients Eta and Beta.}
}
\references{
Zhou, J., Zhang, J. and Lu, W. (2016+). EM Algorithm for the Generalized Odds Rate Mixture Cure Model with Interval 
Censored Data. Submitted to Journal of Computational and Graphical Statistics.
}

\note{
The estimated hessian matrix can be very large and sometimes not invertable. In which case, we try the QR decomposition,
g-inverse or even numerical methods to get the covariance matrix. Different values of hess in the ParVcov indicating the
 different cases.
hess=0:the hessian matrix is invertable;
hess=1:the QR decomposition is applied to solve the hessian matrix;
hess=2:the g-inverce is applied to the hessian matrix;
hess=3:the hessian matrix is obtained from numerical methods.
The variance estimates may be unreliable for the cases when hess>0.
}

\examples{
data(Hemophilia)
head(Hemophilia)
# Set Left/Right Interval End Points as NA
Hemophilia$L[Hemophilia$d1==1]<-Hemophilia$R[Hemophilia$d3==1]<-NA

# Fit PHMC Model (r=0)
fit<-GORMC(survfun=Surv(L,R)~Low+Medium+High,curefun=~Low+Medium+High,
data=Hemophilia,r=0)
summary(fit)

# Fit POMC Model (r=1)
# fit<-GORMC(survfun=Surv(L,R)~Low+Medium+High,curefun=~Low+Medium+High,
# data=Hemophilia,r=1)
# summary(fit)

# Predict Cure Rate and Survival Curve for a New Individual
# Specify coveriate vectors for new.z and new.x
pred1<-predict(fit,new.z=c(1,0,0,0),new.x=c(0,0,0))
pred2<-predict(fit,new.z=c(1,1,0,0),new.x=c(1,0,0))
pred3<-predict(fit,new.z=c(1,0,1,0),new.x=c(0,1,0))
pred4<-predict(fit,new.z=c(1,0,0,1),new.x=c(0,0,1))

# Obtain Cure Rates
pred1$CureRate
pred2$CureRate
pred3$CureRate
pred4$CureRate

# Plot the Survival Curves
plot(pred1,xlab="Time",ylab="Survival Probability",ylim=c(0,1))
lines(pred2$Survival,col=2)
lines(pred3$Survival,col=3)
lines(pred4$Survival,col=4)
legend(0,0.3,c("None","Low","Medium","High"),lty=1,col=1:4)

# Not run: Grid Search r
# rr<-seq(0,2,0.2)
# logl<-numeric()
# for(i in 1:length(rr)){
# fit<-GORMC(survfun=Surv(L,R)~Low+Medium+High,curefun=~Low+Medium+High,
# data=Hemophilia,r=rr[i])
# logl[i]<-fit$ParEst$loglik
# }
# plot(rr,logl,type="l",xlab="r",ylab="Log-likelihood")
}
\keyword{ GORMC model }
\keyword{ Interval censoring }
