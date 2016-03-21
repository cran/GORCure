\name{GORCure-package}
\alias{GORCure-package}
\alias{GORCure}
\docType{package}
\title{
Fit Generalized Odds Rate Mixture Cure Model with Interval Censored Data
}
\description{
Generalized Odds Rate Mixture Cure (GORMC) model is a flexible cure model which includes the Proportional Hazards Mixture
Cure (PHMC) model and the Proportional Odds Mixture Cure (POMC) model as special cases. When the survival data is interval
censored, this package provides an efficient solution for the GORMC model based on a gamma-possion data augmentation.
}
\details{
\tabular{ll}{
Package: \tab GORCure\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-03-07\cr
License: \tab GPL-2\cr
}
The main function is GORMC, which produces the estimates for coefficients and the covariance matrix. The summary function 
can be used to get the test results of the coefficients. The predict function can be used to obtain the estimated cure rate
and survival curve for an individual with specified covariates for the cure rate and survival parts separatly. The plot 
function can be applied to generate the predicted survival curve for that individual.
}
\author{
Jie Zhou, Jiajia Zhang, Wenbin Lu

Maintainer: Jie Zhou <zhoujie02569@gmail.com>
}
\references{
Zhou, J., Zhang, J. and Lu, W. (2016+). EM Algorithm for the Generalized Odds Rate Mixture Cure Model with Interval 
Censored Data. Submitted to Journal of Computational and Graphical Statistics.
}
\keyword{ package, linear transformation model, mixture cure model, interval censoring }


