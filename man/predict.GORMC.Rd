\name{predict.GORMC}
\alias{predict.GORMC}
\title{
Predict cure rate and survival curve for a new subject.
}
\description{
Based on the predicted model, cure rate and survival probabilities are calculated for a new individual with specific 
covariate values. The new.x specifies the covariate vector of the survival part and the new.z specifies that in the 
cure rate part. If no values are specified, zeros will be used and the baseline survival probabilities are returned.
}
\usage{
\method{predict}{GORMC}(object, ...)
}
\arguments{
  \item{object}{
An object from the function GORMC.
}
  \item{\dots}{
Other argumetns including new.x and new.z. See details.
}
}
\details{
For a new subject, new.x and new.z correspond to the covariate vector in the survival part and the cure rate part. The 
first input of new.z is always 1 corresponding to the intercept. Another argument len defines the number of equally 
spaced time points for calculating the survival curve, default is 100.
}
\value{
\item{CureRate}{Estimated cure rate of the new subject.}
\item{Survival}{A data frame including the time points and survival probabilities.}
}

