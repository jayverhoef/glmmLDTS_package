\name{glmmLDTS}
\alias{glmmLDTS}
\title{
	Fitting Generalized Linear Models for Large Data sets with Temporal 
	Autocorrelation using pseudo data and quasi-likelihood
}
\description{
	This function takes data frames with a time column and fits generalized linear
	models with temporal autocorrelation using pseudo data and quasi-likelihood.
}
\usage{
	glmmLDTS(fixed.formula,
		random.formula,
		data,
		timecol,
		EstMeth = "REML",
		distribution = "binomial",
		trialscol = NULL,
		group.vec = NULL,
		ridge.reg = NULL,
		lambda = 0)
}
\arguments{
	\item{fixed.formula}{
		an object of class "formula" (or one that can be coerced to that class):
		a symbolic description of the model to be fitted. The details of model
		specification are given under 'Details'.
	}
	\item{random.formula}{
		an object of class "formula" (or one that can be coerced to that class):
		a symbolic description of the model to be fitted. The details of model
		specification are given under 'Details'.
    }
	\item{data}{
		a data.frame containing the data.
	}
	\item{timecol}{
		column containing the time of each event.  They must be
		ordered sequentially within subject.
	}
	\item{EstMeth}{
		Estimation method; either "ML" for maximum likelihood, or "REML" for 
		restricted maximum likelihood (default).
	}
	\item{distribution}{
		the error distribution and link function to be used in the
		model. This is a character string that is either "Gaussian", "Poisson", or
		"Binomial."
	}
	\item{trialscol}{
		name of the variable in the data data.frame that contains
		the sample size when a binomial distribution is used.  If NULL, 
		a sample size of 1 is assumed, and the response variable must be 
		binary (0 or 1).
	}
	\item{group.vec}{
		The vector containing grouping.
	}
	\item{ridge.reg}{
    	The use of ridge regression.  If NULL (default), no ridge regression
		is used.  Other options are "local" and "global."
	}
	\item{lambda}{
		Parameter for the ridge regression
	}
}
\details{
Models for glmssn are specified symbolically, similar to \code{\link{lm}}
and other models in R.  A typical model has the form response ~ terms
where response is the (numeric) response vector and terms is a series of
fixed effect linear predictors for the response.  A terms specification of
the form first + second indicates all the terms in first together with
all the terms in second with duplicates removed.  A specification of the
form first:second indicates the set of terms obtained by taking the
interactions of all terms in first with all terms in second.  The
specification first*second indicates the cross of first and second.
This is the same as first + second + first:second. See
\code{\link{model.matrix}} for further details.  The terms in the
formula will be re-ordered so that main effects come first, followed by
the interactions, all second-order, all third-order and so on.  A
formula has an implied intercept term. To remove this use either y ~ x -
1 or y ~ 0 + x.  See \code{\link{formula}} for more details of allowed
formulae.
}
\value{
  \command{glmmLDTS} returns a list with the following items:
	\item{dataset}{
		The original data set.
	}
	\item{WARNINGS}{
		any warnings issued by the program
	}
	\item{fixed.formula}{
   		reproducing the fixed.formula argument
	}
	\item{random.formula}{
   		reproducing the random.formula argument
	}
	\item{sample.size}{
   		total number of records used in the data set
	}
	\item{timecol}{
   		reproducing the timecol argument
	}
	\item{trialscol}{
   		reproducing the trialscol argument
	}
	\item{group.vec}{
   		reproducing the group.vec argument
	}
	\item{ridge.reg}{
   		reproducing the ridge.reg argument
	}
	\item{lambda}{
   		reproducing the lambda argument
	}
	\item{start.time}{
   		starting time of the function, taken from Sys.time()
	}
	\item{end.time}{
   		ending time of the function, taken from Sys.time()
	}
	\item{R.cov.parameters}{
   		estimated R-side covariance parameters
	}
	\item{G.cov.parameters}{
   		estimated G-side covariance parameters
	}
	\item{fixed.effects}{
   		fixed effects table
	}
	\item{random.effects}{
   		estimated random effects
	}
	\item{typeIII.hypoth}{
   		SAS type III test of factors
	}
	\item{fit.table}{
   		data.frame of fitted values
	}
	\item{covb}{
   		estimated covariance matrix of fixed effects
	}
	\item{outer.iterations}{
   		number of outer iterations
	}
	\item{inner.iterations}{
   		number of inner iterations
	}

}

\references{
\cite{Ver Hoef, J.M., London, J.M., and Boveng, P.L. (2010) Fast computing of 
	some generalized linear mixed pseudo-models with temporal autocorrelation. 
	Computational Statistics 25(1): 39 – 55. DOI: 10.1007/s00180-009-0160-1}

}

\author{
Jay Ver Hoef \email{jay.verhoef@noaa.gov}
}

\examples{

library(glmmLDTS)
# use ML for inner optimization to compare to glmmPQL
data(ValidData)
glmmLDTSfitML <- glmmLDTS(HOStatus ~ MFAug15 + MFA15.2 + hour2,
    random.formula = HOStatus ~ SpeNo,
    EstMeth = "ML",
    data = ValidData,
	  timecol = "time.vec", group.vec = "SpeNo")
glmmLDTSfitML$fixed.formula
glmmLDTSfitML$random.formula
glmmLDTSfitML$sample.size
glmmLDTSfitML$timecol
glmmLDTSfitML$trialscol
glmmLDTSfitML$group.vec
glmmLDTSfitML$ridge.reg
glmmLDTSfitML$lambda
glmmLDTSfitML$start.time
glmmLDTSfitML$end.time
glmmLDTSfitML$R.cov.parameters
glmmLDTSfitML$G.cov.parameters
glmmLDTSfitML$typeIII.hypoth
glmmLDTSfitML$fixed.effects
glmmLDTSfitML$random.effects
glmmLDTSfitML$outer.iterations
glmmLDTSfitML$inner.iterations2

# compare to glmmPQL in MASS package
attach(ValidData)
glmmPQLfit <- glmmPQL(HOStatus ~ MFAug15 + MFA15.2 + hour2, random = ~ 1 | SpeNo,
    data = ValidData,
    correlation = corExp(form = ~ time.vec),
    family = binomial)
summary(glmmPQLfit)

# use REML for inner optimization to compare to SAS PROC GLIMMIX (REML is the 
# default inner optimization method).  PROC GLIMMIX can also use ML for inner 
# optimization and be compared to above fits)
glmmLDTSfitREML <- glmmLDTS(HOStatus ~ MFAug15 + MFA15.2 + hour2,
    random.formula = HOStatus ~ SpeNo,
    EstMeth = "REML",
    data = ValidData,
	  timecol = "time.vec", group.vec = "SpeNo")
glmmLDTSfitREML$fixed.formula
glmmLDTSfitREML$random.formula
glmmLDTSfitREML$sample.size
glmmLDTSfitREML$timecol
glmmLDTSfitREML$trialscol
glmmLDTSfitREML$group.vec
glmmLDTSfitREML$ridge.reg
glmmLDTSfitREML$lambda
glmmLDTSfitREML$start.time
glmmLDTSfitREML$end.time
glmmLDTSfitREML$R.cov.parameters
glmmLDTSfitREML$G.cov.parameters
glmmLDTSfitREML$typeIII.hypoth
glmmLDTSfitREML$fixed.effects
glmmLDTSfitREML$random.effects
glmmLDTSfitREML$outer.iterations
glmmLDTSfitREML$inner.iterations2

}

