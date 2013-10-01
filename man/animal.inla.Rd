\name{animal.inla}
\alias{animal.inla}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fixed version of animal.inla() function from AnimalINLA package
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Updated some INLA commands. Namely, changed every occurrence of
control.data = list(hyper = list(theta = ...))
by
control.family = list(hyper = list(prec = ...)),
inla.marginal.transform by inla.tmarginal
and 
inla.expectation by inla.emarginal

}
\usage{
animal.inla(response, fixed, random, genetic, Ainverse, type.data, data, standardize = FALSE, E = NULL, lambda = NULL, sigma.e = FALSE, Ntrials = NULL, linear.comb = FALSE, linear.comb.name = NULL, verbose = FALSE, dic = FALSE, only.hyperparam = FALSE)
}
\references{
%% ~put references to the literature/web site here ~
http://www.r-inla.org/related-projects/animalinla
}
\author{
%%  ~~who you are~~
Facundo Muñoz <facundo.munoz@orleans.inra.fr>
}
