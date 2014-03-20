myAnimalINLAfixes
=================

Fixed version of animal.inla() function from AnimalINLA package

Updated some INLA syntax. Namely, changed every occurrence of
- control.data = list(hyper = list(theta = ...))
- inla.marginal.transform
- inla.expectation

by
- control.family = list(hyper = list(prec = ...)),
- inla.tmarginal
- inla.emarginal

respectively

## References
http://www.r-inla.org/related-projects/animalinla
