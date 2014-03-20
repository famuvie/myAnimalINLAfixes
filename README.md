myAnimalINLAfixes
=================

### Update: March 20 2014
These fixes have been incorporated into AnimalINLA v1.3,
which makes this package unnecessary.
See reference below.


### Fixed version of animal.inla() function from AnimalINLA package

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
