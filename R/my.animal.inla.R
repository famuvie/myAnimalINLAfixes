animal.inla <- function (response, fixed, random, genetic, Ainverse, type.data, 
          data, standardize = FALSE, E = NULL, lambda = NULL, sigma.e = FALSE, 
          Ntrials = NULL, linear.comb = FALSE, linear.comb.name = NULL, 
          verbose = FALSE, dic = FALSE, only.hyperparam = FALSE) 
{
  if (class(Ainverse) == "ped") {
    if (!is.list(Ainverse)) 
      stop("'Ainverse' is not a list")
    if (!length(Ainverse) == 2) 
      stop("'Ainverse' with class='ped' have not the not the correct number of elements")
  }
  else {
    if (!class(Ainverse) == "dgCMatrix") 
      stop("'Ainverse' is not a list")
  }
  if (class(Ainverse) == "ped") {
    Ainv = Ainverse$Ainverse
    map = Ainverse$map
    Cmatrix = sparseMatrix(i = Ainv[, 1], j = Ainv[, 2], 
                           x = Ainv[, 3])
  }
  else {
    Cmatrix = Ainverse
  }
  if (!(sum(search() == "package:INLA")) == 1) 
    stop("Packages 'INLA' is not found")
  if (class(Ainverse) == "ped") {
    call <- match.call()
    if (!is.null(call$genetic)) {
      names.genetic = eval(call$genetic)
    }
    IndexA = map[, 2][match(eval(parse(text = paste("data$", 
                                                    names.genetic[1], sep = ""))), map[, 1])]
    eval(parse(text = paste("data$", names.genetic[1], "=IndexA")))
  }
  if (type.data == "gaussian" & standardize == TRUE) {
    call <- match.call()
    eval(parse(text = paste("data$", call$response, "=(data$", 
                            call$response, "- mean(data$", call$response, "))/sd(data$", 
                            call$response, ")")))
  }
  else {
    call <- match.call()
    eval(parse(text = paste("data$", call$response, "=data$", 
                            call$response)))
  }
  if (linear.comb == TRUE & is.null(linear.comb.name)) {
    stop("linear.comb.name is missing")
  }
  if (linear.comb == TRUE) {
    call <- match.call()
    if (!is.null(call$genetic)) {
      names.genetic = eval(call$genetic)
    }
    lineardata <- lincomb(linear.comb.name, data, names.genetic)
    lin = lincombinput(lineardata)
  }
  prec.fixed.generic2 = list(initial = 10, fixed = TRUE)
  prec.A = list(param = c(0.5, 0.5), fixed = FALSE)
  prec.e.generic0 = list(param = c(0.5, 0.5), fixed = FALSE)
  prec.Fixed = list(initial = -10, fixed = TRUE)
  prec.Random = list(param = c(1, 0.001), fixed = FALSE)
  if (type.data == "gaussian") {
    call <- match.call()
    formula = paste(call$response, "~")
    if (!is.null(call$fixed)) {
      names.fixed = eval(call$fixed)
      for (i in 1:length(names.fixed)) formula = paste(formula, 
                                                       "f(", names.fixed[i], ",model = \"iid\",constr=TRUE, hyper = list(theta = prec.Fixed))", 
                                                       "+")
    }
    if (!is.null(call$random)) {
      names.random = eval(call$random)
      for (i in 1:length(names.random)) formula = paste(formula, 
                                                        "f(", names.random[i], ",model = \"iid\", constr=TRUE,hyper = list(theta =  prec.Random))", 
                                                        "+")
    }
    if (!is.null(call$genetic)) {
      names.genetic = eval(call$genetic)
      formula = paste(formula, "f(", names.genetic[1], 
                      ",model = \"generic2\",hyper = list(theta1 = prec.A, theta2 = prec.A),constr=TRUE", 
                      paste(",Cmatrix=Cmatrix)"), "+")
    }
    formula = paste(formula, "1")
    formula = as.formula(formula)
    model.animal = inla(formula, family = type.data, data = data, 
                        verbose = verbose, control.compute = list(dic = FALSE), 
                        control.family = list(hyper = list(prec = prec.fixed.generic2)), 
                        only.hyperparam = only.hyperparam)
    sigma.u = inla.tmarginal(function(x) 1/x, eval(parse(text = paste(paste("model.animal$marginals.hyperpar$\n\t\t\"Precision-cmatrix for", 
                                                                                     names.genetic[1]), "\"", sep = ""))))
    h2 = eval(parse(text = paste(paste("model.animal$marginals.hyperpar$\n\t\t\"h2 for", 
                                       names.genetic[1]), "\"", sep = "")))
    breedingvalues = eval(parse(text = paste("model.animal$marginals.random$", 
                                             names.genetic[1])))[-c(1:dim(Cmatrix)[1])]
    summary.breedingvalues = eval(parse(text = paste("model.animal$summary.random$", 
                                                     names.genetic[1])))[-c(1:dim(Cmatrix)[1]), ]
    row.names(summary.breedingvalues) <- NULL
  }
  if (type.data == "gaussian" & dic == TRUE | type.data == 
        "gaussian" & sigma.e == TRUE | type.data == "gaussian" & 
        linear.comb == TRUE | type.data == "gaussian" & dic == 
        TRUE & sigma.e == TRUE & linear.comb == TRUE | type.data == 
        "gaussian" & sigma.e == TRUE & dic == TRUE | type.data == 
        "gaussian" & sigma.e == TRUE & linear.comb == TRUE | 
        type.data == "gaussian" & dic == TRUE & linear.comb == 
        TRUE) {
    call <- match.call()
    formulaG0 = paste(call$response, "~")
    if (!is.null(call$fixed)) {
      names.fixed = eval(call$fixed)
      for (i in 1:length(names.fixed)) formulaG0 = paste(formulaG0, 
                                                         "f(", names.fixed[i], ",model = \"iid\", constr=TRUE,hyper = list(theta = prec.Fixed))", 
                                                         "+")
    }
    if (!is.null(call$random)) {
      names.random = eval(call$random)
      for (i in 1:length(names.random)) formulaG0 = paste(formulaG0, 
                                                          "f(", names.random[i], ",model = \"iid\", constr=TRUE,hyper = list(theta =  prec.Random))", 
                                                          "+")
    }
    if (!is.null(call$genetic)) {
      names.genetic = eval(call$genetic)
      formulaG0 = paste(formulaG0, "f(", names.genetic[1], 
                        ",model = \"generic0\", hyper = list(theta = prec.A),constr=TRUE", 
                        paste(",Cmatrix=Cmatrix)"), "+")
    }
    formulaG0 = paste(formulaG0, "1")
    formulaG0 = as.formula(formulaG0)
    if (!(linear.comb == TRUE)) {
      model.animalG0 = inla(formulaG0, family = type.data, 
                            data = data, verbose = verbose, control.compute = list(dic = dic), 
                            control.family = list(hyper = list(prec = prec.e.generic0)), 
                            only.hyperparam = only.hyperparam)
    }
    if (linear.comb == TRUE) {
      model.animalG0 = inla(formulaG0, family = type.data, 
                            data = data, verbose = verbose, control.compute = list(dic = dic), 
                            control.family = list(hyper = list(prec = prec.e.generic0)), 
                            only.hyperparam = only.hyperparam, lincomb = lin)
      linear.combinations <- model.animalG0$summary.lincomb.derived
      row.names(linear.combinations) = NULL
      linear.combinations$ID = colnames(lineardata)
    }
    if (sigma.e == TRUE) {
      sigmae = inla.tmarginal(function(x) 1/x, 
                                       model.animalG0$marginals.hyperpar$"Precision for the Gaussian observations")
    }
    if (dic == TRUE) {
      dic2 = model.animalG0$dic
    }
  }
  if (type.data == "binomial") {
    call <- match.call()
    formula = paste(call$response, "~")
    if (!is.null(call$fixed)) {
      names.fixed = eval(call$fixed)
      for (i in 1:length(names.fixed)) formula = paste(formula, 
                                                       "f(", names.fixed[i], ",model = \"iid\", constr=TRUE,hyper = list(theta = prec.Fixed))", 
                                                       "+")
    }
    if (!is.null(call$random)) {
      names.random = eval(call$random)
      for (i in 1:length(names.random)) formula = paste(formula, 
                                                        "f(", names.random[i], ",model = \"iid\", constr=TRUE,hyper = list(theta =  prec.Random))", 
                                                        "+")
    }
    if (!is.null(call$genetic)) {
      names.genetic = eval(call$genetic)
      formula = paste(formula, "f(", names.genetic[1], 
                      ",model = \"generic0\", hyper = list(theta = prec.A), constr=TRUE", 
                      paste(",Cmatrix=Cmatrix)"), "+")
    }
    formula = paste(formula, "1")
    formula = as.formula(formula)
    if (!is.null(call$Ntrials)) {
      Ntrials = eval(parse(text = paste("data$", call$Ntrials)))
      if (all(Ntrials == 1)) {
        warning("Binary data, number of trials are 1, is not accurate!!!")
      }
      else {
        Ntrials = eval(parse(text = paste("data$", call$Ntrials)))
      }
    }
    if (!(linear.comb == TRUE)) {
      model.animal = inla(formula, family = type.data, 
                          data = data, verbose = verbose, control.compute = list(dic = dic), 
                          Ntrials = Ntrials, only.hyperparam = only.hyperparam)
    }
    if (linear.comb == TRUE) {
      model.animal = inla(formula, family = type.data, 
                          data = data, verbose = verbose, control.compute = list(dic = dic), 
                          Ntrials = Ntrials, only.hyperparam = only.hyperparam, 
                          lincomb = lin)
      linear.combinations <- model.animal$summary.lincomb.derived
      row.names(linear.combinations) = NULL
      linear.combinations$ID = colnames(lineardata)
    }
    dic2 = model.animal$dic
    sigma.u = inla.tmarginal(function(x) 1/x, eval(parse(text = paste(paste("model.animal$marginals.hyperpar$\n\t\t  \"Precision for", 
                                                                                     names.genetic[1]), "\"", sep = ""))))
    hbin = inla.tmarginal(function(x) 1/x/(1/x + 
                                                      (pi^2)/3), eval(parse(text = paste(paste("model.animal$marginals.hyperpar$\n\t\t \"Precision for", 
                                                                                               names.genetic[1]), "\"", sep = ""))))
    breedingvalues = eval(parse(text = paste("model.animal$marginals.random$", 
                                             names.genetic[1])))
    summary.breedingvalues = eval(parse(text = paste("model.animal$summary.random$", 
                                                     names.genetic[1])))
    row.names(summary.breedingvalues) <- NULL
  }
  if (type.data == "poisson" | type.data == "zeroinflatedpoisson1" | 
        type.data == "zeroinflatedpoisson0") {
    call <- match.call()
    formula = paste(call$response, "~")
    if (!is.null(call$fixed)) {
      names.fixed = eval(call$fixed)
      for (i in 1:length(names.fixed)) formula = paste(formula, 
                                                       "f(", names.fixed[i], ",model = \"iid\",constr=TRUE,hyper = list(theta = prec.Fixed))", 
                                                       "+")
    }
    if (!is.null(call$random)) {
      names.random = eval(call$random)
      for (i in 1:length(names.random)) formula = paste(formula, 
                                                        "f(", names.random[i], ",model = \"iid\", constr=TRUE,hyper = list(theta =  prec.Random))", 
                                                        "+")
    }
    if (!is.null(call$genetic)) {
      names.genetic = eval(call$genetic)
      formula = paste(formula, "f(", names.genetic[1], 
                      ",model = \"generic0\", hyper = list(theta = prec.A), constr=TRUE", 
                      paste(",Cmatrix=Cmatrix)"), "+")
    }
    if (!is.null(call$E)) {
      E = eval(parse(text = paste("data$", call$E)))
    }
    formula = paste(formula, "1")
    formula = as.formula(formula)
    if (!(linear.comb == TRUE)) {
      model.animal = inla(formula, family = type.data, 
                          data = data, verbose = verbose, control.compute = list(dic = dic), 
                          E = E, only.hyperparam = only.hyperparam)
    }
    if (linear.comb == TRUE) {
      model.animal = inla(formula, family = type.data, 
                          data = data, verbose = verbose, control.compute = list(dic = dic), 
                          E = E, only.hyperparam = only.hyperparam, lincomb = lin)
      linear.combinations <- model.animal$summary.lincomb.derived
      row.names(linear.combinations) = NULL
      linear.combinations$ID = colnames(lineardata)
    }
    dic2 = model.animal$dic
    sigma.u = inla.tmarginal(function(x) 1/x, eval(parse(text = paste(paste("model.animal$marginals.hyperpar$\"Precision for", 
                                                                                     names.genetic[1]), "\"", sep = ""))))
    if (!is.null(lambda)) {
      hpois = inla.tmarginal(function(x) 1/x/(1/x + 
                                                         lambda), eval(parse(text = paste(paste("model.animal$marginals.hyperpar$\"Precision for", 
                                                                                                names.genetic[1]), "\"", sep = ""))))
    }
    breedingvalues = eval(parse(text = paste("model.animal$marginals.random$", 
                                             names.genetic[1])))
    summary.breedingvalues = eval(parse(text = paste("model.animal$summary.random$", 
                                                     names.genetic[1])))
    row.names(summary.breedingvalues) <- NULL
  }
  get.variance <- function(prec) {
    expect = inla.emarginal(function(x) 1/x, prec)
    stdev = sqrt((inla.emarginal(function(x) 1/x^2, prec)) - 
                   (inla.emarginal(function(x) 1/x^1, prec))^2)
    q = inla.qmarginal(c(0.025, 0.5, 0.975), inla.tmarginal(function(x) 1/x, 
                                                                     prec))
    summarytab <- cbind(expect, stdev, q[1], q[2], q[3])
    colnames(summarytab) = c("mean", "sd", "0.025quant", 
                             "0.5quant", "0.975quant")
    return(summarytab)
  }
  if (type.data == "gaussian") {
    summarytab = get.variance(prec = eval(parse(text = paste(paste("model.animal$marginals.hyperpar$\n\t\t\t\"Precision-cmatrix for", 
                                                                   names.genetic[1]), "\"", sep = ""))))
    rownames(summarytab) = paste("Variance for ", names.genetic[1], 
                                 sep = "")
    summarytab = rbind(Heritability = eval(parse(text = paste(paste(paste("model.animal$\n\t\t\tsummary.hyperpar[\"h2 for", 
                                                                          names.genetic[1]), "\"", sep = ""), ",c(\"mean\",\n\t\t\"sd\",\"0.025quant\",\"0.5quant\",\"0.975quant\")]"))), 
                       summarytab)
  }
  if (sigma.e == TRUE & type.data == "gaussian") {
    sigmaev = get.variance(eval(parse(text = paste("model.animalG0$marginals.hyperpar$\n\t\t\t\"Precision for the Gaussian observations\""))))
    rownames(sigmaev) = paste("Variance for e")
    summarytab = rbind(summarytab, sigmaev)
  }
  if (type.data == "binomial" | type.data == "poisson" | type.data == 
        "zeroinflatedpoisson1" | type.data == "zeroinflatedpoisson0") {
    summarytab = get.variance(eval(parse(text = paste(paste("model.animal$marginals.hyperpar$\n\t\t\t\"Precision for", 
                                                            names.genetic[1]), "\"", sep = ""))))
    rownames(summarytab) = paste("Variance for ", names.genetic[1], 
                                 sep = "")
  }
  call <- match.call()
  if (!is.null(call$random)) {
    names.random = eval(call$random)
    tab1 <- matrix(NA, ncol = 5, nrow = length(names.random))
    namerow = rep(NA, length(names.random))
    for (i in 1:(length(names.random))) {
      random.var = get.variance(eval(parse(text = paste(paste("model.animal$marginals.hyperpar$\n\t\t\"Precision for", 
                                                              names.random[i]), "\"", sep = ""))))
      tab1[i, ] <- random.var
      namerow[i] = paste("Variance for ", names.random[i], 
                         sep = "")
    }
    rownames(tab1) <- namerow
    summarytab <- rbind(summarytab, tab1)
  }
  output <- list()
  output$formula = formula
  output$inla.call = model.animal$call
  output$Time.used = model.animal$cpu.used
  output$family = model.animal$family
  output$model.covariates = model.animal$model.random
  output$sigma.u = sigma.u
  output$summary.hyperparam = summarytab
  output$Cmatrix = Cmatrix
  output$dataset = data
  output$summary.covariates = model.animal$summary.random
  output$name.covariates = names(model.animal$summary.random)
  output$intercept = model.animal$summary.fixed
  output$prec.fixed.generic2 = paste(list(prec.fixed.generic2))
  output$prec.A = paste(list(prec.A))
  output$prec.e.generic0 = paste(list(prec.e.generic0))
  output$prec.Fixed = paste(list(prec.Fixed))
  output$prec.Random = paste(list(prec.Random))
  call <- match.call()
  if (!is.null(call$fixed)) {
    output = c(output, fixed.names = list(names.fixed))
  }
  if (!is.null(call$random)) {
    output = c(output, random.names = list(names.random))
  }
  if (!is.null(call$genetic)) {
    output = c(output, genetic.names = list(names.genetic))
  }
  if (class(Ainverse) == "ped") {
    names(breedingvalues) = map[, 1]
    summary.breedingvalues$ID = map[, 1]
    output = c(output, breedingvalues = list(breedingvalues), 
               summary.breedingvalues = list(summary.breedingvalues), 
               Ainversemapping = list(map))
  }
  if (class(Ainverse) != "ped") {
    output = c(output, breedingvalues = list(breedingvalues), 
               summary.breedingvalues = list(summary.breedingvalues))
  }
  if (dic == TRUE) {
    output = c(output, dic = list(dic2$dic))
  }
  if (type.data == "gaussian" & dic == TRUE | type.data == 
        "gaussian" & sigma.e == TRUE | type.data == "gaussian" & 
        linear.comb == TRUE | type.data == "gaussian" & dic == 
        TRUE & sigma.e == TRUE & linear.comb == TRUE | type.data == 
        "gaussian" & sigma.e == TRUE & dic == TRUE | type.data == 
        "gaussian" & sigma.e == TRUE & linear.comb == TRUE | 
        type.data == "gaussian" & dic == TRUE & linear.comb == 
        TRUE) {
    output = c(output, formula.gaussianDIC = formulaG0, inla.call.gaussianDIC = list(model.animalG0$call), 
               model.covariatesDIC = list(model.animalG0$model.random))
  }
  if (type.data == "binomial") {
    output = c(output, binary.h = list(hbin))
  }
  if (linear.comb == TRUE) {
    output = c(output, linear.combinations = list(linear.combinations))
  }
  if (!is.null(lambda) & type.data == "poisson" | !is.null(lambda) & 
        type.data == "zeroinflatedpoisson1" | !is.null(lambda) & 
        type.data == "zeroinflatedpoisson0") {
    output = c(output, poisson.h = list(hpois))
  }
  if (type.data == "gaussian" & standardize == TRUE) {
    output = c(output, Message = paste("Response variable is standardized"))
  }
  if (type.data == "gaussian" & sigma.e == TRUE) {
    output = c(output, sigma.e = list(sigmae), gaussian.h = list(h2))
  }
  class(output) = "Animalinla"
  return(output)
}
