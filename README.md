JM: Joint Models for Longitudinal and Survival Data using Maximum Likelihood
================
[![Travis-CI Build Status](https://travis-ci.org/drizopoulos/JM.svg?branch=master)](https://travis-ci.org/drizopoulos/JM) [![CRAN status](http://www.r-pkg.org/badges/version/JM)](https://cran.r-project.org/package=JM) [![](https://cranlogs.r-pkg.org/badges/grand-total/JM)](https://CRAN.R-project.org/package=JM) [![Download counter](http://cranlogs.r-pkg.org/badges/JM)](https://cran.r-project.org/package=JM)
[![Research software impact](http://depsy.org/api/package/cran/JM/badge.svg)](http://depsy.org/package/r/JM)

Description
------------

This repository contains the source files for the R package <strong>JM</strong>. 
This package fits joint models for longitudinal and time-to-event data using maximum 
likelihood. These models are applicable in mainly two settings. First, when focus
is on the survival outcome and we wish to account for the effect of an endogenous 
(aka internal) time-dependent covariates measured with error. Second, when focus is on the
longitudinal outcome and we wish to correct for nonrandom dropout.

The basic joint-model-fitting function of the package is `jointModel()`. This accepts as
main arguments a linear mixed model fitted by function `lme()` from the 
[**nlme**](https://CRAN.R-project.org/package=nlme) package and a Cox model fitted using
function `coxph()` from the [**survival**](https://CRAN.R-project.org/package=survival) 
package.

Basic Features
------------

- It can fit joint models for a single continuous longitudinal outcome and a time-to-event
outcome. 

- For the survival outcome a relative risk models is assumed. The `method` argument of 
`jointModel()` can be used to define the type of baseline hazard function. Options are a 
B-spline approximation, a piecewise-constant function, the Weibull hazard and a completely
unspecified function (i.e., a discrete function with point masses at the unique event 
times).

- The user has now the option to define custom transformation functions for the terms of 
the longitudinal submodel that enter into the linear predictor of the survival submodel 
(arguments `derivForm`, `parameterization`). For example, the current value of the 
longitudinal outcomes, the velocity of the longitudinal outcome (slope), the area under
the longitudinal profile. From the aforementioned options, in each model up to two terms 
can be included. In addition, using argument `InterFact` interactions terms can be 
considered.

Dynamic predictions
------------

* Function `survfitJM()` computes dynamic survival probabilities.

* Function `predict()` computes dynamic predictions for the longitudinal outcome.

* Function `aucJM()` calculates time-dependent AUCs for joint models, and function 
`rocJM()` calculates the corresponding time-dependent sensitivities and specifies.

* Function `prederrJM()` calculates prediction errors for joint models.
