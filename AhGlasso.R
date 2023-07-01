library(glasso)

CV_ahglasso <-
  function(x = NULL, # expression data, nrow: features; ncol: subjects.
           known_ppi = NULL, # previous known PPI
           lambda = NULL,
           nlambda = 15,
           lambda.min.ratio = 0.1, # ratio for minial lambda in a grid search
           scr = TRUE, # optional screening to speed up
           rho = NULL,
           weight = 0.1,
           eta = 0,
           verbose = FALSE,
           eps = 1e-08,
           seed = NULL, # random seed for K fold sample split
           K = 5, # K fold to optimize lambda
           crit.cv = c("PBIC", "BIC", "loglik", "ploglik", "AIC",  "EBIC"),
           trace = c("progress", "print", "none"),
           ...) {
    X = t(x)
    X = scale(X, center = TRUE, scale = TRUE)
    S = cor(X)
    if (!is.null(lambda))
      nlambda = length(lambda)
    if (is.null(lambda))
    {
      if (is.null(nlambda))
        nlambda = 15
      if (is.null(lambda.min.ratio))
        lambda.min.ratio = 0.1
      lambda.max = max(max(S - diag(p)), -min(S - diag(p)))
      lambda.min = lambda.min.ratio * lambda.max
      lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
      rm(lambda.max, lambda.min, lambda.min.ratio)
      gc()
    }
    # match values
    crit.cv = match.arg(crit.cv)
    trace = match.arg(trace)
    lambda = sort(lambda)
    CV_nLogLikelihood = CV_pNLogLikelihood = CV_Pbic.score = matrix(0, length(lambda), K)
    CV_ebic.score = CV_aic.score = CV_bic.score =matrix(0, length(lambda), K)
    CV_density = CV_nLogLikelihood
    # set progress bar
    if (trace == "progress") {
      progress = txtProgressBar(max = K * length(lambda), style = 3)
    }
    # No need to create folds if K = 1
    if (K == 1) {
      # set sample size
      n = nrow(X)
      X.train = X
      empcov.train = S
      k = 1
      # loop over all tuning parameters
      for (loop.lambda in 1:length(lambda)) {
        # compute the penalized likelihood precision matrix
        # estimator
        fit = cGGM_BM_scr(
          t(X.train),
          known_ppi,
          lambda[loop.lambda],
          scr = scr,
          rho = rho,
          weight = weight,
          eta = eta,
          verbose = FALSE
        )
        siginv = fit$invcov
        no.edge = (sum(abs(siginv) > eps)) / 2
        CV_density[loop.lambda, k] = spa(fit$Adj)
        # compute the observed negative validation loglikelihood
        CV_nLogLikelihood[loop.lambda, k] = matTr(empcov.train %*% siginv) - 
                                  determinant(siginv, logarithm = T)$modulus
        CV_pNLogLikelihood[loop.lambda, k] = CV_nLogLikelihood[loop.lambda, k] + fit$lambda *
          sum(abs((1 - diag(p)) * siginv))
        CV_bic.score[loop.lambda, k] = matTr(empcov.train %*% siginv) -
          determinant(siginv, logarithm = T)$modulus + log(n) * no.edge / n
        CV_Pbic.score[loop.lambda, k] = CV_pNLogLikelihood[loop.lambda, k] + log(n) * no.edge / n
        CV_ebic.score[loop.lambda, k] = CV_bic.score[loop.lambda, k] + 4 * log(p) * no.edge / n
        CV_aic.score[loop.lambda, k] = CV_nLogLikelihood[loop.lambda, k] + 2 * no.edge / n
        
        # update progress bar
        if (trace == "progress") {
          setTxtProgressBar(progress, loop.lambda + (k - 1) *
                              length(lambda))
          
          # if not quiet, then print progress lambda
        } else if (trace == "print") {
          cat("\nFinished lambda = ", paste(loop.lambda, sep = ""))
        }
      }
      
      # if not quiet, then print progress kfold
      if (trace == "print") {
        cat("\nFinished fold ", paste(k, sep = ""))
      }
    }
    else {
      # designate folds and shuffle -- ensures randomized folds
      n = nrow(X)
      if (!is.null(seed)) {
        set.seed(seed)
      }
      ind = sample(n)
      
      # parse data into folds and perform CV
      for (k in 1:K) {
        # training set
        leave.out = ind[(1 + floor((k - 1) * n / K)):floor(k * n / K)]
        X.train = X[-leave.out, , drop = FALSE]
        # # sample covariances
        X.train = scale(X.train, center = TRUE, scale = TRUE)
        empcov.train = cor(X.train)
        
        # loop over all tuning parameters
        for (loop.lambda in 1:length(lambda)) {
          # compute the penalized likelihood precision matrix
          # estimator
          fit = cGGM_BM_scr(
            t(X.train),
            known_ppi,
            lambda[loop.lambda],
            scr = scr,
            rho = rho,
            weight = weight,
            eta = eta,
            verbose = FALSE
          )
          siginv = fit$invcov
          
          no.edge = (sum(abs(siginv) > eps)) / 2
          CV_density[loop.lambda, k] = spa(fit$Adj)
          # compute the observed negative validation loglikelihood
          CV_nLogLikelihood[loop.lambda, k] = matTr(empcov.train %*% siginv) - 
                                    determinant(siginv, logarithm = T)$modulus
          CV_pNLogLikelihood[loop.lambda, k] = CV_nLogLikelihood[loop.lambda, k] + 
                              fit$lambda *csum(abs((1 - diag(p)) * siginv))
          CV_bic.score[loop.lambda, k] = matTr(empcov.train %*% siginv) - 
            determinant(siginv, logarithm = T)$modulus + log(n) * no.edge / n
          CV_Pbic.score[loop.lambda, k] = CV_pNLogLikelihood[loop.lambda, k] + 
                                          log(n) * no.edge / n
          CV_ebic.score[loop.lambda, k] = CV_bic.score[loop.lambda, k] + 4 *
                                          log(p) * no.edge / n
          CV_aic.score[loop.lambda, k] = CV_nLogLikelihood[loop.lambda, k] + 
                                          2 * no.edge / n
          # update progress bar
          if (trace == "progress") {
            setTxtProgressBar(progress, loop.lambda + (k - 1) * length(lambda))
            # if not quiet, then print progress lambda
          } else if (trace == "print") {
            cat("\nFinished lambda = ", paste(loop.lambda, sep = ""))
          }
        }
        # if not quiet, then print progress kfold
        if (trace == "print") {
          cat("\nFinished fold ", paste(k, sep = ""))
        }
      }
    }
    
    if (crit.cv == "loglik") {
      CV_errors = CV_nLogLikelihood
    }
    if (crit.cv == "ploglik") {
      CV_errors = CV_pNLogLikelihood
    }
    if (crit.cv == "PBIC") {
      CV_errors = CV_Pbic.score
    }
    if (crit.cv == "BIC") {
      CV_errors = CV_bic.score
    }
    if (crit.cv == "EBIC") {
      CV_errors = CV_ebic.score
    }
    if (crit.cv == "AIC") {
      CV_errors = CV_aic.score
    }
    
    allMetrics = cbind(
      lambda,
      apply(CV_density, 1, mean),
      apply(CV_nLogLikelihood, 1, mean),
      apply(CV_pNLogLikelihood, 1, mean),
      apply(CV_bic.score, 1, mean),
      apply(CV_Pbic.score, 1, mean),
      apply(CV_ebic.score, 1, mean),
      apply(CV_weightedebic.score, 1, mean),
      apply(CV_aic.score, 1, mean)
    )
    colnames(allMetrics) <-
      c(
        "lambda",
        "density",
        "nLogLikelihood",
        "pNLogLikelihood",
        "BIC",
        "PBIC",
        "EBIC",
        "wEBIC",
        "AIC"
      )
    
    # determine optimal tuning parameters
    AVG = apply(CV_errors, 1, mean)
    STD = apply(CV_errors, 1, sd)
    if (K == 1) {
      bestLambda = lambda[which.min(AVG)]
      minerror = min(AVG)
    } else {
      minerror = min(AVG)
      minerror_sd = STD[which.min(AVG)]
      bestLambda = lambda[which.min(abs(AVG - (minerror + minerror_sd)))]
    }
    
    results = list(
      lambda = lambda,
      bestLambda = bestLambda,
      min.error = minerror,
      CV_density = CV_density,
      avg.error = AVG,
      cv.error = CV_errors,
      CV_nLogLikelihood = CV_nLogLikelihood,
      CV_pNLogLikelihood = CV_pNLogLikelihood,
      CV_bic.score = CV_bic.score,
      CV_Pbic.score = CV_Pbic.score,
      CV_ebic.score = CV_ebic.score,
      CV_weightedebic.score = CV_weightedebic.score,
      CV_aic.score = CV_aic.score,
      allMetrics = allMetrics
    )
    class(results) = "CVahglasso"
    return(results)
    
  }


cGGM_BM_scr <-
  function (x, # expression data, nrow: features; ncol: subjects.
            known_ppi = NULL, # previously known PPI 
            scr = TRUE, # optional screening to speed up
            lambda, 
            rho = NULL,
            weight = 0.1,
            eta = 0,
            verbose = FALSE,
            eps = 1e-08) {
    p <- nrow(x)
    n <- ncol(x)
    Adj = matrix(0, p, p)
    Ip = diag(rep(1, p))
    
    if (!is.null(known_ppi)) {
      one = known_ppi
      one[one > 0] = 1
      rownames(one) = rownames(x)
      penaltyfactor = 1 - known_ppi
    }
    else {
      one = matrix(0, p, p)
      rownames(one) = rownames(x)
      penaltyfactor = matrix(1, p, p)
    }
    
    if (abs(lambda) < eps) {
      stop("The penalty parameter lambda needs to be greater than zero!")
    }
    if (n < 10) {
      warning("The sample size is too small! Network estimate may be unreliable!")
    }
    
    if (is.null(rho)) {
      rhoM = matrix(0.1 * sqrt(log(p) / n), p, p)
      # rhoM = lambda*(1-known_ppi)
    }
    else if (is.matrix(rho)) {
      if (length(rho) != p * p)
        stop("The input matrix for \"rho\" must be of size ",
             p, " by ", p)
      rhoM = rho
    }
    else {
      rhoM = matrix(rho, p, p)
    }
    X = t(x)
    X = scale(X, center = TRUE, scale = TRUE)
    
    if (!is.null(weight)) {
      if (weight < -1e-16) {
        stop("Negative weight parameter detected! Please double check!")
      }
    }
    
    if (scr) {
      Sconv = cov(X)
      scr_index = matrix(1, p, p)
      scr_index[abs(Sconv) <= lambda] = 0
      diag(scr_index) = 0
    }
    
    for (i in 1:p) {
      Y = matrix(X[, i], ncol = 1)
      Xmat = X[,-i]
      infoInd = one[i,-i]
      beta = matrix(0, p - 1, 1)
      
      if (sum(infoInd == 1) == 0) {
        if (verbose) {
          cat("Incomplete information: no known ppi! \n ")
        }
        if (scr) {
          if (verbose) {
            cat("Screening: filter uncorrelated edges! \n ")
          }
          scr_mat = scr_index[i,-i]
          #print (paste0("No. of selected features: ", sum(scr_mat)))
          if (sum(scr_mat) > 0) {
            Xmat_scr = Xmat[, (scr_mat == 1), drop = FALSE]
            beta[(scr_mat == 1),] = glmnet.soft(Xmat_scr, Y, lambda = lambda)
          }
        }
        else{
          beta = glmnet.soft(Xmat, Y, lambda = lambda)
        }
      }
      
      else if (sum(infoInd == 1) == (p - 1)) {
        if (verbose) {
          cat("Complete information: known ppi! \n ")
        }
        Pmat = penaltyfactor[i,-i]
        beta = glmnet.ppi(Xmat,
                          Y,
                          lambda = lambda * weight,
                          penalty.factor = Pmat)
      }
      else {
        if (verbose) {
          cat("Incomplete information: some known ppi! \n ")
        }
        Xmat1 = matrix(Xmat[, (infoInd == 1)],
                       ncol = sum(infoInd == 1))
        Pmat = penaltyfactor[i,-i]
        Pmat1 = Pmat[(infoInd == 1)]
        Xmat2 = matrix(Xmat[, (infoInd == 0)],
                       ncol = sum(infoInd == 0))
        
        beta[(infoInd == 1),] = glmnet.ppi(Xmat1,
                                           Y,
                                           lambda = lambda * weight,
                                           penalty.factor = Pmat1)
        tmp = as.matrix(beta[(infoInd == 1),])
        res.glm = Y - Xmat1 %*% tmp
        if (verbose) {
          cat("Step 1 done, move on the regression on residual! \n ")
        }
        
        if (scr) {
          if (verbose) {
            cat("Screening: filter uncorrelated edges on step 2! \n ")
          }
          scr_mat = scr_index[i,-i]
          scr_mat_infoInd = (scr_mat == 1 & infoInd == 0)
          
          # print (paste0("No. of selected features: ", sum(scr_mat)))
          if (sum(scr_mat_infoInd) > 0) {
            Xmat_scr = matrix(Xmat[, scr_mat_infoInd],
                              ncol = sum(scr_mat_infoInd))
            
            #Xmat[, (scr_mat==1), drop = FALSE]
            
            beta[scr_mat_infoInd,] = glmnet.soft(Xmat_scr, res.glm, 
                                                 lambda = lambda)
          }
        }
        else{
          # print(dim(Xmat2))
          beta[(infoInd == 0),] = glmnet.soft(Xmat2,
                                              res.glm, lambda = lambda)
        }
      }
      Adj[i,-i] = as.vector(beta)
    }
    oriAdj = Adj
    Adj = (Adj + t(Adj)) / 2
    Adj = (abs(Adj) > eps)
    diag(Adj) = 0
    empcov = cov(X)
    if (kappa(empcov) > 1000) {
      empcov = empcov + eta * diag(p)
    }
    BIG = 1e+10
    rhoM[which(Adj == 0)] = BIG
    diag(rhoM) = 0.01 * sqrt(log(p) / n)
    obj <- glassoFast(empcov, rho = rhoM)
    siginv <- obj$wi
    partialCor <- Ip - cov2cor(siginv)
    rownames(partialCor) <- rownames(x)
    colnames(partialCor) <- rownames(x)
    return(
      list(
        Adj = abs(partialCor),
        invcov = siginv,
        lambda = lambda,
        one = one,
        weight = weight,
        oriAdj = oriAdj,
        empcov = empcov,
        rhoM = rhoM,
        Adj_index = Adj
      )
    )
  }
