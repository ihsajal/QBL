
QBL.CC <- function (dat, numSNPs, maxMissingGenos = 1, allelic = TRUE,
                    haplo.baseline = "missing", cov.baseline = "missing", a = 20, b = 20,
                    start.beta = 0.01, gamma = 0.01, lambda = 1, D = 0, e = 0.1,
                    seed = NULL, burn.in = NULL, num.it = NULL, start.u = 0.01,
                    sigma_sq_u = 1, sigma_sq_b = NULL, sigma_sq_c = NULL,
                    e_allHap = 0.4, print.freq.ci = FALSE, print.lambda.ci = FALSE,
                    print.D.ci = FALSE, print.sigma_sq_u.ci = FALSE,
                    print.sigma_sq_b.ci = FALSE, print.sigma_sq_c.ci = FALSE)
{
  haplos.list <- pre.hapassoc(dat, numSNPs = numSNPs, pooling.tol = 0,
                              allelic = allelic, maxMissingGenos = maxMissingGenos)

  n.Y <- dim(haplos.list$nonHaploDM)[2] #number of Y
  temp = cbind(haplos.list$ID, haplos.list$wt, haplos.list$nonHaploDM,
               haplos.list$haploDM) #temporary data
  temp = temp[order(temp[, 1]), ] # Order wrt ID
  colnames(temp) = c("ID", "wt", "Yb", "Yc", colnames(temp[, 5:ncol(temp)]))
  haplos.list$ID = temp[, 1] # ID in haplotype list
  haplos.list$wt = temp[, 2] # weights in haplotype list
  haplos.list$nonHaploDM = temp[, 3:(3 + n.Y - 1)] # response Y in haplotype list
  haplos.list$haploDM = temp[,(3 + n.Y):ncol(temp)] # Haplotype data matrix in haplotype list
  haplo.names <- names(haplos.list$initFreq)
  Freq <- haplos.list$initFreq
  # Setting the maximum frequency as baseline
  if(haplo.baseline == "missing"){haplo.baseline <- haplo.names[which.max(Freq)]}
  column.subset <- colnames(haplos.list$haploDM) != haplo.baseline
  # Design matrix [X_ir]
  haplo.mat <- as.matrix(haplos.list$haploDM[, column.subset, drop = FALSE])
  hap.name <- colnames(haplo.mat)
  hdat <- cbind(haplos.list$nonHaploDM, haplos.list$haploDM[, column.subset])
  ID <- haplos.list$ID
  N <- length(unique(ID))
  yb = as.numeric(hdat[, 1])
  yb = (yb - mean(yb))/sd(yb)
  yc = as.numeric(hdat[, 2])
  yc = (yc - mean(yc))/sd(yc)
  # calculation of correlation
  y_cor = cor(yb, yc)

  # Haplotype frequencies
  freq.new <- Freq[c(colnames(haplo.mat), haplo.baseline)]

  num.haplo.id <- as.vector(table(ID))
  x.length <- as.integer(dim(haplo.mat)[2]) # number of columns in design matrix
  haplo.map <- matrix(rep(NA, 2 * (dim(haplo.mat)[1])), ncol = 2)
  for (i in 1:dim(haplo.mat)[1]) {
    for (j in 1:dim(haplo.mat)[2]) {
      if (haplo.mat[i, j] == 2) {
        haplo.map[i, 1] = j
        haplo.map[i, 2] = j
        break
      }
      if (haplo.mat[i, j] == 1) {
        if (is.na(haplo.map[i, 1]) == TRUE)
          haplo.map[i, 1] = j
        else {
          haplo.map[i, 2] = j
          break
        }
      }
    }
    if (is.na(haplo.map[i, 1]) == TRUE) {
      haplo.map[i, 1] = dim(haplo.mat)[2] + 1
      haplo.map[i, 2] = dim(haplo.mat)[2] + 1
    }
    if (is.na(haplo.map[i, 2]) == TRUE)
      haplo.map[i, 2] = dim(haplo.mat)[2] + 1
  }
  uniq.mat <- unique(haplo.mat)
  uniq.map <- unique(haplo.map)
  if (is.null(sigma_sq_b) == TRUE)
    sigma_sq_b = var(dat[,1])
  if (is.null(sigma_sq_c) == TRUE)
    sigma_sq_c = var(dat[, 2])
  if (is.null(seed) == FALSE)
    set.seed(seed)
  if (is.null(burn.in) == TRUE)
    burn.in = 50000
  if (is.null(num.it) == TRUE)
    num.it = 1e+05
  beta_b = rep(start.beta, x.length + 1)
  beta_c = rep(start.beta, x.length + 1)
  u = rep(start.u, N)
  sigma_u = sqrt(sigma_sq_u)
  beta_b.out <- numeric((num.it - burn.in) * (1 + x.length))
  beta_c.out <- numeric((num.it - burn.in) * (1 + x.length))
  lambda_b.out <- numeric(num.it - burn.in)
  lambda_c.out <- numeric(num.it - burn.in)
  freq.out <- numeric((num.it - burn.in) * (x.length + 1))
  D.out <- numeric(num.it - burn.in)
  sigma_u.out <- numeric(num.it - burn.in)
  sigma_sq_b.out <- numeric(num.it - burn.in)
  sigma_sq_c.out <- numeric(num.it - burn.in)
  out = .C("mcmc_BiCC", as.integer(haplo.map), as.integer(dim(haplo.mat)[1]),
           as.double(yb), as.double(yc), as.integer(N), as.integer(num.haplo.id),
           as.integer(x.length), as.double(freq.new), as.double(D),
           as.double(beta_b[1]), as.double(beta_b),
           as.double(beta_c[1]), as.double(beta_c),
           as.double(a), as.double(b), as.double(a), as.double(b), as.integer(t(uniq.map)),
           as.integer(dim(uniq.mat)[1]), as.double(lambda),
           as.double(lambda), as.double(u), as.double(sigma_u), as.double(sigma_sq_b),
           as.double(sigma_sq_c), as.integer(num.it), as.integer(burn.in),
           beta_b.out = as.double(beta_b.out), beta_c.out = as.double(beta_c.out),
           lambda_b.out = as.double(lambda_b.out), lambda_c.out = as.double(lambda_c.out),
           freq.out = as.double(freq.out), D.out = as.double(D.out),
           sigma_u.out = as.double(sigma_u.out),
           sigma_sq_b.out = as.double(sigma_sq_b.out),
           sigma_sq_c.out = as.double(sigma_sq_c.out), PACKAGE = "QBL")

  beta_b.out <- matrix(out$beta_b.out, nrow = num.it - burn.in, byrow = TRUE)
  beta_c.out <- matrix(out$beta_c.out, nrow = num.it - burn.in, byrow = TRUE)
  freq.out <- matrix(out$freq.out, nrow = num.it - burn.in, byrow = TRUE)
  D.out <- matrix(out$D.out, nrow = num.it - burn.in, byrow = TRUE)
  lambda_b.out <- matrix(out$lambda_b.out, nrow = num.it - burn.in, byrow = TRUE)
  lambda_c.out <- matrix(out$lambda_c.out, nrow = num.it - burn.in, byrow = TRUE)
  sigma_sq_b.out <- matrix(out$sigma_sq_b.out, nrow = num.it - burn.in, byrow = TRUE)
  sigma_sq_c.out <- matrix(out$sigma_sq_c.out, nrow = num.it - burn.in, byrow = TRUE)
  sigma_u.out <- matrix(out$sigma_u.out, nrow = num.it - burn.in, byrow = TRUE)

  hap.name <- colnames(hdat[, -c(1, 2)])
  n.hap <- length(hap.name)
  colnames(beta_b.out) <- c("beta0", hap.name)
  colnames(beta_c.out) <- c("beta0", hap.name)
  beta.out <- cbind(beta_b.out[, hap.name], beta_c.out[, hap.name])
  post.mean.beta_b <- numeric(1 + x.length)
  post.mean.beta_c <- numeric(1 + x.length)
  post.mean.freq <- numeric(x.length + 1)
  ci.beta_b.95 <- numeric((1 + x.length) * 2)
  ci.beta_c.95 <- numeric((1 + x.length) * 2)
  post.mean.beta_b <- round(post.mean.beta_b, 3)
  post.mean.beta_c <- round(post.mean.beta_c, 3)
  names(post.mean.beta_b) <- c("beta0", hap.name)
  names(post.mean.beta_c) <- c("beta0", hap.name)
  k <- 1
  for (i in 1:(1 + x.length)) {
    ci.beta_b.95[k] <- quantile(beta_b.out[, i], probs = 0.025)
    ci.beta_b.95[k + 1] <- quantile(beta_b.out[, i], probs = 0.975)
    ci.beta_c.95[k] <- quantile(beta_c.out[, i], probs = 0.025)
    ci.beta_c.95[k + 1] <- quantile(beta_c.out[, i], probs = 0.975)
    post.mean.beta_b[i] <- mean(beta_b.out[, i])
    post.mean.beta_c[i] <- mean(beta_c.out[, i])
    post.mean.freq[i] <- mean(freq.out[, i])
    k <- k + 2
  }
  ci.beta_b.95 <- matrix(ci.beta_b.95, nrow = 1 + x.length, ncol = 2, byrow = TRUE)
  ci.beta_c.95 <- matrix(ci.beta_c.95, nrow = 1 + x.length, ncol = 2, byrow = TRUE)
  post.mean.beta_b <- round(post.mean.beta_b, 3)
  post.mean.beta_c <- round(post.mean.beta_c, 3)
  ci.beta_b.95 <- round(ci.beta_b.95, 2)
  ci.beta_c.95 <- round(ci.beta_c.95, 2)
  names(post.mean.beta_b) <- c("beta0", hap.name)
  names(post.mean.beta_c) <- c("beta0", hap.name)
  ci.beta_b.95 <- data.frame(c("beta0", hap.name), ci.beta_b.95)
  colnames(ci.beta_b.95) <- c("Hap", "Lower", "Upper")
  ci.beta_c.95 <- data.frame(c("beta0", hap.name), ci.beta_c.95)
  colnames(ci.beta_c.95) <- c("Hap", "Lower", "Upper")
  names(post.mean.freq) <- c(hap.name, haplo.baseline)
  BF_bivariate_hap = rep(0, x.length)
  max_all_beta = apply(abs(beta.out), 1, max)
  prior_H0 = Sprob(x.length * 2, 0, a, b, e_allHap)
  prior_Ha = 1 - prior_H0
  post_H0 = length(max_all_beta[max_all_beta <= e_allHap])/length(max_all_beta)
  post_Ha = 1 - post_H0
  all_BF = (post_Ha/post_H0)/(prior_Ha/prior_H0)
  for (j in 1:x.length) {
    ebeta = beta.out[, c(j, j + x.length)]
    maxbeta = apply(abs(ebeta), 1, max)
    prior_H0 = Sprob(2, 0, a, b, e)
    prior_Ha = 1 - prior_H0
    post_H0 = length(maxbeta[maxbeta <= e])/length(maxbeta)
    post_Ha = 1 - post_H0
    BF_bivariate_hap[j] = round((post_Ha/post_H0)/(prior_Ha/prior_H0), 4)
  }
  names(BF_bivariate_hap) <- hap.name
  if (print.lambda.ci == TRUE) {
    ci.lambda_b <- numeric(2)
    ci.lambda_c <- numeric(2)
    ci.lambda_b <- c(quantile(out$lambda_b.out, probs = 0.025),
                     quantile(out$lambda_b.out, probs = 0.975))
    ci.lambda_c <- c(quantile(out$lambda_c.out, probs = 0.025),
                     quantile(out$lambda_c.out, probs = 0.975))
  }
  else {
    ci.lambda_b <- ""
    ci.lambda_c <- ""
  }
  if (print.D.ci == TRUE) {
    ci.D <- numeric(2)
    ci.D <- c(quantile(out$D.out, probs = 0.025), quantile(out$D.out, probs = 0.975))
  }
  else ci.D <- ""
  if (print.sigma_sq_u.ci == TRUE) {
    ci.sigma_u <- numeric(2)
    ci.sigma_u <- c(quantile(out$sigma_u.out, probs = 0.025),
                    quantile(out$sigma_u.out, probs = 0.975))
  }
  else ci.sigma_u <- ""
  if(print.sigma_sq_b.ci==TRUE){
    ci.sigma_sq_b <- numeric(2)
    ci.sigma_sq_b <- c(quantile(out$sigma_sq_b.out, probs = 0.025),
                       quantile(out$sigma_sq_b.out, probs = 0.975))
  }
  else ci.sigma_sq_b <- ""
  if (print.sigma_sq_c.ci == TRUE) {
    ci.sigma_sq_c <- numeric(2)
    ci.sigma_sq_c <- c(quantile(out$sigma_sq_c.out, probs = 0.025),
                       quantile(out$sigma_sq_c.out, probs = 0.975))
  }
  else ci.sigma_sq_c <- ""
  if (print.freq.ci == TRUE) {
    ci.freq <- numeric((x.length + 1) * 2)
    k <- 1
    for (i in 1:(1 + x.length)) {
      ci.freq[k] <- quantile(freq.out[, i], probs = 0.025)
      ci.freq[k + 1] <- quantile(freq.out[, i], probs = 0.975)
      k <- k + 2
    }
    ci.freq <- matrix(ci.freq, nrow = (x.length + 1), ncol = 2, byrow = TRUE)
    ci.freq <- data.frame(c(hap.name, haplo.baseline), ci.freq)
    colnames(ci.freq) <- c("Hap", "Lower", "Upper")
  }
  else {ci.freq <- ""}

  ans <- list(BF_bivariate_hap = BF_bivariate_hap, BF_bivariate_allHap = all_BF,
              beta_b = post.mean.beta_b, beta_c = post.mean.beta_c,
              CI.beta_b = ci.beta_b.95, CI.beta_c = ci.beta_c.95,
              freq = post.mean.freq, CI.freq = ci.freq, CI.lambda_b = ci.lambda_b,
              CI.lambda_c = ci.lambda_c, CI.D = ci.D, CI.sigma_u = ci.sigma_u,
              CI.sigma_sq_b = ci.sigma_sq_b, CI.sigma_sq_c = ci.sigma_sq_c)

  return(ans)
}
