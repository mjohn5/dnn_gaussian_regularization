

library(ncpen)

tb <- true.beta <- c(seq(0.01, 0.99, 0.01), seq(1,10))

n <- 1000; p <- length(tb)
err.sigma <- 0.1 #0.1, 1




for(ov.iter in 1:1000){ print(ov.iter)

  ####################
  ### Data Generation
  ####################

  X.mat <- matrix(rnorm(n*p), nrow = n, ncol = p)
  tb.vec <- matrix(tb, nrow = p, ncol = 1)

  y <- X.mat%*%tb.vec + rnorm(n, 0, err.sigma)

  #########################################
  ### Estimation with various penalties using the CCP-MLQA algorithm in 'ncpen' R package. Optimal \lambda is determined via CV and beta's corresponding to this CV-optimal lambda are extracted
  #########################################

  scad.fit <- ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "scad", intercept = "FALSE")
  scad.cv.fit <- cv.ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "scad", intercept = "FALSE")

  mcp.cv.fit <- cv.ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "mcp", intercept = "FALSE")

  lasso.cv.fit <- cv.ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "lasso", intercept = "FALSE")

  tlp.cv.fit <- cv.ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "tlp", intercept = "FALSE")

  classo.cv.fit <- cv.ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "classo", intercept = "FALSE")

  ridge.cv.fit <- cv.ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "ridge", intercept = "FALSE")

  sridge.cv.fit <- cv.ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "sridge", intercept = "FALSE")

  mbridge.cv.fit <- cv.ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "mbridge", intercept = "FALSE")

  mlog.cv.fit <- cv.ncpen(y.vec = y, x.mat = X.mat, family = "gaussian", penalty = "mlog", intercept = "FALSE")


  ############################################################
  #### Gradient Descent based estimation with Gaussian penalty
  ############################################################

  b.init <- array(0,p)

  tol <- 1e-12; max.iter <- 100000; lr <- c(rep(0.1, max.iter/5), rep(0.05, max.iter/5), rep(0.025, max.iter/5), rep(0.01, max.iter/5), rep(0.001, max.iter/5))

  lambda <- scad.fit$lambda
  kappa <- 1

  for(L in 1:length(lambda)){ #print(L)
   for(m in 1:max.iter){ #print(m)
        if(m == 1){ b.curr <- matrix(b.init, nrow = p, ncol = 1) }
        if(m > 1){  b.curr <- b.new }
                  g.vec <- (2/n)*t(X.mat)%*%((X.mat%*%b.curr) - y) + lambda[L]*2*kappa*b.curr*exp(-kappa*(b.curr^2))
          	       #print(g.vec[1:10])
	          b.new <- b.curr - lr[m]*g.vec
	               }
    if(L == 1){ b.mat <- b.new }
    if(L > 1){ b.mat <- cbind(b.mat, b.new) } }

################################################
#### CV for Gaussian penalty estimation to determine the optimal lambda (lambda.opt)
###########################################################

  prederr <- array(, length(lambda))

  for(L in 1:length(lambda)){ cv.err <- array(,5)
     for(cv.iter in 1:5){ X.train <- X.mat[-((((cv.iter-1)*200) + 1):(cv.iter*200)),]
                          X.test  <- X.mat[(((cv.iter-1)*200) + 1):(cv.iter*200),]
                          y.train <- y[-((((cv.iter-1)*200) + 1):(cv.iter*200))]
                          y.test  <- y[(((cv.iter-1)*200) + 1):(cv.iter*200)] 
      for(m in 1:max.iter){ #print(m)
          if(m == 1){ b.curr <- matrix(b.init, nrow = p, ncol = 1) }
          if(m > 1){  b.curr <- b.new }
                    g.vec <- (2/n)*t(X.train)%*%((X.train%*%b.curr) - y.train) + lambda[L]*2*kappa*b.curr*exp(-kappa*(b.curr^2))          	   
	            b.new <- b.curr - lr[m]*g.vec
	                   }
                    b.train.fit <- b.new
		    y.test.fitted <- X.test%*%b.train.fit
		    cv.err[cv.iter] <- sum((y.test - y.test.fitted)^2) }
		    prederr[L] <- mean(cv.err) }

  prederr[is.nan(prederr)] <- max(prederr[!is.nan(prederr)])

  lambda.opt     <- lambda[prederr == min(prederr)]
  lambda.indx    <- c(1:length(lambda))
  lambda.b.opt.L <- lambda.indx[lambda == lambda.opt]



##########################################################
#### Gradient Descent based estimation with Ridge penalty
##########################################################


  rb.init <- array(0,p)
  
  tol <- 1e-12; max.iter <- 100000; lr <- c(rep(0.1, max.iter/5), rep(0.05, max.iter/5), rep(0.025, max.iter/5), rep(0.01, max.iter/5), rep(0.001, max.iter/5))

  lambda <- scad.fit$lambda


  for(L in 1:length(lambda)){ #print(L)
    for(m in 1:max.iter){ #print(m)
        if(m == 1){ rb.curr <- matrix(rb.init, nrow = p, ncol = 1) }
        if(m > 1){  rb.curr <- rb.new }
                  g.vec <- (2/n)*t(X.mat)%*%((X.mat%*%rb.curr) - y) + lambda[L]*2*rb.curr
          	       #print(g.vec[1:10])
	          rb.new <- rb.curr - lr[m]*g.vec
	               }
    if(L == 1){ rb.mat <- rb.new }
    if(L > 1){ rb.mat <- cbind(rb.mat, rb.new) } }


################################################
#### CV for Ridge penalty estimation to determine the optimal lambda (lambda.opt)
###########################################################


  prederr.rb <- array(, length(lambda))

  for(L in 1:length(lambda)){ cv.err.rb <- array(,5)
     for(cv.iter in 1:5){ X.train <- X.mat[-((((cv.iter-1)*200) + 1):(cv.iter*200)),]
                          X.test  <- X.mat[(((cv.iter-1)*200) + 1):(cv.iter*200),]
                          y.train <- y[-((((cv.iter-1)*200) + 1):(cv.iter*200))]
                          y.test  <- y[(((cv.iter-1)*200) + 1):(cv.iter*200)] 
      for(m in 1:max.iter){ #print(m)
          if(m == 1){ rb.curr <- matrix(rb.init, nrow = p, ncol = 1) }
          if(m > 1){  rb.curr <- rb.new }
                    g.vec <- (2/n)*t(X.train)%*%((X.train%*%rb.curr) - y.train) + lambda[L]*2*kappa*rb.curr          	   
	            rb.new <- rb.curr - lr[m]*g.vec
	                   }
                    rb.train.fit <- rb.new
		    y.rb.test.fitted <- X.test%*%rb.train.fit
		    cv.err.rb[cv.iter] <- sum((y.test - y.rb.test.fitted)^2) }
		    prederr.rb[L] <- mean(cv.err.rb) }

  prederr.rb[is.nan(prederr.rb)] <- max(prederr.rb[!is.nan(prederr.rb)])
		  
  lambda.opt.rb <- lambda[(prederr.rb == min(prederr.rb))]
  lambda.indx <- c(1:length(lambda))
  lambda.rb.opt.L <- lambda.indx[lambda == lambda.opt.rb]



   ############################
   #### Results ###############
   ############################

  bias.ratio.gauss  <- sum(abs((b.mat[,lambda.b.opt.L] - tb)/tb))
  bias.ratio.ridge  <- sum(abs((rb.mat[,lambda.rb.opt.L] - tb)/tb))
  bias.ratio.scad   <- sum(abs((coef(scad.cv.fit)$beta - tb)/tb))
  bias.ratio.mcp    <- sum(abs((coef(mcp.cv.fit)$beta - tb)/tb))
  bias.ratio.lasso  <- sum(abs((coef(lasso.cv.fit)$beta - tb)/tb))
  bias.ratio.tlp    <- sum(abs((coef(tlp.cv.fit)$beta - tb)/tb))
  bias.ratio.classo <- sum(abs((coef(classo.cv.fit)$beta - tb)/tb))
  #sum(abs((coef(ridge.cv.fit)$beta - tb)/tb))
  bias.ratio.sridge  <- sum(abs((coef(sridge.cv.fit)$beta - tb)/tb))
  bias.ratio.mbridge <- sum(abs((coef(mbridge.cv.fit)$beta - tb)/tb))
  bias.ratio.mlog    <- sum(abs((coef(mlog.cv.fit)$beta - tb)/tb))

  b.fit <- matrix(b.mat[,lambda.b.opt.L], nrow = p, ncol = 1);    y.fit <- X.mat%*%b.fit; mse.gauss   <- (1/n)*sum((y - y.fit)^2)
  b.fit <- matrix(rb.mat[,lambda.rb.opt.L], nrow = p, ncol = 1);  y.fit <- X.mat%*%b.fit; mse.ridge   <- (1/n)*sum((y - y.fit)^2)
  b.fit <- matrix(coef(scad.cv.fit)$beta, nrow = p, ncol = 1);    y.fit <- X.mat%*%b.fit; mse.scad    <- (1/n)*sum((y - y.fit)^2)
  b.fit <- matrix(coef(mcp.cv.fit)$beta, nrow = p, ncol = 1);     y.fit <- X.mat%*%b.fit; mse.mcp     <- (1/n)*sum((y - y.fit)^2)
  b.fit <- matrix(coef(lasso.cv.fit)$beta, nrow = p, ncol = 1);   y.fit <- X.mat%*%b.fit; mse.lasso   <- (1/n)*sum((y - y.fit)^2)
  b.fit <- matrix(coef(tlp.cv.fit)$beta, nrow = p, ncol = 1);     y.fit <- X.mat%*%b.fit; mse.tlp     <- (1/n)*sum((y - y.fit)^2)
  b.fit <- matrix(coef(classo.cv.fit)$beta, nrow = p, ncol = 1);  y.fit <- X.mat%*%b.fit; mse.classo  <- (1/n)*sum((y - y.fit)^2)
  b.fit <- matrix(coef(sridge.cv.fit)$beta, nrow = p, ncol = 1);  y.fit <- X.mat%*%b.fit; mse.sridge  <- (1/n)*sum((y - y.fit)^2)
  b.fit <- matrix(coef(mbridge.cv.fit)$beta, nrow = p, ncol = 1); y.fit <- X.mat%*%b.fit; mse.mbridge <- (1/n)*sum((y - y.fit)^2)
  b.fit <- matrix(coef(mlog.cv.fit)$beta, nrow = p, ncol = 1);    y.fit <- X.mat%*%b.fit; mse.mlog     <- (1/n)*sum((y - y.fit)^2)

  bias.gauss   <- sum(abs((b.mat[,lambda.b.opt.L] - tb)))
  bias.ridge   <- sum(abs((rb.mat[,lambda.rb.opt.L] - tb)))
  bias.scad    <- sum(abs((coef(scad.cv.fit)$beta - tb)))
  bias.mcp     <- sum(abs((coef(mcp.cv.fit)$beta - tb)))
  bias.lasso   <- sum(abs((coef(lasso.cv.fit)$beta - tb)))
  bias.tlp     <- sum(abs((coef(tlp.cv.fit)$beta - tb)))
  bias.classo  <- sum(abs((coef(classo.cv.fit)$beta - tb)))
  #sum(abs((coef(ridge.cv.fit)$beta - tb)))
  bias.sridge  <- sum(abs((coef(sridge.cv.fit)$beta - tb)))
  bias.mbridge <- sum(abs((coef(mbridge.cv.fit)$beta - tb)))
  bias.mlog    <- sum(abs((coef(mlog.cv.fit)$beta - tb)))

  bias.gauss.101.109   <- sum(abs((b.mat[101:109, lambda.b.opt.L] - tb[101:109])))
  bias.ridge.101.109   <- sum(abs((rb.mat[101:109,lambda.rb.opt.L] - tb[101:109])))
  bias.scad.101.109    <- sum(abs((coef(scad.cv.fit)$beta[101:109]  - tb[101:109])))
  bias.mcp.101.109     <- sum(abs((coef(mcp.cv.fit)$beta[101:109]   - tb[101:109])))
  bias.lasso.101.109   <- sum(abs((coef(lasso.cv.fit)$beta[101:109] - tb[101:109])))
  bias.tlp.101.109     <- sum(abs((coef(tlp.cv.fit)$beta[101:109] - tb[101:109])))
  bias.classo.101.109  <- sum(abs((coef(classo.cv.fit)$beta[101:109] - tb[101:109])))
  bias.sridge.101.109  <- sum(abs((coef(sridge.cv.fit)$beta[101:109] - tb[101:109])))
  bias.mbridge.101.109 <- sum(abs((coef(mbridge.cv.fit)$beta[101:109] - tb[101:109])))
  bias.mlog.101.109    <- sum(abs((coef(mlog.cv.fit)$beta[101:109] - tb[101:109])))

  bias.gauss.75.100   <- sum(abs((b.mat[75:100, lambda.b.opt.L] - tb[75:100])))
  bias.ridge.75.100   <- sum(abs((rb.mat[75:100,lambda.rb.opt.L] - tb[75:100])))
  bias.scad.75.100    <- sum(abs((coef(scad.cv.fit)$beta[75:100]  - tb[75:100])))
  bias.mcp.75.100     <- sum(abs((coef(mcp.cv.fit)$beta[75:100]   - tb[75:100])))
  bias.lasso.75.100   <- sum(abs((coef(lasso.cv.fit)$beta[75:100] - tb[75:100])))
  bias.tlp.75.100     <- sum(abs((coef(tlp.cv.fit)$beta[75:100] - tb[75:100])))
  bias.classo.75.100  <- sum(abs((coef(classo.cv.fit)$beta[75:100] - tb[75:100])))
  bias.sridge.75.100  <- sum(abs((coef(sridge.cv.fit)$beta[75:100] - tb[75:100])))
  bias.mbridge.75.100 <- sum(abs((coef(mbridge.cv.fit)$beta[75:100] - tb[75:100])))
  bias.mlog.75.100    <- sum(abs((coef(mlog.cv.fit)$beta[75:100] - tb[75:100])))

  bias.gauss.50.75   <- sum(abs((b.mat[50:75,lambda.b.opt.L] - tb[50:75])))
  bias.ridge.50.75   <- sum(abs((rb.mat[50:75,lambda.rb.opt.L] - tb[50:75])))
  bias.scad.50.75    <- sum(abs((coef(scad.cv.fit)$beta[50:75]  - tb[50:75])))
  bias.mcp.50.75     <- sum(abs((coef(mcp.cv.fit)$beta[50:75]   - tb[50:75])))
  bias.lasso.50.75   <- sum(abs((coef(lasso.cv.fit)$beta[50:75] - tb[50:75])))
  bias.tlp.50.75     <- sum(abs((coef(tlp.cv.fit)$beta[50:75] - tb[50:75])))
  bias.classo.50.75  <- sum(abs((coef(classo.cv.fit)$beta[50:75] - tb[50:75])))
  bias.sridge.50.75  <- sum(abs((coef(sridge.cv.fit)$beta[50:75] - tb[50:75])))
  bias.mbridge.50.75 <- sum(abs((coef(mbridge.cv.fit)$beta[50:75] - tb[50:75])))
  bias.mlog.50.75    <- sum(abs((coef(mlog.cv.fit)$beta[50:75] - tb[50:75])))


  bias.gauss.25.50   <- sum(abs((b.mat[25:50, lambda.b.opt.L] - tb[25:50])))
  bias.ridge.25.50   <- sum(abs((rb.mat[25:50, lambda.rb.opt.L] - tb[25:50])))
  bias.scad.25.50    <- sum(abs((coef(scad.cv.fit)$beta[25:50]  - tb[25:50])))
  bias.mcp.25.50     <- sum(abs((coef(mcp.cv.fit)$beta[25:50]   - tb[25:50])))
  bias.lasso.25.50   <- sum(abs((coef(lasso.cv.fit)$beta[25:50] - tb[25:50])))
  bias.tlp.25.50     <- sum(abs((coef(tlp.cv.fit)$beta[25:50] - tb[25:50])))
  bias.classo.25.50  <- sum(abs((coef(classo.cv.fit)$beta[25:50] - tb[25:50])))
  bias.sridge.25.50  <- sum(abs((coef(sridge.cv.fit)$beta[25:50] - tb[25:50])))
  bias.mbridge.25.50 <- sum(abs((coef(mbridge.cv.fit)$beta[25:50] - tb[25:50])))
  bias.mlog.25.50    <- sum(abs((coef(mlog.cv.fit)$beta[25:50] - tb[25:50])))


  bias.gauss.1.25   <- sum(abs((b.mat[1:25, lambda.b.opt.L] - tb[1:25])))
  bias.ridge.1.25   <- sum(abs((rb.mat[1:25, lambda.rb.opt.L] - tb[1:25])))
  bias.scad.1.25    <- sum(abs((coef(scad.cv.fit)$beta[1:25]  - tb[1:25])))
  bias.mcp.1.25     <- sum(abs((coef(mcp.cv.fit)$beta[1:25]   - tb[1:25])))
  bias.lasso.1.25   <- sum(abs((coef(lasso.cv.fit)$beta[1:25] - tb[1:25])))
  bias.tlp.1.25     <- sum(abs((coef(tlp.cv.fit)$beta[1:25] - tb[1:25])))
  bias.classo.1.25  <- sum(abs((coef(classo.cv.fit)$beta[1:25] - tb[1:25])))
  bias.sridge.1.25  <- sum(abs((coef(sridge.cv.fit)$beta[1:25] - tb[1:25])))
  bias.mbridge.1.25 <- sum(abs((coef(mbridge.cv.fit)$beta[1:25] - tb[1:25])))
  bias.mlog.1.25    <- sum(abs((coef(mlog.cv.fit)$beta[1:25] - tb[1:25])))


  #coef(tlp.cv.fit)
  #coef(classo.cv.fit)
  #coef(ridge.cv.fit)
  #coef(sridge.cv.fit)
  #coef(mbridge.cv.fit)
  #coef(mlog.cv.fit)



rslts <- c(mse.gauss, mse.ridge, mse.scad,  mse.mcp,   mse.lasso,  mse.tlp,   mse.classo, mse.sridge, mse.mbridge, mse.mlog, bias.ratio.gauss, bias.ratio.ridge, bias.ratio.scad,  bias.ratio.mcp,  bias.ratio.lasso, bias.ratio.tlp,  bias.ratio.classo, bias.ratio.sridge,  bias.ratio.mbridge, bias.ratio.mlog, bias.gauss, bias.ridge, bias.scad,  bias.mcp,   bias.lasso,  bias.tlp,   bias.classo, bias.sridge, bias.mbridge, bias.mlog, bias.gauss.101.109,  bias.ridge.101.109,  bias.scad.101.109,   bias.mcp.101.109,   bias.lasso.101.109,  bias.tlp.101.109,    bias.classo.101.109, bias.sridge.101.109,  bias.mbridge.101.109, bias.mlog.101.109, bias.gauss.75.100, bias.ridge.75.100, bias.scad.75.100,  bias.mcp.75.100,   bias.lasso.75.100,  bias.tlp.75.100,   bias.classo.75.100, bias.sridge.75.100, bias.mbridge.75.100, bias.mlog.75.100, bias.gauss.50.75, bias.ridge.50.75, bias.scad.50.75,  bias.mcp.50.75,   bias.lasso.50.75,  bias.tlp.50.75,   bias.classo.50.75, bias.sridge.50.75, bias.mbridge.50.75, bias.mlog.50.75, bias.gauss.25.50, bias.ridge.25.50, bias.scad.25.50,  bias.mcp.25.50,   bias.lasso.25.50,  bias.tlp.25.50,   bias.classo.25.50, bias.sridge.25.50, bias.mbridge.25.50, bias.mlog.25.50, bias.gauss.1.25, bias.ridge.1.25, bias.scad.1.25,  bias.mcp.1.25,   bias.lasso.1.25,  bias.tlp.1.25,   bias.classo.1.25, bias.sridge.1.25, bias.mbridge.1.25, bias.mlog.1.25)    

c.names <- c("mse.gauss", "mse.ridge", "mse.scad",  "mse.mcp",   "mse.lasso",  "mse.tlp",   "mse.classo", "mse.sridge", "mse.mbridge", "mse.mlog", "bias.ratio.gauss", "bias.ratio.ridge", "bias.ratio.scad",  "bias.ratio.mcp",  "bias.ratio.lasso", "bias.ratio.tlp",  "bias.ratio.classo", "bias.ratio.sridge",  "bias.ratio.mbridge", "bias.ratio.mlog", "bias.gauss", "bias.ridge", "bias.scad",  "bias.mcp",   "bias.lasso",  "bias.tlp",   "bias.classo", "bias.sridge", "bias.mbridge", "bias.mlog", "bias.gauss.101.109",  "bias.ridge.101.109",  "bias.scad.101.109",   "bias.mcp.101.109",   "bias.lasso.101.109",  "bias.tlp.101.109",    "bias.classo.101.109", "bias.sridge.101.109",  "bias.mbridge.101.109", "bias.mlog.101.109", "bias.gauss.75.100", "bias.ridge.75.100", "bias.scad.75.100",  "bias.mcp.75.100",   "bias.lasso.75.100",  "bias.tlp.75.100",   "bias.classo.75.100", "bias.sridge.75.100", "bias.mbridge.75.100", "bias.mlog.75.100", "bias.gauss.50.75", "bias.ridge.50.75", "bias.scad.50.75",  "bias.mcp.50.75",   "bias.lasso.50.75",  "bias.tlp.50.75",   "bias.classo.50.75", "bias.sridge.50.75", "bias.mbridge.50.75", "bias.mlog.50.75", "bias.gauss.25.50", "bias.ridge.25.50", "bias.scad.25.50",  "bias.mcp.25.50",   "bias.lasso.25.50",  "bias.tlp.25.50",   "bias.classo.25.50", "bias.sridge.25.50", "bias.mbridge.25.50", "bias.mlog.25.50", "bias.gauss.1.25", "bias.ridge.1.25", "bias.scad.1.25",  "bias.mcp.1.25",   "bias.lasso.1.25",  "bias.tlp.1.25",   "bias.classo.1.25", "bias.sridge.1.25", "bias.mbridge.1.25", "bias.mlog.1.25")    


    if(ov.iter == 1){ rslts.mat <- rslts }
    if(ov.iter > 1){  rslts.mat <- rbind(rslts.mat, rslts) 

                      rslts.df1 <- rslts.mat
		      colnames(rslts.df1) <- c.names

                      rslts.df <- as.data.frame(rslts.df1)


                     #write.table(rslts.df, "...pathname.../sim.109.csv", sep = ",", row.names = FALSE)
                   }

           }


		  
