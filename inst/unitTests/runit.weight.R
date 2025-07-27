test.mdw <- function() {

    library("RUnit")
    library("mdw")
    library("MASS")
    library("Rmosek")
    RNGkind("Mersenne-Twister", "Inversion")    
    tolerance=1e-1
    # R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
    if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu","x86_64, darwin15.6.0")) tolerance=1e-6 
    print(tolerance)

    ##############################
    # maximum entropy weights
    # generate cotinuous dataset X, bandwidth h 
    set.seed(1)
    X <- mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    # specified bandwidth
    h <- 1
    w.e <- entropy.weight(X,h)
    v.e <- asym.v.e(X,w.e,h)
    # confirm res
    checkEqualsNumeric(w.e, c(0.3518879, 0.3436133, 0.3044988), tolerance=tolerance)
    checkEqualsNumeric(v.e, matrix(c(0.04372809, -0.01454785, -0.01454785, 0.05235954), byrow = T, ncol = 2), tolerance=tolerance)
    # bandwidth selection by bw.nrd
    h <- get.bw(X,'nrd','na')
    w.e <- entropy.weight(X,h)
    v.e <- asym.v.e(X,w.e,h)
    # confirm res
    checkEqualsNumeric(h, 0.5153582, tolerance=tolerance)
    checkEqualsNumeric(w.e, c(0.3467989, 0.3367580, 0.3164431), tolerance=tolerance)
    checkEqualsNumeric(v.e, matrix(c(0.03285821, -0.01358718, -0.01358718, 0.04022874), byrow = T, ncol = 2), tolerance=tolerance)
    # generate categorical dataset dat
    set.seed(1)
    tmp=mvrnorm(n=10,mu=c(0,0,0),Sigma = diag(3))
    dat=t(apply(tmp, 1, function(x) cut(x,c(-Inf,-0.5,0.5,Inf),labels=1:3)))
    w.e=entropy.weight(dat,'na')
    # confirm res
    checkEqualsNumeric(w.e, c(0.4080125, 0.1840237, 0.4079638), tolerance=tolerance)
    ##############################
    # maximum variance weights
    # generate X
    set.seed(1)
    X <- mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    w.v.o <- var.weight(X)
    w.v.m=var.weight(X,'mosek')
    v.v <- asym.v.v(X,w.v)
    # confirm res
    checkEqualsNumeric(w.v.o, c(0.3593430, 0.3451024, 0.2955546), tolerance=tolerance)
    checkEqualsNumeric(w.v.m, c(0.3593872, 0.3451175, 0.2954953), tolerance=tolerance)
    checkEqualsNumeric(v.v, matrix(c(0.06705141, -0.01533853, -0.01533853, 0.10214535), byrow = T, ncol = 2), tolerance=tolerance)
        



}
