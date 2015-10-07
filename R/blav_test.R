blav_model_test <- function(lavmodel       = NULL, 
                            lavpartable    = NULL, 
                            lavsamplestats = NULL, 
                            lavoptions     = NULL, 
                            x              = NULL, 
                            VCOV           = NULL, 
                            lavcache       = NULL,
                            lavdata        = NULL,
                            lavjags        = NULL,
                            control        = list()) {


    TEST <- list()

    ## marginal log-likelihood approximation
    ## needs original partable with rhos
    mll <- margloglik(lavpartable, lavmodel, lavoptions, 
                      lavsamplestats, lavdata, lavcache, lavjags)
    ## TODO this will be TEST[[2]] once fitMeasures gives us
    ## DIC; DIC restricted to fitMeasures?
    ppp <- postpred(lavpartable, lavmodel, lavoptions,
                    lavsamplestats, lavdata, lavcache, lavjags)$ppval

    TEST[[1]] <- list(test="mloglik",
                      stat=as.numeric(mll),
                      stat.group=as.numeric(NA),
                      df=as.integer(NA),
                      refdistr="NA",
                      pvalue=as.numeric(NA))

    TEST[[2]] <- list(test="ppp",
                      ## DIC: 2*ll(theta_hat) - 4*mean(ll(theta_samp))
                      stat=as.numeric(ppp),
                      stat.group=as.numeric(NA),
                      df=as.integer(NA),
                      refdistr="NA",
                      pvalue=as.numeric(NA))

    TEST
}
