## predictions from blavaan object; similar to lavPredict, but lavPredict is never called
## overload standard R function `predict'
setMethod("predict", "blavaan",
function(object, newdata = NULL) {
    blavPredict(object, newdata = newdata)
})

blavPredict <- function(object, newdata = NULL, type = "lv", level = 1L) {

  stopifnot(inherits(object, "blavaan"))
  blavmodel <- object@Model
  blavpartable <- object@ParTable
  blavsamplestats <- object@SampleStats
  blavdata <- object@Data
  standata <- object@external$mcmcdata
  
  type <- tolower(type)
  if(type %in% c("latent", "lv", "factor", "factor.score", "factorscore"))
      type <- "lv"
  if(type %in% c("ov","yhat"))
      type <- "yhat"
  if(type %in% c("ypred", "ydist"))
      type <- "ypred"
  if(type %in% c("ymis", "ovmis")){
      type <- "ymis"
      if(all(!is.na(unlist(blavdata@X)))) stop("blavaan ERROR: No missing data are present.", call. = FALSE)
  }

  lavopt <- lavInspect(object, "options")
  stantarget <- lavopt$target == "stan"

  if(lavInspect(object, "categorical") & type == "ymis") stop("blavaan ERROR: ymis is not yet implemented for ordinal models.", call. = FALSE)

  if(level == 2L){
    if(all(unlist(lavInspect(object, "nclusters")) == 1)) stop("blavaan ERROR: level 2 was requested but this does not appear to be a 2-level model.", call. = FALSE)
    if(type %in% c("yhat", "ypred", "ymis")) stop("blavaan ERROR: option", type, "is not yet implemented for two-level models.", call. = FALSE)
  }
  
  if(!is.null(newdata)) stop("blavaan ERROR: posterior predictions for newdata are not currently supported")
  
  ## lv: posterior dist of lvs (use blavInspect functionality); matrix frame
  ## lvmeans: use blavInspect functionality; matrix
  ## yhat: posterior expected value of ovs conditioned on lv samples; mcmc list
  ## ypred: posterior predictive distribution of ovs conditioned on lv samples; mcmc list
  ## ymis: posterior predictive distribution of missing values conditioned on observed values; matrix
  if(type == "lv") {
    FS <- do.call("rbind", blavInspect(object, 'lvs', level = level))

    ## N and latent variable names, to set dimensions
    lvmn <- lavInspect(object, "mean.lv")
    if(!inherits(lvmn, "list")){
      lvmn <- list(lvmn)
    }
    if(level == 1L){
      nlv <- length(lvmn[[1]])
      N <- sum(lavInspect(object, "ntotal"))
      etas <- names(lvmn[[1]])
    } else {
      nlv <- length(lvmn[[2]])
      N <- sum(unlist(lavInspect(object, "nclusters")))
      etas <- names(lvmn[[2]])
    }

    out <- lapply(1:NROW(FS), function(i) {
      rowmat <- matrix(FS[i,], N, nlv)
      colnames(rowmat) <- etas
      rowmat } )
  } else if(type == "lvmeans") {
    out <- blavInspect(object, 'lvmeans')
  } else if(type %in% c("yhat", "ypred", "ymis")) {
    if(!stantarget) stop(paste0("blavaan ERROR: '", type, "' is only supported for target='stan'"))

    if(type %in% c("yhat", "ypred")) {
      if(is.null(object@external$stanlvs)) stop("blavaan ERROR: for predictions, save.lvs must be TRUE during model estimation")
      lavmcmc <- make_mcmc(blavInspect(object, 'mcobj'), object@external$stanlvs)
      itnums <- sampnums(object@external$mcmcout, thin = 1)
      nsamps <- length(itnums)
      nchain <- length(lavmcmc)
      ng <- blavInspect(object, 'ngroups')

      tmpres <- vector("list", nchain)
      for(j in 1:nchain) {
        loop.args <- list(X = 1:nsamps, FUN = function(i, j){
          cond_moments(lavmcmc[[j]][itnums[i],],
                       blavmodel,
                       blavpartable,
                       blavsamplestats,
                       blavdata,
                       object)}, j = j, future.seed = TRUE)
        tmpres[[j]] <- do.call("future_lapply", loop.args)
      }
      tmpres <- unlist(tmpres, recursive = FALSE)

      if(type == "ypred") {
        ## use mean and cov from each entry of tmpres to randomly sample
        tmpres <- lapply(tmpres, function(x){
          lapply(1:ng, function(g){
            sigchol <- chol(x$cov[[g]])
            t(apply(x$mean[[g]], 1, function(y) mnormt::rmnorm(n=1, mean=y, sqrt=sigchol)))
          })
        })
      } else {
        tmpres <- lapply(tmpres, function(x) x$mean)
      }

      ## these are now lists by group; rearrange to match original data
      cids <- unlist(blavInspect(object, 'case.idx'))
      cnms <- lavNames(object)
      yres <- lapply(tmpres, function(x) do.call("rbind", x)[cids,])
      
      out <- yres
    }

    if(type == "ymis") {
      out <- samp_data(object@external$mcmcout, blavmodel, blavpartable, standata, blavdata)
    }
  } else {
    stop("blavaan ERROR: unknown type supplied; use lv lvmeans yhat ypred ymis")
  }
  
  out
}
