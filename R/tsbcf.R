# Fit tsbcf.  Wrapper function for calling tsbartFit.cpp and tsbartProbit.cpp.

### For argument validation.
.ident <- function(...){
   # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
   args <- c(...)
   if( length( args ) > 2L ){
      #  recursively call ident()
      out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
   }else{
      out <- identical( args[1] , args[2] )
   }
   return( all( out ) )
}

### Main tsbcf wrapper function.
tsbcf <- function(y, pihat, z, tgt, x_control, x_moderate,
                  pihatpred=NULL,
                  zpred=NULL,
                  tpred=NULL,
                  xpred_control=NULL, xpred_moderate=NULL,
                  nburn=100, nsim=1000,
                  ntree_control=200, ntree_moderate=50,
                  lambda=NULL, sigq=.9, sighat=NULL, nu=3,
                  base_control=.95, power_control=2,
                  base_moderate=.25, power_moderate=3,
                  sd_control=2*sd(y), sd_moderate=sd(y),
                  treatment_init = rep(1,length(unique(tgt))),
                  use_muscale=T, use_tauscale=T,
                  ecross_control=1, ecross_moderate=1,
                  ecross_control_candidates = NULL,
                  ecross_moderate_candidates = NULL,
                  ecross_tune_nsim=100, ecross_tune_nburn=1000,
                  pihat_in_trt=F,
                  probit=FALSE, yobs=NULL, set_probit_scales=F,
                  verbose=T, mh=F, save_inputs=T){

   ################################################################
   # Capture key arguments.
   ################################################################

   options(expressions=10000)

   inputs = cbind.data.frame(
      'arg'=c('nburn','nsim','ntree_control','ntree_moderate','lambda','sigq','sighat','nu','base_control','power_control',
              'base_moderate','power_moderate','sd_control','sd_moderate','treatment_init','use_muscale','use_tauscale',
              'ecross_control','ecross_moderate','pihat_in_trt','probit','verbose','mh'),
      'value'=c(nburn,nsim,ntree_control,ntree_moderate,ifelse(is.null(lambda),"NULL",lambda),sigq,ifelse(is.null(sighat),"NULL",sighat),nu,base_control,power_control,
                base_moderate,power_moderate,sd_control,sd_moderate, paste0(treatment_init,'',collapse=','), use_muscale,use_tauscale,
                ecross_control,ecross_moderate,pihat_in_trt,probit,verbose,mh)
   )

   ################################################################
   # Validate inputs.
   ################################################################

   #---------------------------------------------------------------
   # If not predicting, set pred objects to first three obs.
   # These are not output.
   #---------------------------------------------------------------
   predict = 1

   # If some, but not all, prediction objects supplied - stop!
   if(sum(is.null(pihatpred), is.null(zpred), is.null(tpred), is.null(xpred_control), is.null(xpred_moderate)) %in% c(1:5)){

      if(sum(is.null(pihatpred), is.null(zpred), is.null(tpred), is.null(xpred_control), is.null(xpred_moderate))==5){
         predict = 0
         pihatpred = pihat[1:min(3,length(pihat))]
         zpred = z[1:min(3,length(z))]
         tpred = tgt[1:min(3,length(tgt))]
         xpred_control = x_control[1:min(3,nrow(x_control)),,drop=F]
         xpred_moderate = x_moderate[1:min(3,nrow(x_moderate)),,drop=F]
      }
      else{
         stop('When supplying data for prediction, must input all of the following:\n',
              'tpred, zpred, pihatpred, xpred_control, xpred_moderate.\n')
      }
   }

   #---------------------------------------------------------------
   # Data size.
   #---------------------------------------------------------------

   # Check data size match.
   if( !.ident(length(y), length(pihat), length(tgt), nrow(x_control), nrow(x_moderate))){

      stop("Data size mismatch. The following should all be equal:
           length(y): ", length(y), "\n",
           "length(pihat): ", length(pihat), "\n",
           "length(tgt): ", length(tgt), "\n",
           "nrow(x_control): ", nrow(x_control), "\n",
           "nrow(x_moderate): ", nrow(x_moderate), "\n")
   }

   # Check out-of-sample data size match.
   if( !.ident(length(tpred), nrow(xpred_control), nrow(xpred_moderate), length(pihatpred))){

      stop("Data size mismatch. The following should all be equal:
           length(tpred): ", length(tpred), "\n",
           "length(pihatpred): ", length(pihatpred), "\n",
           "nrow(xpred_control): ", nrow(xpred_control), "\n",
           "nrow(xpred_control): ", nrow(xpred_control), "\n")
   }

   # Check trt_init length.
   if( !.ident(length(treatment_init), length(unique(tgt)))){
      stop('tgt_init must have length matching unique number of values in tgt.')
   }

   #---------------------------------------------------------------
   # Probit checks.
   #---------------------------------------------------------------
   if(probit==TRUE){

      # Data size match including yobs.
      if( !.ident(length(y), length(yobs), length(tgt), nrow(x_moderate), nrow(x_control))){
         stop("Data size mismatch. The following should all be equal:
              length(y): ", length(y), "\n",
              "length(yobs): ", length(yobs), "\n",
              "length(tgt): ", length(tgt), "\n",
              "nrow(x_control): ", nrow(x_control), "\n",
              "nrow(x_moderate): ", nrow(x_moderate), "\n")
      }

      #Yobs must be only 0/1.  Y must not be only 0/1.
      if(length(unique(y))>2) warning("In probit case,
                                      y should contain initial latent variables,
                                      and yobs should contain observed binary values.")

      if(length(unique(yobs))>2) stop("In probit case,
                                      y should contain initial latent variables,
                                      and yobs should contain observed binary values.")
      if(is.null(yobs)) stop("yobs must be populated when probit=TRUE, and should contain observed binary values of 0/1.")

      # Warn user that manually input lambda and sighat values are ignored in probit case.
      if(!is.null(lambda) || !is.null(sighat)) warning("lambda and sighat inputs are ignored in probit case, as prior
                                                       tuning for sigma^2 is not applicable.")
   }

   if(probit==FALSE){
      if(!is.null(yobs)) stop("yobs is only for probit=TRUE case. Must be NULL when probit=FALSE.")
   }

   #---------------------------------------------------------------
   # Other inputs.
   #---------------------------------------------------------------

   if(any(is.na(y))) stop("Missing values in y")
   if(any(is.na(tgt))) stop("Missing values in tgt")
   if(any(is.na(tpred))) stop("Missing values in tpred")
   if(any(is.na(x_control))) stop("Missing values in x_control")
   if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
   if(any(is.na(xpred_control))) stop("Missing values in xpred_control")
   if(any(is.na(xpred_moderate))) stop("Missing values in xpred_moderate")

   if(length(unique(y))<5 && probit==FALSE) warning("y appears to be discrete")
   # Ok for probit bc might be initializing 0/1 to two latent values.

   if(any(pihat>1 || pihat<0)) stop("pihat must be in 0-1 range.")
   if(any(pihatpred>1 || pihatpred<0)) stop("pihatpred must be in 0-1 range.")

   if(nburn<0) stop("nburn must be positive")
   if(nsim<0) stop("nsim must be positive")
   if(ecross_tune_nburn<0) stop("ecross_tune_nburn must be positive")
   if(ecross_tune_nsim<0) stop("ecross_tune_nsim must be positive")
   if(ecross_control<=0) stop("ecross_control must be positive")
   if(class(ecross_control)=="character" & ecross_control!="tune") stop("ecross_control must be a positive value or set to 'tune'.")
   if(ecross_moderate<=0) stop("ecross_moderate must be positive")
   if(class(ecross_moderate)=="character" & ecross_moderate!="tune") stop("ecross_moderate must be a positive value or set to 'tune'.")

   # For tuning parameter inputs.
   if(ecross_control!="tune" & !is.null(ecross_control_candidates)) warning("ecross_control_candidates will be ignored; ecross_control value has been provided.")
   if(ecross_moderate!="tune" & !is.null(ecross_moderate_candidates)) warning("ecross_moderate_candidates will be ignored; ecross_moderate value has been provided.")

   if(ecross_control=="tune" & is.null(ecross_control_candidates)) warning("Default ecross_control_candidates will be used: {1, 2.5, 5}")
   if(ecross_moderate=="tune" & is.null(ecross_moderate_candidates)) warning("Default ecross_moderate_candidates will be used: {1, 2.5, 5, 7.5, 10}")

   ################################################################
   # Scale/center response y. (For non-probit case.)
   ################################################################
   ybar = ifelse(probit==FALSE, mean(y), 0)
   sdy = ifelse(probit==FALSE, sd(y), 1)
   y_scale = (y - ybar) / sdy

   ################################################################
   # Create model matrix and set up hyperparameters.
   ################################################################

   # Add propensity score to x_control and xpred_control matrices.
   x_control[,ncol(x_control)+1] = pihat
   xpred_control[,ncol(xpred_control)+1] = pihatpred

   # If indicated by user, add propensity score to x_moderate and xpred_moderate matrices.
   if(pihat_in_trt==TRUE){
      x_moderate[,ncol(x_moderate)+1] = pihat
      xpred_moderate[,ncol(xpred_moderate)+1] = pihatpred
   }

   # Model matrices.
   xx_control = tsbart::makeModelMatrix(x_control)
   xxpred_control = tsbart::makeModelMatrix(xpred_control)
   cutpoints_control = tsbart::makeCutpoints(xx_control)

   xx_moderate = tsbart::makeModelMatrix(x_moderate)
   xxpred_moderate = tsbart::makeModelMatrix(xpred_moderate)
   cutpoints_moderate = tsbart::makeCutpoints(xx_moderate)

   # Sighat and lambda calibration.
   if(is.null(sighat)){
      df = cbind.data.frame(y_scale, tgt, xx_control)
      lmf = lm(y_scale ~ ., data=df)
      sighat = sigma(lmf)
   }

   if(is.null(lambda)){
      qchi = qchisq(1-sigq, nu)
      lambda = (sighat * sighat * qchi) / nu
   }

   ################################################################
   # Set up probit parameters if probit=T.
   ################################################################
   offset=0

   if(probit==TRUE){
      phat = mean(unlist(yobs))
      offset = qnorm(phat)
   }

   ################################################################
   # Probit-scaled default control_sd and moderate_sd, based on
   # estimates of baseline risk and relative risk in data, if indicated.
   ################################################################
   if(probit==TRUE & set_probit_scales==T){

      phat = sum(yobs==1)/length(yobs)
      rrhat = (sum(yobs==1 & z==1) / sum(z==1)) /
         (sum(yobs==1 & z==0) / sum(z==0))

      sd_control = abs(qnorm(phat))
      sd_moderate = abs(tau_calc(phat, rrhat))

      print(paste0('Setting probit scales:'))
      print(paste0('sd_control: ', sd_control, ', sd_moderate: ', sd_moderate))

      # Update inputs.
      levels(inputs$value) <- c(levels(inputs$value), sd_control, sd_moderate)
      inputs$value[which(inputs$arg=="sd_control")] = sd_control
      inputs$value[which(inputs$arg=="sd_moderate")] = sd_moderate

   }

   ################################################################
   # Set up ordering, so that z=1 is first in in-samp and out-of-samp datasets.
   ################################################################

   perm = order(z, decreasing=TRUE)
   perm_oos = order(zpred, decreasing=TRUE)

   ################################################################
   # Parameter tuning if necessary.
   ################################################################

   # Set up ecross_candidates for tuning.
   if(ecross_control=="tune" & is.null(ecross_control_candidates)){
      ecross_control_candidates = c(1,2.5,5)
   }
   if(ecross_control!="tune"){
      ecross_control_candidates = ecross_control
   }

   if(ecross_moderate=="tune" & is.null(ecross_moderate_candidates)){
      ecross_moderate_candidates = c(1,2.5,5,7.5,10)
   }
   if(ecross_moderate!="tune"){
      ecross_moderate_candidates = ecross_moderate
   }


   # Perform tuning.
   if(ecross_control=="tune" || ecross_moderate=="tune"){

      tuned = tuneEcrossCausal(ecross_control_candidates=ecross_control_candidates,
                               ecross_moderate_candidates=ecross_moderate_candidates,
                               y, pihat, z, tgt, x_control, x_moderate,
                               pihatpred, zpred, tpred, xpred_control, xpred_moderate,
                               nburn=ecross_tune_nburn, nsim=ecross_tune_nsim, ntree_control, ntree_moderate,
                               lambda, sigq, sighat, nu,
                               base_control, power_control, base_moderate, power_moderate,
                               sd_control, sd_moderate, treatment_init,
                               use_muscale, use_tauscale, pihat_in_trt, probit, yobs)

      ecross_control = as.numeric(tuned$ecross_control)
      ecross_moderate = as.numeric(tuned$ecross_moderate)
   }

   ################################################################
   # Set up some info about t; vector of counts for each time point (nt).
   ################################################################
   nt = as.numeric(table( factor(tgt, levels = min(c(tgt,tpred)):max(tgt,tpred))))
   ntpred = as.numeric(table( factor(tpred, levels = min(c(tgt,tpred)):max(tgt,tpred))))

   # Lookup values for matching
   tgrid = unique(sort(c(tgt,tpred)))
   tgt_idx = match(tgt, tgrid)
   tpred_idx = match(tpred,tgrid)


   ################################################################
   # Call tsbartFit.cpp or tsbartProbit.cpp
   ################################################################
   out = NULL

   if(probit==FALSE){
      out = tsbcfFit(y = y_scale[perm], z = z[perm], zpred = zpred[perm_oos], tgt = tgt[perm],
                     tpred = tpred[perm_oos],
                     x_con = t(xx_control[perm,]), x_mod = t(xx_moderate[perm,]),
                     xpred_con = t(xxpred_control[perm_oos,]), xpred_mod = t(xxpred_moderate[perm_oos,]),
                     xinfo_list_con = cutpoints_control, xinfo_list_mod = cutpoints_moderate,
                     nburn = nburn, nsim = nsim, ntree_con = ntree_control, ntree_mod = ntree_moderate,
                     lambda=lambda, sigq=sigq, sigma=sighat, nu=nu,
                     base_con=base_control, power_con=power_control,
                     base_mod=base_moderate, power_mod=power_moderate,
                     ecross_con=ecross_control, ecross_mod=ecross_moderate,
                     con_sd=sd_control, mod_sd=sd_moderate,
                     trt_init = treatment_init,
                     use_muscale=use_muscale, use_tauscale=use_tauscale,
                     treef_name_="tsbtrees.txt", save_trees=FALSE, silent_mode=!verbose)
   } else{

      out =  tsbcfProbit(y = y_scale[perm], yobs=yobs[perm], z = z[perm], zpred = zpred[perm_oos], tgt = tgt[perm],
                         tpred = tpred[perm_oos],
                         x_con = t(xx_control[perm,]), x_mod = t(xx_moderate[perm,]),
                         xpred_con = t(xxpred_control[perm_oos,]), xpred_mod = t(xxpred_moderate[perm_oos,]),
                         xinfo_list_con = cutpoints_control, xinfo_list_mod = cutpoints_moderate,
                         nburn = nburn, nsim = nsim, ntree_con = ntree_control, ntree_mod = ntree_moderate,
                         offset=offset,
                         lambda=lambda, sigq=sigq, nu=nu,
                         base_con=base_control, power_con=power_control,
                         base_mod=base_moderate, power_mod=power_moderate,
                         ecross_con=ecross_control, ecross_mod=ecross_moderate,
                         con_sd=sd_control, mod_sd=sd_moderate,
                         trt_init = treatment_init,
                         use_muscale=use_muscale, use_tauscale=use_tauscale,
                         treef_name_="tsbtrees.txt", save_trees=FALSE, silent_mode=!verbose)
   }

   ################################################################
   # Rescale output and restore to correct order.
   ################################################################

   # In-sample
   yhat = ybar + out$yhat[,order(perm)] * sdy
   mu = ybar + out$mu[,order(perm)] * sdy
   tau = out$tau[,order(perm)] * sdy
   sig_rescaled = out$sigma * sdy

   # Out-of-sample
   yhat_oos = ybar + out$yhat_oos[,order(perm_oos)] * sdy
   mu_oos = ybar + out$mu_oos[,order(perm_oos)] * sdy
   tau_oos = out$tau_oos[,order(perm_oos)] * sdy

   ################################################################
   # Adjust alpha's and add accept/reject indicator.
   # Note: bd function returns:
   #     alpha in (0,1) for births which are accepted
   #     -alpha in (-1,0) for deaths which are accepted
   #     10+alpha for rejected births
   #     -alpha-10 for rejected deaths.
   ################################################################

   if(mh){
      # Con metropolis info.
      bd_con = ifelse(out$alpha_con<=0,0,1)             # 1 = birth, 0 = death
      accept_con = ifelse(abs(out$alpha_con)<10,1,0)   # 1 = accepted, 0 = rejected
      alpha_con = ifelse(accept_con==1, abs(out$alpha_con), abs(out$alpha_con)-10)

      # Mod metropolis info.
      bd_mod = ifelse(out$alpha_mod<=0,0,1)             # 1 = birth, 0 = death
      accept_mod = ifelse(abs(out$alpha_mod)<10,1,0)   # 1 = accepted, 0 = rejected
      alpha_mod = ifelse(accept_mod==1, abs(out$alpha_mod), abs(out$alpha_mod)-10)

      # Assemble dataframes and convert bd to character.
      # Note: alpha is (nburn+nsim x ntree).
      metrop_con = cbind.data.frame(
         'iter' = rep(1:(nburn+nsim), times=ntree_control),
         'tree' = rep(1:ntree_control, each=nburn+nsim),
         'accept' = as.numeric(accept_con),
         'alpha' = as.numeric(alpha_con),
         'bd' = as.numeric(bd_con)
      )

      metrop_mod = cbind.data.frame(
         'iter' = rep(1:(nburn+nsim), times=ntree_moderate),
         'tree' = rep(1:ntree_moderate, each=nburn+nsim),
         'accept' = as.numeric(accept_mod),
         'alpha' = as.numeric(alpha_mod),
         'bd' = as.numeric(bd_mod)
      )

      metrop_con$bd = ifelse(metrop_con$bd==1,'birth','death')
      metrop_mod$bd = ifelse(metrop_mod$bd==1,'birth','death')

      # Combine into one metrop dataframe.
      metrop = rbind.data.frame(metrop_con, metrop_mod)
      metrop$tree = c(rep("control",nrow(metrop_con)), rep("moderate",nrow(metrop_mod)))
      rm(metrop_con); rm(metrop_mod)
   }

   ################################################################
   # Return output.
   ################################################################

   # Only include out-of-sample info if indicated by user.
   if(predict){
      out = list('tgt'=tgt,
                 'yhat'=yhat,
                 'mu'=mu,
                 'tau'=tau,
                 'tpred'=tpred,
                 'yhat_oos'=yhat_oos,
                 'mu_oos'=mu_oos,
                 'tau_oos'=tau_oos,
                 'sigma'=sig_rescaled,
                 'mu_sd' = out$mu_sd_post,
                 'tau_sd' = out$tau_sd_post,
                 #'mscale' = out$mscale,
                 #'bscale1'=out$bscale1,
                 #'bscale0' = out$bscale0,
                 #'alpha_t'=out$alpha_t*sdy + ybar,
                 'ecross_control'=ecross_control,
                 'ecross_moderate'=ecross_moderate)
   } else{
      out = list('tgt'=tgt,
                 'yhat'=yhat,
                 'mu'=mu,
                 'tau'=tau,
                 'sigma'=sig_rescaled,
                 'mu_sd' = out$mu_sd_post,
                 'tau_sd' = out$tau_sd_post,
                 #'mscale' = out$mscale,
                 #'bscale1'=out$bscale1,
                 #'bscale0' = out$bscale0,
                 'ecross_control'=ecross_control,
                 #'alpha_t'=out$alpha_t*sdy + ybar,
                 'ecross_moderate'=ecross_moderate)
   }

   # Include metropolis info if indicated.
   if(mh){
      out$metrop = metrop
   }

   # Include inputs if indicated.
   if(save_inputs){
      out$inputs = inputs
   }

   # Return output.
   return(out)
   }
