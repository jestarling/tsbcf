tuneEcross = function(ecross_candidates = seq(.25,5,by=.25),
                      y, tgt, tpred, x, xpred, nburn=1000, nsim=1000, ntree=200,
                      lambda=NULL, sigq=.9, sighat=NULL, nu=3,
                      base_tree=.95, power_tree=2,
                      probit=FALSE, yobs=NULL){

   #---------------------------------------------------
   # FUNCTION: Optimizes expected number of crossings across specified grid.
   #---------------------------------------------------

   # Test that candidate list is long enough.
   if(length(ecross_candidates)<3) stop('Try at least 3 candidate values for tuning.')

   # Calculate WAIC for each candidate ecross value.
   waic = rep(NA,length(ecross_candidates))

   for(i in 1:length(ecross_candidates)){

      print(paste0('Iteration ', i, ' of ', length(ecross_candidates)))

      # Fit tsBART model for each ecross candidate.
      fit = tsbart(y=y, tgt=tgt, x=x, tpred=0, xpred=matrix(0,0,0), nburn=nburn, nsim=nsim, ntree=ntree,
                   lambda=lambda, sigq=sigq, sighat=sighat, nu=nu,
                   ecross=ecross_candidates[i], base_tree=base_tree, power_tree=power_tree, use_fscale=T, #sd_control=1,
                   probit=probit, yobs=yobs, verbose=F)


      # Calculate in-sample WAIC.
      check = checkFit(y=y,
                          mcmcdraws = fit[["mcmcdraws"]],
                          sig = fit[["sigma"]],
                          probit=probit,
                          doWaic=TRUE,
                          yobs=yobs)

      # Save WAIC.
      waic[i] = check[["waic"]]
   }

   # Cubic spline fit.
   ec = ecross_candidates
   myfit = loess(waic ~ ec, span=1.25)

   # Calculate sd(resids) from spline fit.
   sd = sd(myfit[["residuals"]],na.rm=T)

   # Data frames for calculation.
   edf = cbind.data.frame(ec,waic)
   sdf = cbind.data.frame('x' = ecross_candidates,'y' = myfit[["fitted"]])

   # Calculate optimal ecross.
   ymin = sdf[which(sdf$y==min(sdf$y)),2]

   # Save ecross value where WAIC is minimized.
   df = cbind.data.frame(ec, waic,ymin+sd)
   exp_cross = ecross_candidates[min(which(df[,2] < df[,3]))]

   if(is.null(exp_cross)){
      exp_cross = min(ecross_candidates)
   }

   # Create plot
   waicplt = ggplot() +
      geom_line(data=edf, aes(x=ec,y=waic, colour='waic',linetype='waic'), size=.8) +
      geom_line(data=sdf, aes(x=ec,y=y, colour='smoothfit', linetype='smoothfit'), size=1.2) +
      geom_line(data=sdf, aes(x=ec, y=ymin+sd, colour='bound', linetype='bound'), size=.8) +

      scale_colour_manual(values = c('waic' = 'grey20', 'smoothfit' = 'dodgerblue4','bound' = 'purple'),
                          labels=c('waic'='WAIC','smoothfit'='Smooth Fit','bound'='Upper bound'), name='') +
      scale_linetype_manual(values=c('waic'=2, 'smoothfit'=1,'bound'=3),
                            labels=c('waic'='WAIC','smoothfit'='Smooth Fit','bound'='Upper bound'), name='') +
      geom_point(data=edf, aes(x=ec,y=waic), colour='black') +
      geom_vline(aes(xintercept=exp_cross), colour='black', show.legend=FALSE) +
      labs(x = 'Candidate expected crossings',
           y = 'WAIC',
           title = paste0('Optimal expected crossings: ', exp_cross)) +
      theme(legend.key.size = unit(0.5,"in"))

   # Return function output.
   return(list('ecross_opt' = exp_cross,
               'waic_plot' = waicplt,
               'waic_grid' = cbind.data.frame('ec' = ecross_candidates, waic)))
}
