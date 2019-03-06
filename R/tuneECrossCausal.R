tuneEcrossCausal = function(ecross_candidates = cbind.data.frame(
                                                   'ecross_control'=rep(c(.5,1,2.5,5),each=4),
                                                   'ecross_moderate'=rep(c(.5,1,2.5,5),times=4)),
                    y, pihat, z, tgt, x_control, x_moderate,
                    pihatpred=0, zpred=0, tpred=0, xpred_control=matrix(0,0,0), xpred_moderate=matrix(0,0,0),
                    nburn=100, nsim=1000, ntree_control=200, ntree_moderate=50,
                    lambda=NULL, sigq=.9, sighat=NULL, nu=3,
                    base_control=.95, power_control=2,
                    base_moderate=.25, power_moderate=3,
                    sd_control=2*sd(y), sd_moderate=sd(y),
                    treatment_init = rep(1,length(unique(tgt))),
                    use_muscale=T, use_tauscale=T,
                    pihat_in_trt=F,
                    probit=FALSE, yobs=NULL){

   #---------------------------------------------------
   # FUNCTION: Optimizes expected number of crossings across specified grid.
   #---------------------------------------------------

   require(tidyverse)
   require(ggthemes)

   # Test that candidate list is long enough.
   if(nrow(ecross_candidates)<3) stop('Try at least 3 candidate values for tuning.')

   # Ensure correct col names on candidate df.
   colnames(ecross_candidates) = c('ecross_control','ecross_moderate')

   # Calculate WAIC for each candidate ecross value.
   waic = rep(NA,nrow(ecross_candidates))

   for(i in 1:nrow(ecross_candidates)){

      print(paste0('Iteration ', i, ' of ', nrow(ecross_candidates)))

      # Fit tsbcf model for each ecross candidate.
      fit = tsbcf(y, pihat, z, tgt, x_control, x_moderate,
             pihatpred=0, zpred=0, tpred=0,
             xpred_control=matrix(0,0,0), xpred_moderate=matrix(0,0,0),
             nburn=nburn, nsim=nsim, ntree_control=200, ntree_moderate=50,
             lambda=NULL, sigq=.9, sighat=NULL, nu=3,
             base_control=.95, power_control=2,
             base_moderate=.25, power_moderate=3,
             sd_control=2*sd(y), sd_moderate=sd(y),
             treatment_init = rep(1,length(unique(tgt))),
             use_muscale=T, use_tauscale=T,
             ecross_control=ecross_candidates$ecross_control[i],
             ecross_moderate=ecross_candidates$ecross_moderate[i],
             pihat_in_trt=F,
             probit=FALSE, yobs=NULL, verbose=F, mh=F, save_inputs=F)

         # Calculate in-sample WAIC.
         check = checkFit(y=y,
                          mcmcdraws = fit[["yhat"]],
                          sig = fit[["sigma"]],
                          probit=probit,
                          doWaic=TRUE,
                          yobs=yobs)

      # Save WAIC.
      waic[i] = check$waic
   }

   # Cubic spline fit.  If only tuning one variable, loess fit over that variable only.
   ec = ecross_candidates
   n_candidates_con = length(unique(ec$ecross_control)) # number unique ec_con candidates.
   n_candidates_mod = length(unique(ec$ecross_moderate))# Number unique ec_mod candidates.

   if(n_candidates_con >= 3 & n_candidates_mod >= 3){ #both params being tuned
      myfit = loess(waic ~ ec$ecross_control + ec$ecross_moderate, span=1.25)
   } else if(n_candidates_con >= 3 & n_candidates_mod <3){ # Only control being tuned.
      myfit = loess(waic ~ ec$ecross_control, span=1.25)
   } else{ # Only moderate being tuned.
      myfit = loess(waic ~ ec$ecross_moderate, span=1.25)
   }


   # Calculate sd(resids) from spline fit.
   sd = sd(myfit$residuals)

   # Data frames for calculation.
   edf = cbind.data.frame(ec,waic)
   sdf = cbind.data.frame('x' =ecross_candidates,'y' = myfit$fitted)

   # Calculate optimal ecross.
   ymin = sdf %>% filter(y==min(y)) %>% select(y)

   # Save ecross value where WAIC is minimized.
   df = cbind.data.frame(ec, waic,ymin+sd)
   df$norm = df$ecross_control^2 + df$ecross_moderate^2

   if(length(which(df[,3]<=df[,4]))>0){
      exp_cross = df[which(df[,3]<=df[,4]),]
      exp_cross = exp_cross[which(exp_cross$norm==min(exp_cross$norm)),]
      exp_cross = exp_cross[,1:2]
   } else{
      exp_cross = data.frame(matrix(0,nrow=1,ncol=2))
      colnames(exp_cross) = c('ecross_control','ecross_moderate')
      exp_cross$ecross_control = min(ecross_candidates$ecross_control)
      exp_cross$ecross_moderate = min(ecross_candidates$ecross_moderate)
   }

   if(n_candidates_con >= 3 & n_candidates_mod >= 3){ #both params being tuned

      # Create plot.  Need different plot if only one parameter being tuned.
      df$linesize = 1
      df$linesize[which(df$ecross_moderate==as.numeric(exp_cross$ecross_moderate))] = 2

      waicplt=ggplot(df, aes(x=ecross_control, y=waic, colour=factor(ecross_moderate),
                             linetype=factor(linesize), size=factor(linesize))) +
         geom_line() +
         geom_hline(aes(yintercept=y), colour='grey') +
         geom_vline(aes(xintercept=exp_cross$ecross_control),colour='red',linetype=6) +
         scale_colour_colorblind(name='Ec_moderate') +
         scale_linetype_manual(values=c(1,6),name='Optimal',labels=c('No','Yes'))+
         scale_size_manual(values=c(.8,1.2),name='Optimal',labels=c('No','Yes'), guide=F)

   } else if(n_candidates_con >= 3 & n_candidates_mod <3){ # Only control being tuned.

      waicplt=ggplot(df, aes(x=ecross_control, y=waic)) +
         geom_line() +
         geom_hline(aes(yintercept=y), colour='grey') +
         geom_vline(aes(xintercept=exp_cross$ecross_control),colour='red',linetype=6) +
         scale_colour_colorblind(name='Ec_control')

   } else{ # Only moderate being tuned.
      waicplt=ggplot(df, aes(x=ecross_moderate, y=waic)) +
         geom_line() +
         geom_hline(aes(yintercept=y), colour='grey') +
         geom_vline(aes(xintercept=exp_cross$ecross_moderate),colour='red',linetype=6) +
         scale_colour_colorblind(name='Ec_moderate')
   }



   return(list('ecross_opt' = exp_cross,
               'waic_plot' = waicplt,
               'waic_grid' = cbind.data.frame('ec' = ecross_candidates, waic)))
}
