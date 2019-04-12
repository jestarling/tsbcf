########################################################################
# getCausalEstimands: Function to get posterior causal estimands from tsbcf() object.
#-----------------------------------------------------------------------
# INPUTS: tsbcf_output: A list; the result of the tsbcf() function from tsbart package.
#         subgroups:    An optional n-length vector of subgroups for which to estimate ATE over time, and overall.
#         probit:       A probit indicator.  If yes, transforms counterfactuals to probit scale.
#         relrisk:      Indicator for returning output as a relative risk or an ate on prob scale.
#         indiv:        Indicator for whether to include individual ATE mean/lb/ub estimates, for
#                        in-sample obs, and for out-of-sample obs if predictions are inlucded in tsbcf_output obejct.
#
#-----------------------------------------------------------------------
# OUTPUTS:
# If subgroups = NULL:
#        ate_post:     A nsim-length vector, posterior draws of ATE.
#        ate_t_post:   A (nsim xnt) matrix, posterior draws of ATE at each time (col).
#        ate_t_hat:    A dataframe with one row for each t; gives posterior mean, and credible interval founds.
# If subgroups != NULL
#        ate_post:    A (nsim x ngrps) matrix, posterior draws of CATE for each subgroup.
#        ate_t_post:  A (nsim x ngrps*nt) matrix, posterior draws of CATE for each subgroup/time combo.
#        ate_t_hat:    A dataframe with one row for each t/subgroup combo; gives posterior mean, and credible interval founds.
# If indiv = T
#        ate_indiv:    A dataframe with n rows and 4 cols: obs number, posterior mean indivual ate, lb, ub.
#        ate_indiv_oos: A dataframe with npred rows and 4 cols: oos obs number, posterior mean indivual ate, lb, ub.
#
#-----------------------------------------------------------------------
# Note about probit outcomes:
#  For probit outcome, the counterfactual probabilities are:
#  w_{it}(0) = Phi(mu(x,t))
#  w_{it}(1) = Phi(mu(x,t) + tau(x,t))
#  Treatment Effect = w_{it}(1) - w_{it}(0)
#  Relative Risk = w_{it}(1) / w_{it}(0) where Phi() = the standard normal cdf, ie pnorm.
########################################################################

getCausalEstimands = function(tsbcf_output, probit=F, relrisk=F, indiv=F, subgroups=NULL, subgroups_pred=NULL){

   # Error checking.
   if(!is.null(subgroups)){
      if(length(subgroups)!=ncol(tsbcf_output$tau)){
         stop('subgroups must be an n-length vector.')
      }
   }

   if(!is.null(subgroups)){
      if(length(subgroups)!=ncol(tsbcf_output$tau_oos)){
         stop('subgroups must be an npred-length vector.')
      }
   }

   if(relrisk==T & probit==F){
      warning('relrisk=T ignored unless probit=T.')
   }

   # Set up grids of times and groups.
   tgrid = sort(unique(tsbcf_output$tgt))
   grps = NA
   if(!is.null(subgroups)){
      grps = sort(unique(subgroups))
   }

   # Number of times and groups.
   nt = length(tgrid)
   ng = length(unique(subgroups))

   # Empty list for output.
   out = list()

   # Indicator for whether out-of-sample data is included in the fit.
   pred = 0
   if("tau_oos" %in% names(tsbcf_output)){
      pred=1
   }

   ########################################################################
   # Get counterfactuals, and if probit, transform.
   ########################################################################

   #-----------------------------------------------------------------------
   # In-Sample.
   #-----------------------------------------------------------------------
   if(probit==FALSE){
      my_matrix = tsbcf_output$tau

   } else{
      w1  =  pnorm(tsbcf_output$mu + tsbcf_output$tau)
      w0 =  pnorm(tsbcf_output$mu)

      if(relrisk==FALSE){
         my_matrix = w1-w0
      } else{
         my_matrix = w1/w0
      }
   }

   #-----------------------------------------------------------------------
   # Out-of-Sample.
   #-----------------------------------------------------------------------
   if(pred){

      if(probit==FALSE){
         my_matrix_oos = tsbcf_output$tau_oos

      } else{
         w1_oos  =  pnorm(tsbcf_output$mu_oos + tsbcf_output$tau_oos)
         w0_lls =  pnorm(tsbcf_output$mu_oos)

         if(relrisk==FALSE){
            my_matrix_oos = w1_oos-w0_oos
         } else{
            my_matrix_oos = w1_oos/w0_oos
         }
      }

   }

   ########################################################################
   # No Subgroups Case:  Calculate ATE and ATE_t.
   ########################################################################
   if(is.null(subgroups)){

      #-----------------------------------------------------------------------
      # a. Posterior draws for overall ATE.
      #    Sets up nsim-length vector of posterior ate draws.
      #-----------------------------------------------------------------------

      # In-sample.
      ate = apply(my_matrix,1,function(x) mean(x,na.rm=T)) #rowMeans(my_matrix)
      out$ate_post = ate

      # Out-of-sample.
      if(pred){
         ate_oos = apply(my_matrix_oos,1,function(x) mean(x,na.rm=T)) #rowMeans(my_matrix)
         out$ate_post_oos = ate_oos
      }

      #-----------------------------------------------------------------------
      # b. Posterior draws for ATE at each time.
      #    Sets up nsim x length(tgrid) matrix of posterior ate draws for each tgt value (column).
      #-----------------------------------------------------------------------

      # In-sample.
      ate_t = as.data.frame(matrix(0,nrow=nrow(my_matrix),ncol=length(tgrid)))
      colnames(ate_t) = tgrid

      for(t in 1:length(tgrid)){
         ate_t[,t] = rowMeans(my_matrix[,which(tsbcf_output$tgt==tgrid[t])])
      }

      out$ate_t_post = ate_t

      # Out-of-sample.
      if(pred){

         ate_t_oos = as.data.frame(matrix(0,nrow=nrow(my_matrix_oos),ncol=length(tgrid)))
         colnames(ate_t_oos) = tgrid

         for(t in 1:length(tgrid)){
            ate_t_oos[,t] = rowMeans(my_matrix_oos[,which(tsbcf_output$tgt==tgrid[t])])
         }

         out$ate_t_post_oos = ate_t_oos

      }

      ########################################################################
      # c. Posterior means and credible intervals for ATE at each time.
      ########################################################################

      # In-sample.
      out$ate_t_hat = cbind.data.frame('tgt'=tgrid,
                                       'ate'=apply(out$ate_t_post,2,function(x) mean(x,na.rm=T)),
                                       'lb'=apply(out$ate_t_post,2,function(x) quantile(x,.025,na.rm=T)),
                                       'ub'=apply(out$ate_t_post,2,function(x) quantile(x,.975,na.rm=T)))

      # Out-of-sample.
      if(pred){
         out$ate_t_hat_oos = cbind.data.frame('tgt'=tgrid,
                                              'ate'=apply(out$ate_t_post_oos,2,function(x) mean(x,na.rm=T)),
                                              'lb'=apply(out$ate_t_post_oos,2,function(x) quantile(x,.025,na.rm=T)),
                                              'ub'=apply(out$ate_t_post_oos,2,function(x) quantile(x,.975,na.rm=T)))
      }

   }

   #-----------------------------------------------------------------------
   # SUBGROUPS: CATE and CATE_t
   #-----------------------------------------------------------------------
   if(!is.null(subgroups)){

      #-----------------------------------------------------------------------
      # a. Posterior draws for overall CATE (using subgroups).
      #-----------------------------------------------------------------------

      # In-sample.
      cate = as.data.frame(matrix(0,nrow=nrow(my_matrix),ncol=ng))
      colnames(cate) = grps

      for(g in 1:ng){
         cate[,g] =  rowMeans(my_matrix[,which(subgroups==grps[g])])
      }

      out$ate_post = cate

      # Out-of-sample.
      if(pred){

         cate_oos = as.data.frame(matrix(0,nrow=nrow(my_matrix_oos),ncol=ng))
         colnames(cate_oos) = grps

         for(g in 1:ng){
            cate_oos[,g] =  rowMeans(my_matrix_oos[,which(subgroups_oos==grps[g])])
         }

         out$ate_post_oos = cate_oos
      }

      #-----------------------------------------------------------------------
      # . Posterior draws for CATE at each time (using subgroups.)
      #-----------------------------------------------------------------------

      # In-sample.
      cate_t = as.data.frame(matrix(0,nrow=nrow(my_matrix),ncol=ng*nt))

      groups_and_times = expand.grid(grps,tgrid)
      colnames(cate_t) = do.call(paste0, expand.grid(grps,'-',tgrid))

      for(tg in 1:(ng*nt)){
         temp = my_matrix[,which(subgroups==groups_and_times[tg,1]
                                 & tsbcf_output$tgt==groups_and_times[tg,2])]
         cate_t[,tg] =  apply(temp,1,function(x) mean(x,na.rm=T))
      }

      out$ate_t_post = cate_t

      # Out-of-sample.
      if(pred){

         cate_t_oos = as.data.frame(matrix(0,nrow=nrow(my_matrix_oos),ncol=ng*nt))

         groups_and_times = expand.grid(grps,tgrid)
         colnames(cate_t_oos) = do.call(paste0, expand.grid(grps,'-',tgrid))

         for(tg in 1:(ng*nt)){
            temp = my_matrix_oos[,which(subgroups==groups_and_times[tg,1]
                                    & tsbcf_output$tgt_oos==groups_and_times[tg,2])]
            cate_t_oos[,tg] =  apply(temp,1,function(x) mean(x,na.rm=T))
         }

         out$ate_t_post_oos = cate_t_oos

      }

      #-----------------------------------------------------------------------
      # c. Posterior means and credible intervals for ATE at each time.
      #-----------------------------------------------------------------------

      # In-sample.
      out$ate_t_hat = cbind.data.frame('tgt'=rep(tgrid, each=ng),
                                       'subgroup'=rep(grps, times=nt),
                                       'ate'=apply(out$ate_t_post,2,function(x) mean(x,na.rm=T)),
                                       'lb'=apply(out$ate_t_post,2,function(x) quantile(x,.025,na.rm=T)),
                                       'ub'=apply(out$ate_t_post,2,function(x) quantile(x,.975,na.rm=T)))

      # Out-of-sample.
      if(pred){
         out$ate_t_hat_oos = cbind.data.frame('tgt'=rep(tgrid, each=ng),
                                          'subgroup'=rep(grps, times=nt),
                                          'ate'=apply(out$ate_t_post_oos,2,function(x) mean(x,na.rm=T)),
                                          'lb'=apply(out$ate_t_post_oos,2,function(x) quantile(x,.025,na.rm=T)),
                                          'ub'=apply(out$ate_t_post_oos,2,function(x) quantile(x,.975,na.rm=T)))
      }

   }

   ########################################################################
   # Individual ATE estimates/intervals.
   ########################################################################
   if(indiv){

      # Estimate in-sample individual treatment effects.
      out$ate_indiv = cbind.data.frame(
         'obs'=1:ncol(my_matrix),
         'mean'=apply(my_matrix,2,function(x) mean(x,na.rm=T)),
         'lb'=apply(my_matrix,2,function(x) quantile(x,.025,na.rm=T)),
         'ub'=apply(my_matrix,2,function(x) quantile(x,.975,na.rm=T)))

      # Estimate out-of-sample individual treatment effects.
      if(pred){
         out$ate_indiv_oos = cbind.data.frame(
            'obs'=1:ncol(my_matrix_oos),
            'mean'=apply(my_matrix_oos,2,function(x) mean(x,na.rm=T)),
            'lb'=apply(my_matrix_oos,2,function(x) quantile(x,.025,na.rm=T)),
            'ub'=apply(my_matrix_oos,2,function(x) quantile(x,.975,na.rm=T)))
      }

   }

   ########################################################################
   # Return function output.
   ########################################################################
   return(out)
}
