makePredGrid = function(df_new, timevar, timegrid){
   #---------------------------------------------------
   # FUNCTION: Creates expanded data frame df_new (xpred), where
   #           each obs in df appears at all time points
   #           over grid defined by incr.
   #---------------------------------------------------
   # INPUTS:   df_new   = the data frame containing new observations at which
   #                      to predict.  Should not contain time points; if it does,
   #                      they are dropped since all obs are predicted over the grid.
   #           timevar  = name of the time variable, for formatting output data set.
   #           timegrid = grid of timepoints at which to predict.
   #---------------------------------------------------
   # OUTPUTS: predgrid = data frame containing each obs in df_new over the time grid.
   #           Note that all covariates not part of df_orig are dropped.
   #---------------------------------------------------

   require(data.table)

   # Add id variable to data set.
   df_new["obs_id"] = 1:nrow(df_new)

   #---------------------------------------------------
   # Expand data frame.
   #---------------------------------------------------

   nreps = length(timegrid)
   predgrid = setDT(df_new)
   predgrid = predgrid[rep(1:nrow(predgrid), nreps),]

   # Add timevar column.
   predgrid[,paste0(timevar)] = rep(timegrid, each=nrow(df_new))

   return(data.frame(predgrid))
}
