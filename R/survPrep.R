survPrep = function(df, timevar='t', statusvar='status'){
   #---------------------------------------------------
   # FUNCTION: Expands data for survival analysis; adds
   #           one row per time point until event, with status=0 until
   #           time of event.  If event occurs at the observed time,
   #           status should be 1 in df.  If not, satus should be 0 in df.
   #---------------------------------------------------
   # INPUTS:   df          = The data frame to be expanded.  Can include variables
   #                          for id, test/train split, etc.
   #           timevar     = The name of the variable containing time info.  Defaults to 't'.
   #           statusvar   = The name of the variable containing status (response) info.
   #                          Defaults to 'status.
   #---------------------------------------------------
   # OUTPUTS: surv         = Expanded survival data frame, with one row per obs per time until event.
   #                          Includes an 'id' variable to identify each observation from orig df.
   #---------------------------------------------------

   require(data.table)

   # Check that all times are > 0.
   if(sum(df[,paste0(timevar)]<=0) > 0){
      stop('Time points exist less than zero.')
   }

   # Add id variable to data set.
   df[,"obs_id"] = 1:nrow(df)

   #---------------------------------------------------
   # Expand data frame.
   #---------------------------------------------------

   nreps = df[,paste0(timevar)]
   surv = setDT(df)
   surv = surv[rep(1:nrow(surv), nreps),]

   #---------------------------------------------------
   # Set response to 0 for all times prior to observed event.
   #---------------------------------------------------

   # Calculate event rows
   surv[,"obs_id_lag1"] = rep(1,nrow(surv))
   surv[2:nrow(surv),"obs_id_lag1"] = surv[1:(nrow(surv)-1),"obs_id"]
   events = c(which(surv[,"obs_id"] != surv[,"obs_id_lag1"]) -1, nrow(surv))

   # Set all non-event responses to zero.
   surv[-events,paste0(statusvar)] = 0

   #---------------------------------------------------
   # Set time var for each obs_id to sequence, up to event.
   #---------------------------------------------------

   t = surv[c(which(surv[,"obs_id"] != surv[,"obs_id_lag1"])-1, nrow(surv)), paste0(timevar)]

   # Add correct time points.
   vecSeq <- Vectorize(seq.default, vectorize.args = c("from", "to"))
   new_times = unlist(vecSeq(1,t))
   surv[,paste0(timevar)] = as.numeric(new_times)

   #---------------------------------------------------
   # Return output. (Drop rownums temp variable during return.)
   #---------------------------------------------------

   # Drop temp rownums variable.
   surv[,"obs_id_lag1"] = NULL
   return(data.frame(surv))

}
