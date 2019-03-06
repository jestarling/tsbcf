# This function sorts on z, binary (0/1) treatment variable.  1's first, 0's second.

dfSort = function(df, trt_var){
   perm = rev(order(df[paste0(trt_var)]))
   return(df[perm,])
}

