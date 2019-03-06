# From the nnet package.  Creates

makeModelMatrix = function(df){
   #---------------------------------------------------
   # FUNCTION: Creates model matrix for all covariates in df.
   #           Is overparameterized; 1 col for each level
   #           of categorical (factor) variables.
   #---------------------------------------------------
   # INPUTS:   df = A data frame.
   #           Note: If class of column is factor,
   #              is treated as categorical.  Else, treated
   #              as continuous.  (Characters are cast as factors.)
   #---------------------------------------------------
   # OUTPUTS: A matrix containing cols of indicators for cl levels.
   #---------------------------------------------------

   #====================================================================
   # Define the class.ind helper function.
   #====================================================================

   #Description
   #Generates a class indicator function from a given factor.
   #Usage
   #class.ind(cl)
   #Arguments
   #cl factor or vector of classes for cases. Value
   #a matrix which is zero except for the column corresponding to the class.
   #References
   #Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer. Examples

   class.ind <- function(cl){
      #---------------------------------------------------
      # FUNCTION: Creates indicators matrix for a vector
      #           containing a factor.
      #---------------------------------------------------
      # INPUTS:   cl = A vector containing a factor.
      #---------------------------------------------------
      # OUTPUTS: A matrix containing cols of indicators for cl levels.
      #---------------------------------------------------
      n <- length(cl)
      cl <- as.factor(cl)
      x <- matrix(0, n, length(levels(cl)) )
      x[(1:n) + n*(unclass(cl)-1)] <- 1
      dimnames(x) <- list(names(cl), levels(cl))
      x
   }

   #====================================================================
   # Set up model matrix.
   #====================================================================

   X = matrix(NA, nrow=nrow(df), ncol=0)

   # Loop through columns.
   for(j in 1:ncol(df)){

      # Create indicators for column.
      temp = NULL

      if(class(df[,j])=='factor'){
         temp = class.ind(df[,j])
         colnames(temp) = paste0(colnames(df)[j],'.',colnames(temp))
      } else if(class(df[,j])=='character'){
         temp = class.ind(factor(df[,j]))
         colnames(temp) = paste0(colnames(df)[j],'.',colnames(temp))
      } else{
         temp = df[,j, drop=FALSE]
         colnames(temp) = colnames(df[,j,drop=FALSE])
      }

      X = cbind(X,temp)
   }

   # Clean up empty columns, where there were no observed values.
   #X = X[which(colSums(X)>0)]

   return(X)
}
