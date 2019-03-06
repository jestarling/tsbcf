make_lookup_ind = function(cx, type) {
  pos = 1
  if(type[1]<3) {
    ta = cbind(1:1, c(pos,pos))
    pos = pos+1
  } else {
    ta = cbind(rep(1,cx[1]), pos:(pos+cx[1]-1))
    pos = pos+cx[1]
  }
  
  
  for(j in 2:length(cx)) {
    if(type[j]<3) {
      ta = rbind(ta, cbind(c(j,j), c(pos,pos)))
      pos = pos+1
    } else {
      ta = rbind(ta, cbind(rep(j,cx[j]), pos:(pos+cx[j]-1)))
      pos = pos+cx[j]
    }
  }
  
  ta
}

make_ind = function(X, cx, type) {
  if(type[1] <3) {
    out = X[,1]
  } else {
    out = as.numeric(X[,1]==0)
    for(k in 2:cx[1]) {
      out = cbind(out, as.numeric(X[,1]==(k-1)))
    }
  }
  
  for(j in 2:ncol(X)) {
    if(type[j]<3) {
      out = cbind(out, X[,j])
    } else {
      out = cbind(out, as.numeric(X[,j]==0))
      for(k in 2:cx[j]) {
        out = cbind(out, as.numeric(X[,j]==(k-1)))
      }
    }
  }
  out
}