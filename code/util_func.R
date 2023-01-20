make_design_matrix=function(x,grp){
    x2=matrix(nrow=length(x),ncol=length(unique(grp)))
    for(i in 1:length(unique(grp))){
      x2[,i]=ifelse(grp==levels(factor(grp))[i],1,0)*x
    }
  return(x2)
}