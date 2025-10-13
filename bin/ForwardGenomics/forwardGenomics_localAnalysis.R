# Xavier Prudent, 2015

######################################################
## Initialization of lists
######################################################
#logcat <- function(..., file = "logfile.log", append = TRUE) {
#  cat(..., file = file, append = append)
#}

listInitial = function(){

  ## Distance to the expected value
  totBranchesDist  <<- c()
  neutBranchesDist <<- c()
  selBranchesDist  <<- c()

  ## Weight 
  brWeights <<- c()
  neutBranchesWeights <<- c()

  ## Status
  selBranchesStatus <<- c()
  brStatus <<- c()

  ## For the breanch-per-branch method
  ## Distance to the expected value
  BpB.totBranchesDist<<- c()
  BpB.neutBranchesDist <<- c()

  ## Weight 
  BpB.brWeights <<- c()
  BpB.neutBranchesWeights <<- c()

  ## Status
  BpB.brStatus <<- c()

}


######################################################
## branch-based analysis
######################################################
treeBased_analysis = function( elID, thisCatalog, thisTree, thisvcvClade, thisvcvCladeSel ){

  listInitial()
  allOutputToNA=FALSE

   ## ===================
  ## Branch per branch analysis
  if( verbose ){
    logcat("\n ==================================\n", font="bold")
    logcat("     Branch per branch analysis\n", font="bold")
    logcat(" ==================================\n", font="bold")
  }
  branch.subTree <<- thisCatalog
  BpB.anl = na.omit(computeExpDist( elID, thisTree ))

  ## Total size of the data for that element
  n = nrow( BpB.anl )
  
  ## If no branches for that element
  if( is.null(n) || n == 0 ){
    logcat("No branches data for that element\n")
    allOutputToNA = TRUE
   # return(FALSE)
  }
  else{

    ## Reject the element if a big indel occured on a conserved branch
    if( is.logical(BpB.anl) ) return(BpB.anl)

    ## Conserved branches
    BpB.anl.sel  = na.omit(subset( BpB.anl, stat == 1 ))
    n1.branch <<- nrow( BpB.anl.sel )
    if( !is.null(n1.branch) ){
      
      ## If no branches
      if( n1.branch == 0 ){
        if( verbose ) logcat("No conserved branches for the analysis, all outputs set to NA\n")
        allOutputToNA = TRUE
      }

      ## Branches under selection
      if( in_collapseClades == "no" ){
        branch.subTree <<- na.omit(subset( thisCatalog, subTree == -1 ))
        selec.subTree = computeExpDist( elID, thisTree )
        selBranchesDist      <<- c( selBranchesDist, selec.subTree$Dist )
      }
      BpB.brWeights       <<- c( BpB.brWeights, rep( times = n1.branch, x = 1 ) )
      BpB.brStatus        <<- c( BpB.brStatus, rep( times = n1.branch, x = 1 ) )
      BpB.totBranchesDist  <<- c( BpB.totBranchesDist, BpB.anl.sel$Dist )

      ## Neutral
      BpB.anl.neut  = na.omit(subset( BpB.anl, stat == 0 ))
      n0.branch <<- nrow( BpB.anl.neut )
      ## If no neutral branches
      if( !is.null(n0.branch) ){      
        if( n0.branch == 0 ){
           if( verbose ) logcat("No neutral branches for the analysis, all outputs set to NA\n")
           allOutputToNA = TRUE        
        }          
        BpB.brWeights       <<- c( BpB.brWeights, BpB.anl.neut$w )
        BpB.brStatus        <<- c( BpB.brStatus, rep( times = n0.branch, x = 0 ) )
        BpB.totBranchesDist  <<- c( BpB.totBranchesDist, BpB.anl.neut$Dist )
        BpB.neutBranchesWeights  <<- c( BpB.neutBranchesWeights, BpB.anl.neut$w )
        BpB.neutBranchesDist      <<- c( BpB.neutBranchesDist, BpB.anl.neut$Dist )     
      }
    }
  }
  ## Perform statistical tests on the branch sample
  statTestsOnBranches( allOutputToNA )

  return(TRUE)
}


###################################################################
## Perform statistical tests on the branch sample
## p-values are recalculated to get a one-sided test
###################################################################
statTestsOnBranches = function( allOutputToNA ){

## If no data for the signal branches or species
  if( allOutputToNA ){
    n0.branch       <<- "NA"
    n1.branch       <<- "NA"
    wPearson          <<- "NA"
    wPearsonPval      <<- "NA"
    BpB.wPearson          <<- "NA"
    BpB.wPearsonPval      <<- "NA"
  } else {

     ######################################
     ## Subtree, collapsing clades
#     if( in_collapseClades != "no" ){
#       ## Weighted Pearson correlation
#       res          = wtd.cor( x=totBranchesDist, y=brStatus, weight = brWeights )
#       wPearson     <<- res[which(dimnames(res)[[2]]=="correlation")]
#       ## t-value
#       res.t = as.numeric(unlist(strsplit(summary(res)[16],":"))[2])
#       ## Degrees of freedom
#       res.n = sum(brWeights) - 2
#       ## p-value
#       ## 1-sided
#       wPearsonPval <<- pt( q = res.t, df = res.n, lower.tail = FALSE )
#       ## 2-sided
#       ##wPearsonPval <<- res[which(dimnames(res)[[2]]=="p.value")]
#     }

     ######################################
     ## Branch per branch, ignoring clades
     ## Weighted Pearson correlation
     #logcat("before run wtd.cor totBranchesDist:",BpB.totBranchesDist,"\n",file="FG.log",append = TRUE)
     #logcat("before run wtd.cor brStatus:",BpB.brStatus,"\n",file="FG.log",append = TRUE)
     #logcat("before run wtd.cor brWeights:",BpB.brWeights,"\n",file="FG.log",append = TRUE)
     Letmetry<-tryCatch(
     {
     res          = wtd.cor( x = BpB.totBranchesDist, y = BpB.brStatus, weight = BpB.brWeights )
     0
     },error=function(e) NA)
     if(!is.na(Letmetry)){
     BpB.wPearson     <<- res[which(dimnames(res)[[2]]=="correlation")]
     ## t-value
     res.t = as.numeric(unlist(strsplit(summary(res)[16],":"))[2])
     ## Degrees of freedom
     res.n = sum(BpB.brWeights) - 2
     ## p-value
     ## 1-sided
     BpB.wPearsonPval <<- pt( q = res.t, df = res.n, lower.tail = FALSE )
     ## 2-sided
     ##BpB.wPearsonPval <<- res[which(dimnames(res)[[2]]=="p.value")]
     }else{
	BpB.wPearson<<-NA
	BpB.wPearsonPval<<-NA
     }
   }

   ## Verbose
   if( verbose && !is.na(Letmetry)){
      logcat("\n")
      logcat("Results of the correlation analysis:\n", font="bold" )
      logcat("\n")
#      if( in_collapseClades != "no" ){
#         cat("    > Weighted subtree method:\n")    
#         cat(paste("     > p-value =", wPearsonPval, "\n" ))
#         cat("\n")
#      }
      logcat("    > Weighted branch method:\n")    
      logcat(paste("     > p-value =", BpB.wPearsonPval, "\n" ))
      logcat("\n") 
   }
}


##########################################################################
## Compute the distance to the expected value for a clade
##########################################################################
computeExpDist.clade = function( elID, clade, thisTree ){

  ## Data for all the branches of the clade for that element
  el.branch.subTree = subset( in_localPid, br %in% branch.subTree$br & id == elID )
  ## Keep the branches for which there is data for the element
  branch.subTree    = subset( branch.subTree, branch.subTree$br %in% el.branch.subTree$br )
  ## Merge the data for the element with the characteristics of the branches
  el.branch.subTree = merge( branch.subTree, el.branch.subTree, by="br" )
  ## distance to exp value: el(pid) - mean(pid(sim))
  el.branch.subTree$Dist = as.double(el.branch.subTree$pid) - as.double(el.branch.subTree$mPid)
  ## Subtrees: loop over the tips and sum the distances
  for( i in 1:nrow(el.branch.subTree) ){
    if( el.branch.subTree$br[i] %in% thisTree$node.label ) next
    
    ## Along the path to the tips
    Fi = 0
    Li = 0
    for( tip.anc in unlist(el.branch.subTree$path[i]) ){
      j  = which( el.branch.subTree$br == tip.anc )
      ## Sum the distances
      Fj = el.branch.subTree$Dist[j]
      if( length(Fj) == 0 ) next
      Fi = Fi + Fj
      ## Sum the branch lengths
      Lj = el.branch.subTree$len[j]
      Li = Li + Lj
    }
    ## If the total distance is smaller than -1, set it to -1
    el.branch.subTree$Dist[i] = ifelse( Fi < -1, -1, Fi )
    ## Update the length and the neutral weight
    el.branch.subTree$len[i] = Li
    el.branch.subTree$w[i] = fnBrWeights( Li )
  }
  if( verbose ){
    cat("\n")
    cat(paste("    > Clade", clade, "\n" ) )
    print( el.branch.subTree )
  }
  
  return( el.branch.subTree )
}


##########################################################################
## Compute the distance to the expected value branch per branch
##########################################################################
computeExpDist = function( elID, thisTree ){

  ## Data for all the branches for that element, omit branches with NA
  el.branch = na.omit(subset( in_localPid, br %in% branch.subTree$br & id == elID ))
  
  ## Check for duplicates
  el.branch.check = data.frame( el.branch$br, el.branch$id )
  ndupl = length(which(duplicated(el.branch.check)) == TRUE)
  if( ndupl > 0 ){
    logcat( "\nERROR: different values were found in the input data file for given branches and element\n" )
    break
  }  
  ## Keep the branches for which there is data for the element
  branch.subTree    = subset( branch.subTree, branch.subTree$br %in% el.branch$br )
  
  ## Merge the data for the element with the characteristics of the branches
  el.branch = merge( branch.subTree, el.branch, by="br" )

  ## distance: el(pid) - mean(pid(sim))
  el.branch$Dist = el.branch$pid - el.branch$mPid

  ## Threshold for big events (large indels) on conserved internal branches (not tips)
  bigEvent = subset( el.branch, stat == 1 & pid <= in_thresholdConserved  )
  bigDelEvent = FALSE
  internBranch = FALSE
  if( nrow(bigEvent) != 0 ){
    bigDelEvent = TRUE
  for( i in 1:nrow(bigEvent) )if( bigEvent$br[i] %in% thisTree$node.label ) internBranch = TRUE
  }
  
  if( bigDelEvent && internBranch ){
    logcat(paste("Large indel occured on a conserved internal branch:\n") )
    if( verbose ){
      print( bigEvent )
      logcat("\n")
    }
    return( FALSE )
  }
  
  if( verbose ){
    cat("\n")
    cat("    > Subset of branches\n",font="bold" )
    print( el.branch )
  }  
  return( el.branch )
}


##############################################################
## Remove from the covariance matrix the species without data
##############################################################
removeSpecies = function( clade, df, thisvcvClade, thisvcvCladeSel ){

  ## Neutral or selection
  if( clade < 1000 ) cladeMatrix = thisvcvClade else cladeMatrix = thisvcvCladeSel
  if( clade >= 1000 ) clade = clade - 999

  ## Species in the matrix
  colMatrix = colnames( cladeMatrix[[clade]] )
  
  ## Species from the matrix that are not in the data list
  l = which( ! colMatrix %in% df$br )
  
  ## Remove these species
  if( length(l) > 0 ) vcvClade.el <<- cladeMatrix[[clade]][ -l, -l ] else vcvClade.el <<- cladeMatrix[[clade]]
  
  ## Nodes in the data list
  l = which( df$br %in% tree$node.label )
  
  ## We average tips only
  if( length(l) > 0 ) df = df[ -l, ]
  
  ## If more than one species left, keep and order the species according to the matrix
  if( length(vcvClade.el) > 1 ){
    colMatrix = colnames( vcvClade.el )
    df = df[ match(colMatrix, df$br), ]
  } 
  return(df)
}


#####################################
## Complementary error function
#####################################
erfc = function(x){
  2 * pnorm(x * sqrt(2), lower = FALSE)
}

