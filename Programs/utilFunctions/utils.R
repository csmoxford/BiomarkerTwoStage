
kableTable = function(what, dir="h"){

  if(class(what)[1] != "table"){
    what = table(what)
  }

  if(dir=="h"){
    df = data.frame(1)
  for(i in 1:length(what)){
    df[,names(what)[i]] = what[i]
  }
  df$X1 = NULL
  } else {
    df = data.frame(Name=rep("",0), Count=rep(0,0),stringsAsFactors = FALSE)
    for(i in 1:length(what)){
      df[i,] = c(names(what)[i],what[i])
    }
  }



  return(kable(df))
}


getDetailedSummary = function(file,i){


  cat("\n\n### Example",i)

  load(file)

  stopReasons = sapply(res@sims,function(x) x$decision@stopReason)
  nrecruited = sapply(res@sims, function(x) return(sum(!is.na(x$data$outcome))))
  nseen = sapply(res@sims, function(x) return(dim(x$data)[1]))

  print(kable(data.frame(Biomarker=c("No change / normal","Partial loss","Complete loss"),"Simulated Truth" = res@p$biomarkerResp,stringsAsFactors = FALSE)))

  cat("\n\n### Stopping reasons")
  ta = table(stopReasons)/res@nSim
  print(kableTable(ta,"v"))

  cat("\n\n### Stages trail")
  trail = sapply(res@sims, function(x) return(paste(x$decision@trail,collapse = ",")))
  ta = table(trail)/res@nSim
  print(kableTable(ta,"v"))
}


getSummary = function(res){

  coln = names(res@sims[[1]]$data)

  stopReasons = sapply(res@sims,function(x) x$decision@stopReason)

    nrecruited = sapply(res@sims, function(x) return(sum(!is.na(x$data$outcome))))

  if("tested" %in% coln){
    nseen = sapply(res@sims, function(x) return(sum(x$data$tested == 1, na.rm = TRUE)))
  } else {
    nseen = sapply(res@sims, function(x) return(dim(x$data)[1]))
  }
  print(kable(data.frame(Biomarker=paste("Level", 1:length(res@p$biomarkerResp)-1),"Simulated Truth" = res@p$biomarkerResp,stringsAsFactors = FALSE)))

  cat("\n\n### Stopping reasons")
  stopTab = table(stopReasons)/res@nSim
  print(kableTable(stopTab,"v"))

  cat("\n\n### Stages trail")
  trail = sapply(res@sims, function(x) return(paste(x$decision@trail,collapse = ",")))
  trailTab = table(trail)/res@nSim
  print(kableTable(trailTab,"v"))

  cat("\n\n### Final decision")
  recommend = rep("futile",length(stopReasons))
  recommend[grepl("success 012",stopReasons)] = "group 012"
  recommend[grepl("success 12",stopReasons)] = "group 12"
  recommend[grepl("success 2",stopReasons)] = "group 2"
  recommend[grepl("success 01",stopReasons)] = "group 01"
  recommend[grepl("success 1",stopReasons)] = "group 1"

  print(kableTable(recommend,"v"))

  tra = unlist(sapply(res@sims, function(x) return(x$decision@trail)))
  ta = table(tra)/res@nSim

  par(mfrow=c(2,2),mar=c(3,3,2,2))
  barplot(trailTab,yaxs="i",ylim=c(0,1), main="Stages taken")
  hist(nrecruited,xaxs="i",main="Number recruited")
  hist(nseen,yaxs="i",xaxs="i",main="Number screened")
  barplot(stopTab,main="Stop reason",ylim=c(0,1))
}

inv_logit = function(alpha) 1/(1+exp(-alpha))

postProbDist = function(model, level){
  m = model$BUGSoutput$sims.list
  if(level == 0){
    return(inv_logit(m$alpha))
  } else if(level == 1){
    return(inv_logit(m$alpha + m$beta))
  } else if(level == 2){
    return(inv_logit(m$alpha + m$beta + m$delta))
  } else {
    stop("What level?")
  }
}

postProb = function(prior, threshold, side = 0){
  pr = pbeta(threshold, prior[1], prior[2])
  if(side){
    pr = 1 - pr
  }
  return(pr)
}
