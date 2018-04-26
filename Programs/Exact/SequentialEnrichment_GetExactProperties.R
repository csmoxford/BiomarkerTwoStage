


enrichGetPropertiesTable = function(design, designName,pr.positive = 0.5){
  n1.p = design$"n1.positve"
  r1.p = design$"r1.positive"

  n.p = design$"n.positive"
  r.p = design$"r.positive"


  n1.n = design$"n1.negative"
  r1.n = design$"r1.negative"

  n.n = design$"n.negative"
  r.n = design$"r.negative"


  getProperties = function(n1,r1,n,r,p){

    prob.s1 = dbinom(0:n1,n1,p)

    stop.s1 = sum(prob.s1[0:n1 <= r1])

    prob.s1[0:n1 <= r1] = 0

    total = 0
    prob.s2 = rep(0, n+1)
    for(i in 0:n1){
      prob.s2[i + 1 + 0:(n - n1)] = prob.s2[i + 1 + 0:(n - n1)] + prob.s1[i+1] * dbinom(0:(n-n1),n-n1,p)
    }

    stop.s2 = sum(prob.s2[0:n <= r])

    stop.t1 = stop.s1 + stop.s2

    return(list(
      EN = n1 + (n-n1) * (1-stop.s1),
      prStopS1 = stop.s1,
      prEFF = 1 - stop.t1
    ))

  }



  dta = data.frame(matrix(0,0,13))
  names(dta) = c("Design","P(+)","P(b-)","P(b+)","interim stop","interim unselected","interim positive","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested")

  p.p = p0
  p.n = p0
  t1 = getProperties(n1.p,r1.p,n.p,r.p,p.p)
  t2 = getProperties(n1.n,r1.n,n.n,r.n,p.n)

  dta[1,] = c(designName,pr.positive,p.n,p.p,NA,NA,NA,roundWZero(c((1 - t1$prEFF),t1$prEFF*t2$prEFF,t1$prEFF * (1 - t2$prEFF)),3),roundWZero(t1$EN + t1$prEFF*t2$EN,1),ifelse(t1$prEFF*(1-t2$prStopS1) > 0.01, n.n + n.p, n.p + n1.n), roundWZero(t1$EN/pr.positive +t1$prEFF*t2$EN/ (1-pr.positive),1))

  p.p = p1
  p.n = p0
  t1 = getProperties(n1.p,r1.p,n.p,r.p,p.p)
  t2 = getProperties(n1.n,r1.n,n.n,r.n,p.n)
  dta[2,] = c(designName, pr.positive,p.n,p.p,NA,NA,NA,roundWZero(c((1 - t1$prEFF),t1$prEFF*t2$prEFF,t1$prEFF * (1 - t2$prEFF)),3),roundWZero(t1$EN + t1$prEFF*t2$EN,1),ifelse(t1$prEFF*(1-t2$prStopS1) > 0.01, n.n + n.p, n.p + n1.n), roundWZero(t1$EN/pr.positive +t1$prEFF*t2$EN/ (1-pr.positive),1))


  p.p = p1
  p.n = p1
  t1 = getProperties(n1.p,r1.p,n.p,r.p,p.p)
  t2 = getProperties(n1.n,r1.n,n.n,r.n,p.n)
  dta[3,] = c(designName, pr.positive,p.n,p.p,NA,NA,NA,roundWZero(c((1 - t1$prEFF),t1$prEFF*t2$prEFF,t1$prEFF * (1 - t2$prEFF)),3),roundWZero(t1$EN + t1$prEFF*t2$EN,1),ifelse(t1$prEFF*(1-t2$prStopS1) > 0.01, n.n + n.p, n.p + n1.n), roundWZero(t1$EN/pr.positive +t1$prEFF*t2$EN/ (1-pr.positive),1))

  dta
}


