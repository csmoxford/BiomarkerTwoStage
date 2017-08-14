library(binom)


n = 10000

x = 1:999 * 10

res = binom.bayes(x, n,
            conf.level = 0.95,
            type = c("central"),
            prior.shape1 = 0.5,
            prior.shape2 = 0.5,
            tol = .Machine$double.eps^0.5,
            maxit = 1000)


plot(x/n,x/n, type = "l", xaxs = "i", yaxs = "i")

lines(x/n,res$lower)
lines(x/n,res$upper)

plot(x/n, res$upper - res$lower, type = "l")
