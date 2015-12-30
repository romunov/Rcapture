# simulate capture history
library(edgeEffect)
my.ch <- edgeEffect::edgeEffect(a = 1, b = 1, r = 10, arena = 10, world = 40, N = 100, prob = 0.2, J = 5, space.use = "halfnormal",
                    sigma = 1, two.groups = FALSE, seed = 357, filename = "temp", verbose = FALSE, out.ch = TRUE)

library(Rcapture)

descriptive(my.ch)
plot(descriptive(my.ch))

(m1 <- closedp.0(my.ch, dtype = "hist"))

rbindCIs(list(profileCI(my.ch, m = "M0"), 
              profileCI(my.ch, m = "Mh", h = "Poisson", a = 2),
              profileCI(my.ch, m = "Mh", h = "Chao"),
              profileCI(my.ch, m = "Mh", h = "Darroch")
))

m2 <- closedp.t(my.ch, dtype = "hist")
uifit(m2)


