# For example tests, see: 
# https://github.com/WSDeWitt/BDMS/blob/v0.5.0/tests/test_bdms.py
test_that("BirthDeath() reproduces @Celentano2024", {
  lambda <- rbind(fit = c(fit = 1, unfit = 1), unfit = c(0.25, 0.25))
  mu <- c(fit = 0.25, unfit = 0.25)
  gamma <- cbind(fit = c(fit = 0, unfit = 0.5),
                unfit = c(fit = 0.5, unfit = 0))

  BirthDeath(
   pi = c(fit = 0.5, unfit = 0.5),
   lambda = lambda,
   mu = mu,
   psi = c(fit = 0, unfit = 0),
   gamma = gamma,
   tMax = 500,
   nMax = 500
   )

})
