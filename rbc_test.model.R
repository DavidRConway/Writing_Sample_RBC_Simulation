# Generated on 2025-11-29 18:47:05 by gEcon ver. 1.2.3 (2025-04-13)
# http://gecon.r-forge.r-project.org/

# Model name: rbc_test

# info
info__ <- c("rbc_test", "/Users/dcsmac/Desktop/RBC/rbc_test.gcn", "2025-11-29 18:47:05", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("r",
                 "C",
                 "I",
                 "K_s",
                 "L_s",
                 "U",
                 "W",
                 "Y",
                 "Z")

variables_tex__ <- c("r",
                     "C",
                     "I",
                     "K^{\\mathrm{s}}",
                     "L^{\\mathrm{s}}",
                     "U",
                     "W",
                     "Y",
                     "Z")

# shocks
shocks__ <- c("epsilon_Z")

shocks_tex__ <- c("\\epsilon^{\\mathrm{Z}}")

# parameters
parameters__ <- c("alpha",
                  "beta",
                  "delta",
                  "eta",
                  "mu",
                  "psi",
                  "rho")

parameters_tex__ <- c("\\alpha",
                     "\\beta",
                     "\\delta",
                     "\\eta",
                     "\\mu",
                     "\\psi",
                     "\\rho")

# free parameters
parameters_free__ <- c("alpha",
                       "beta",
                       "delta",
                       "eta",
                       "mu",
                       "psi",
                       "rho")

# free parameters' values
parameters_free_val__ <- c(0.33,
                           0.99,
                           0.025,
                           2,
                           0.3,
                           0.8,
                           0.95)

# equations
equations__ <- c("-r[] + alpha * Z[] * K_s[-1]^(-1 + alpha) * L_s[]^(1 - alpha) = 0",
                 "-W[] + Z[] * (1 - alpha) * K_s[-1]^alpha * L_s[]^(-alpha) = 0",
                 "-Y[] + Z[] * K_s[-1]^alpha * L_s[]^(1 - alpha) = 0",
                 "-Z[] + exp(epsilon_Z[] + rho * log(Z[-1])) = 0",
                 "beta * (mu * E[][(r[1] - psi * (-delta + K_s[]^-1 * I[1])^2 + 2 * psi * K_s[]^-1 * I[1] * (-delta + K_s[]^-1 * I[1])) * C[1]^(-1 + mu) * (1 - L_s[1])^(1 - mu) * (C[1]^mu * (1 - L_s[1])^(1 - mu))^(-eta)] - mu * (1 - delta) * E[][(-1 - 2 * psi * (-delta + K_s[]^-1 * I[1])) * C[1]^(-1 + mu) * (1 - L_s[1])^(1 - mu) * (C[1]^mu * (1 - L_s[1])^(1 - mu))^(-eta)]) + mu * (-1 - 2 * psi * (-delta + K_s[-1]^-1 * I[])) * C[]^(-1 + mu) * (1 - L_s[])^(1 - mu) * (C[]^mu * (1 - L_s[])^(1 - mu))^(-eta) = 0",
                 "(-1 + mu) * C[]^mu * (1 - L_s[])^(-mu) * (C[]^mu * (1 - L_s[])^(1 - mu))^(-eta) + mu * W[] * C[]^(-1 + mu) * (1 - L_s[])^(1 - mu) * (C[]^mu * (1 - L_s[])^(1 - mu))^(-eta) = 0",
                 "I[] - K_s[] + K_s[-1] * (1 - delta) = 0",
                 "U[] - beta * E[][U[1]] - (1 - eta)^-1 * (C[]^mu * (1 - L_s[])^(1 - mu))^(1 - eta) = 0",
                 "-C[] - I[] + Y[] - psi * K_s[-1] * (-delta + K_s[-1]^-1 * I[])^2 = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3,
                                 3, 3, 4, 5, 5, 5, 5, 5, 6, 6,
                                 6, 7, 7, 8, 8, 8, 9, 9, 9, 9),
                           j = c(1, 4, 5, 9, 4, 5, 7, 9, 4, 5,
                                 8, 9, 9, 1, 2, 3, 4, 5, 2, 5,
                                 7, 3, 4, 2, 5, 6, 2, 3, 4, 8),
                           x = c(2, 1, 2, 2, 1, 2, 2, 2, 1, 2,
                                 2, 2, 3, 4, 6, 6, 3, 6, 2, 2,
                                 2, 2, 3, 2, 2, 6, 2, 2, 1, 2),
                           dims = c(9, 9))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 9))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(9, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(1, 2, 3, 4, 5, 5, 5, 5, 5, 6,
                                     6, 7, 8, 8, 8, 9, 9),
                               j = c(1, 1, 1, 7, 2, 3, 4, 5, 6, 4,
                                     5, 3, 2, 4, 5, 3, 6),
                               x = rep(1, 17), dims = c(9, 7))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 7))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(4),
                             j = c(1),
                             x = rep(1, 1), dims = c(9, 1))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(9)
    r[1] = -v[1] + pf[1] * v[9] * v[4]^(-1 + pf[1]) * v[5]^(1 - pf[1])
    r[2] = -v[7] + v[9] * (1 - pf[1]) * v[4]^pf[1] * v[5]^(-pf[1])
    r[3] = -v[8] + v[9] * v[4]^pf[1] * v[5]^(1 - pf[1])
    r[4] = -v[9] + exp(pf[7] * log(v[9]))
    r[5] = pf[2] * (pf[5] * (v[1] - pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[5] * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * (1 - pf[3]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])) + pf[5] * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    r[6] = (-1 + pf[5]) * v[2]^pf[5] * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) + pf[5] * v[7] * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    r[7] = v[3] - v[4] + v[4] * (1 - pf[3])
    r[8] = v[6] - pf[2] * v[6] - (1 - pf[4])^-1 * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(1 - pf[4])
    r[9] = -v[2] - v[3] + v[8] - pf[6] * v[4] * (-pf[3] + v[3] * v[4]^-1)^2

    return(r)
}

# calibrating equations
calibr_eq__ <- function(v, pc, pf)
{
    r <- numeric(0)

    return(r)
}

# steady state and calibrating equations Jacobian
ss_calibr_eq_jacob__ <- function(v, pc, pf)
{
    r <- numeric(0)
    jac <- numeric(30)
    jac[1] = -1
    jac[2] = pf[1] * v[9] * (-1 + pf[1]) * v[4]^(-2 + pf[1]) * v[5]^(1 - pf[1])
    jac[3] = pf[1] * v[9] * (1 - pf[1]) * v[4]^(-1 + pf[1]) * v[5]^(-pf[1])
    jac[4] = pf[1] * v[4]^(-1 + pf[1]) * v[5]^(1 - pf[1])
    jac[5] = pf[1] * v[9] * (1 - pf[1]) * v[4]^(-1 + pf[1]) * v[5]^(-pf[1])
    jac[6] = -pf[1] * v[9] * (1 - pf[1]) * v[4]^pf[1] * v[5]^(-1 - pf[1])
    jac[7] = -1
    jac[8] = (1 - pf[1]) * v[4]^pf[1] * v[5]^(-pf[1])
    jac[9] = pf[1] * v[9] * v[4]^(-1 + pf[1]) * v[5]^(1 - pf[1])
    jac[10] = v[9] * (1 - pf[1]) * v[4]^pf[1] * v[5]^(-pf[1])
    jac[11] = -1
    jac[12] = v[4]^pf[1] * v[5]^(1 - pf[1])
    jac[13] = -1 + pf[7] * v[9]^-1 * exp(pf[7] * log(v[9]))
    jac[14] = pf[2] * pf[5] * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    jac[15] = pf[2] * (-pf[4] * pf[5]^2 * (v[1] - pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)) * (v[2]^(-1 + pf[5]))^2 * ((1 - v[5])^(1 - pf[5]))^2 * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4]) + pf[5] * (-1 + pf[5]) * (v[1] - pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-2 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) + pf[4] * pf[5]^2 * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * (1 - pf[3]) * (v[2]^(-1 + pf[5]))^2 * ((1 - v[5])^(1 - pf[5]))^2 * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4]) - pf[5] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * (1 - pf[3]) * v[2]^(-2 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])) - pf[4] * pf[5]^2 * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * (v[2]^(-1 + pf[5]))^2 * ((1 - v[5])^(1 - pf[5]))^2 * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4]) + pf[5] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-2 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    jac[16] = pf[2] * (2 * pf[5] * pf[6] * v[3] * v[4]^-2 * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) + 2 * pf[5] * pf[6] * v[4]^-1 * (1 - pf[3]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])) - 2 * pf[5] * pf[6] * v[4]^-1 * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    jac[17] = pf[2] * (-2 * pf[5] * pf[6] * v[3]^2 * v[4]^-3 * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - 2 * pf[5] * pf[6] * v[3] * v[4]^-2 * (1 - pf[3]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])) + 2 * pf[5] * pf[6] * v[3] * v[4]^-2 * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    jac[18] = pf[2] * (pf[5] * (-1 + pf[5]) * (v[1] - pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[5] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * (1 - pf[3]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5] * (-1 + pf[5]) * (v[1] - pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(-pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4]) + pf[4] * pf[5] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * (1 - pf[3]) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(-pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])) + pf[5] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(-pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])
    jac[19] = pf[5] * (-1 + pf[5]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5]^2 * v[7] * (v[2]^(-1 + pf[5]))^2 * ((1 - v[5])^(1 - pf[5]))^2 * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4]) + pf[5] * v[7] * (-1 + pf[5]) * v[2]^(-2 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5] * (-1 + pf[5]) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(-pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])
    jac[20] = -pf[4] * (-1 + pf[5])^2 * (v[2]^pf[5])^2 * ((1 - v[5])^(-pf[5]))^2 * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4]) + pf[5] * (-1 + pf[5]) * v[2]^pf[5] * (1 - v[5])^(-1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) + pf[5] * v[7] * (-1 + pf[5]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5] * v[7] * (-1 + pf[5]) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(-pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])
    jac[21] = pf[5] * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    jac[22] = 1
    jac[23] = -pf[3]
    jac[24] = -pf[5] * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    jac[25] = -(-1 + pf[5]) * v[2]^pf[5] * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    jac[26] = 1 - pf[2]
    jac[27] = -1
    jac[28] = -1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)
    jac[29] = -pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)
    jac[30] = 1
    jacob <- sparseMatrix(i = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3,
                                3, 3, 4, 5, 5, 5, 5, 5, 6, 6,
                                6, 7, 7, 8, 8, 8, 9, 9, 9, 9),
                          j = c(1, 4, 5, 9, 4, 5, 7, 9, 4, 5,
                                8, 9, 9, 1, 2, 3, 4, 5, 2, 5,
                                7, 3, 4, 2, 5, 6, 2, 3, 4, 8),
                          x = jac, dims = c(9, 9))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(7)
    Atm1x[1] = pf[1] * v[9] * (-1 + pf[1]) * v[4]^(-2 + pf[1]) * v[5]^(1 - pf[1])
    Atm1x[2] = pf[1] * v[9] * (1 - pf[1]) * v[4]^(-1 + pf[1]) * v[5]^(-pf[1])
    Atm1x[3] = pf[1] * v[9] * v[4]^(-1 + pf[1]) * v[5]^(1 - pf[1])
    Atm1x[4] = pf[7] * v[9]^-1 * exp(pf[7] * log(v[9]))
    Atm1x[5] = 2 * pf[5] * pf[6] * v[3] * v[4]^-2 * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    Atm1x[6] = 1 - pf[3]
    Atm1x[7] = -pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)
    Atm1 <- sparseMatrix(i = c(1, 2, 3, 4, 5, 7, 9),
                         j = c(4, 4, 4, 9, 4, 4, 4),
                         x = Atm1x, dims = c(9, 9))

    Atx <- numeric(25)
    Atx[1] = -1
    Atx[2] = pf[1] * v[9] * (1 - pf[1]) * v[4]^(-1 + pf[1]) * v[5]^(-pf[1])
    Atx[3] = pf[1] * v[4]^(-1 + pf[1]) * v[5]^(1 - pf[1])
    Atx[4] = -pf[1] * v[9] * (1 - pf[1]) * v[4]^pf[1] * v[5]^(-1 - pf[1])
    Atx[5] = -1
    Atx[6] = (1 - pf[1]) * v[4]^pf[1] * v[5]^(-pf[1])
    Atx[7] = v[9] * (1 - pf[1]) * v[4]^pf[1] * v[5]^(-pf[1])
    Atx[8] = -1
    Atx[9] = v[4]^pf[1] * v[5]^(1 - pf[1])
    Atx[10] = -1
    Atx[11] = pf[5] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-2 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5]^2 * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-2 + 2 * pf[5]) * (1 - v[5])^(2 - 2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])
    Atx[12] = -2 * pf[5] * pf[6] * v[4]^-1 * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    Atx[13] = pf[2] * (-2 * pf[5] * pf[6] * v[3]^2 * v[4]^-3 * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - 2 * pf[5] * pf[6] * v[3] * v[4]^-2 * (1 - pf[3]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]))
    Atx[14] = pf[5] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(1 - 2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])
    Atx[15] = pf[5] * (-1 + pf[5]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5] * (-1 + pf[5]) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(1 - 2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4]) - pf[4] * pf[5]^2 * v[7] * v[2]^(-2 + 2 * pf[5]) * (1 - v[5])^(2 - 2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4]) + pf[5] * v[7] * (-1 + pf[5]) * v[2]^(-2 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    Atx[16] = -pf[4] * (-1 + pf[5])^2 * v[2]^(2 * pf[5]) * (1 - v[5])^(-2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4]) + pf[5] * (-1 + pf[5]) * v[2]^pf[5] * (1 - v[5])^(-1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) + pf[5] * v[7] * (-1 + pf[5]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5] * v[7] * (-1 + pf[5]) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(1 - 2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])
    Atx[17] = pf[5] * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    Atx[18] = 1
    Atx[19] = -1
    Atx[20] = -pf[5] * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    Atx[21] = (1 - pf[5]) * v[2]^pf[5] * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    Atx[22] = 1
    Atx[23] = -1
    Atx[24] = -1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)
    Atx[25] = 1
    At <- sparseMatrix(i = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4,
                             5, 5, 5, 5, 6, 6, 6, 7, 7, 8,
                             8, 8, 9, 9, 9),
                       j = c(1, 5, 9, 5, 7, 9, 5, 8, 9, 9,
                             2, 3, 4, 5, 2, 5, 7, 3, 4, 2,
                             5, 6, 2, 3, 8),
                         x = Atx, dims = c(9, 9))

    Atp1x <- numeric(5)
    Atp1x[1] = pf[2] * pf[5] * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4])
    Atp1x[2] = pf[2] * (pf[5] * ((-1 + pf[5]) * (v[1] - pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-2 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5] * (v[1] - pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-2 + 2 * pf[5]) * (1 - v[5])^(2 - 2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])) - pf[5] * (1 - pf[3]) * ((-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-2 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * pf[5] * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-2 + 2 * pf[5]) * (1 - v[5])^(2 - 2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])))
    Atp1x[3] = pf[2] * (2 * pf[5] * pf[6] * v[3] * v[4]^-2 * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) + 2 * pf[5] * pf[6] * v[4]^-1 * (1 - pf[3]) * v[2]^(-1 + pf[5]) * (1 - v[5])^(1 - pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]))
    Atp1x[4] = pf[2] * (pf[5] * ((-1 + pf[5]) * (v[1] - pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * (-1 + pf[5]) * (v[1] - pf[6] * (-pf[3] + v[3] * v[4]^-1)^2 + 2 * pf[6] * v[3] * v[4]^-1 * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(1 - 2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])) - pf[5] * (1 - pf[3]) * ((-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + pf[5]) * (1 - v[5])^(-pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-pf[4]) - pf[4] * (-1 + pf[5]) * (-1 - 2 * pf[6] * (-pf[3] + v[3] * v[4]^-1)) * v[2]^(-1 + 2 * pf[5]) * (1 - v[5])^(1 - 2 * pf[5]) * (v[2]^pf[5] * (1 - v[5])^(1 - pf[5]))^(-1 - pf[4])))
    Atp1x[5] = -pf[2]
    Atp1 <- sparseMatrix(i = c(5, 5, 5, 5, 8),
                         j = c(1, 2, 3, 5, 6),
                         x = Atp1x, dims = c(9, 9))

    Aepsx <- numeric(1)
    Aepsx[1] = exp(pf[7] * log(v[9]))
    Aeps <- sparseMatrix(i = c(4),
                         j = c(1),
                         x = Aepsx, dims = c(9, 1))

    return(list(Atm1, At, Atp1, Aeps))
}

ext__ <- list()

# create model object
gecon_model(model_info = info__,
            index_sets = index_sets__,
            variables = variables__,
            variables_tex = variables_tex__,
            shocks = shocks__,
            shocks_tex = shocks_tex__,
            parameters = parameters__,
            parameters_tex = parameters_tex__,
            parameters_free = parameters_free__,
            parameters_free_val = parameters_free_val__,
            equations = equations__,
            calibr_equations = calibr_equations__,
            var_eq_map = vareqmap__,
            shock_eq_map = shockeqmap__,
            var_ceq_map = varcalibreqmap__,
            cpar_eq_map = calibrpareqmap__,
            cpar_ceq_map = calibrparcalibreqmap__,
            fpar_eq_map = freepareqmap__,
            fpar_ceq_map = freeparcalibreqmap__,
            ss_function = ss_eq__,
            calibr_function = calibr_eq__,
            ss_calibr_jac_function = ss_calibr_eq_jacob__,
            pert = pert1__,
            ext = ext__)
