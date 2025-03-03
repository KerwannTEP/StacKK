# include("Args.jl")

const G = 1.0
const M = 1.0

const alpha = -a^2
const gamma = -c^2

const Delta2 = abs(gamma-alpha)
const Delta = sqrt(Delta2)
const sqrtDelta = sqrt(Delta)

const ta = a/(a+c)
const tc = c/(a+c)

const h = (a^2-c^2)/a^2 # Squared eccentricity

const _E0 = -G*M/(a+c)

const Jumax = 5.0
const Jvmax = 5.0
const Lzmax = 5.0

const cutoffEff = 0.01 # Effective anomaly cutoff precision

const err_shell = 10^(-10)
const err_u = 0.00001

const INV_PI = 1.0/pi
const PI = pi


