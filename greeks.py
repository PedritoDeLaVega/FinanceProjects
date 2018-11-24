from distribution import normalcum, normaldens
from numpy import exp, log, sqrt
from PricingLib import BSprice, MCsimulations
from randomNumberGenerator import gaussiangenerator
import numpy as np
from other_functions import indicatrice

class GreekBS(object):
    def __init__(self, T, r, d, K, S, vol):
        self.T = T
        self.r = r
        self.d = d
        self.K = K
        self.S = S
        self.vol = vol

        self.d1 = (log(self.S / self.K) + (self.r - self.d + 0.5 * self.vol * self.vol) * self.T) \
             /(self.vol * sqrt(self.T))
        self.d2 = self.d1 - self.vol * sqrt(self.T)

    def delta(self):
        print("Calculating a call option delta with Black-Scholes. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return exp(-self.d*self.T)*normalcum(self.d1)

    def gamma(self):
        print("Calculating a call option gamma with Black-Scholes. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return self.K*exp(-self.r*self.T)*normaldens(self.d2)/(self.S*self.S*self.vol*sqrt(self.T))

    def vega(self):
        print("Calculating a call option vega with Black-Scholes.S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return exp(-self.d*self.T)*normaldens(self.d1)*self.S*sqrt(self.T)

    def rho(self):
        print("Calculating a call option rho with Black-Scholes. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return self.K*self.T*exp(-self.r*self.T)*normalcum(self.d2)

    def theta(self):
        print("Calculating a call option theta with Black-Scholes. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return -exp(-self.d*self.T)*self.S*self.vol*normaldens(self.d1)/(2*sqrt(self.T)) - \
               self.r*self.K*exp(-self.r*self.T)*normalcum(self.d2) + \
               self.d*self.S*exp(-self.d*self.T)*normalcum(self.d1)

class GreekFD(object):
    def __init__(self, T, r, d, K, S, vol, eps):
        self.T = T
        self.r = r
        self.d = d
        self.K = K
        self.S = S
        self.vol = vol
        self.eps = eps

    def delta(self):
        print("Calculating a call option delta with Finite Difference. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/self.eps)*(BSprice(self.T, self.r, self.d, self.K, self.S + self.eps, self.vol).callprice() -
                          BSprice(self.T, self.r, self.d, self.K, self.S, self.vol).callprice())

    def gamma(self):
        print("Calculating a call option gamma with Finite Difference. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/(self.eps*self.eps)) * \
               (BSprice(self.T, self.r, self.d, self.K, self.S + self.eps, self.vol).callprice() -
                2 * BSprice(self.T, self.r, self.d, self.K, self.S, self.vol).callprice() +
                BSprice(self.T, self.r, self.d, self.K, self.S - self.eps, self.vol).callprice())

    def vega(self):
        print("Calculating a call option vega with Finite Difference. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/self.eps)*(BSprice(self.T, self.r, self.d, self.K, self.S, self.vol + self.eps).callprice() -
                          BSprice(self.T, self.r, self.d, self.K, self.S, self.vol).callprice())

    def rho(self):
        print("Calculating a call option rho with Finite Difference. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/self.eps)*(BSprice(self.T, self.r + self.eps, self.d, self.K, self.S, self.vol).callprice() -
                          BSprice(self.T, self.r, self.d, self.K, self.S, self.vol).callprice())

    def theta(self):
        print("Calculating a call option theta with Finite Difference. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/self.eps)*(BSprice(self.T + self.eps, self.r, self.d, self.K, self.S, self.vol).callprice() -
                          BSprice(self.T, self.r, self.d, self.K, self.S, self.vol).callprice())

class GreekMCL(object):
    def __init__(self, T, r, d, K, S, vol, N, eps, seed):
        self.T = T
        self.r = r
        self.d = d
        self.K = K
        self.S = S
        self.vol = vol
        self.N = N
        self.eps = eps
        self.seed = seed

    def delta(self):
        print("Calculating a call option delta with Monte Carlo. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/self.eps) * \
               (MCsimulations(self.T, self.r, self.d, self.K, self.S + self.eps, self.vol, self.N, self.seed).callprice()[0] -
                MCsimulations(self.T, self.r, self.d, self.K, self.S, self.vol, self.N, self.seed).callprice()[0])

    def gamma(self):
        print("Calculating a call option gamma with Monte Carlo. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/(self.eps*self.eps)) * \
               (MCsimulations(self.T, self.r, self.d, self.K, self.S + self.eps, self.vol, self.N, self.seed).callprice()[0] -
               2 * MCsimulations(self.T, self.r, self.d, self.K, self.S, self.vol, self.N, self.seed).callprice()[0] +
               MCsimulations(self.T, self.r, self.d, self.K, self.S - self.eps, self.vol, self.N, self.seed).callprice()[0])
        # return (1/self.eps) * \
        #        (GreekMCL(self.T, self.r, self.d, self.K, self.S, self.vol, self.N, self.eps, self.seed).delta() -
        #        GreekMCL(self.T, self.r, self.d, self.K, self.S - self.eps, self.vol, self.N, self.eps, self.seed).delta())

    def vega(self):
        print("Calculating a call option vega with Monte Carlo. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/self.eps) * \
               (MCsimulations(self.T, self.r, self.d, self.K, self.S, self.vol + self.eps, self.N, self.seed).callprice()[0] -
                MCsimulations(self.T, self.r, self.d, self.K, self.S, self.vol, self.N, self.seed).callprice()[0])

    def rho(self):
        print("Calculating a call option rho with Monte Carlo. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/self.eps) * \
               (MCsimulations(self.T, self.r + self.eps, self.d, self.K, self.S, self.vol, self.N, self.seed).callprice()[0] -
                MCsimulations(self.T, self.r, self.d, self.K, self.S, self.vol, self.N, self.seed).callprice()[0])

    def theta(self):
        print("Calculating a call option theta with Monte Carlo. S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        return (1/self.eps) * \
               (MCsimulations(self.T + self.eps, self.r, self.d, self.K, self.S, self.vol, self.N, self.seed).callprice()[0] -
                MCsimulations(self.T, self.r, self.d, self.K, self.S, self.vol, self.N, self.seed).callprice()[0])



class GreekMCLlikelihoodRatio(object):
    def __init__(self, T, r, d, K, S, vol, N):
        self.T = T
        self.r = r
        self.d = d
        self.K = K
        self.S = S
        self.vol = vol
        self.N = N
        self.W1, self.W2 = gaussiangenerator(self.N, None, None)
        self.N = 4 * self.N
        self.ST = self.S * exp((self.r - self.d) * self.T - 0.5 * self.vol * self.vol * self.T +
                               self.vol * sqrt(self.T) * np.concatenate((self.W1,
                                                                         self.W2, -self.W1, -self.W2), axis=0))
        self.Z = (log(self.ST/self.S)-(self.r - 0.5*self.vol*self.vol)*self.T)/(self.vol*sqrt(self.T))


    def delta(self):
        V = []
        print("Calculating a call option delta with Likelihood ratio.S = " + str(self.S) + " T = " + str(self.T) +
              " r = " + str(self.r) + " d = " + str(self.d) + " K = " + str(self.K) + " vol = " + str(self.vol))
        for i in range(self.N):
            V.append(max(self.ST[i] - self.K, 0) * self.Z[i] / (self.S * self.vol * sqrt(self.T)))

        return exp(-self.r*self.T) * sum(V) / self.N

