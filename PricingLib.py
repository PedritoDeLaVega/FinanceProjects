from distribution import normalcum
from numpy import exp, log, sqrt, var, std
import numpy as np
from other_functions import indicatrice
from randomNumberGenerator import gaussiangenerator

class BSprice(object):
    def __init__(self, T, r, d, K, S, vol):
        self.T = T
        self.r = r
        self.d = d
        self.K = K
        self.S = S
        self.vol = vol

        self.d1 = (log(self.S / self.K) + (self.r + 0.5 * self.vol * self.vol) * self.T) \
             /(self.vol * sqrt(self.T))
        self.d2 = self.d1 - self.vol * sqrt(self.T)


    def forwardprice(self):
        return exp(-self.r*self.T)*(exp((self.r-self.d)*self.T)*self.S-self.K)

    def callprice(self):
        return self.S*normalcum(self.d1)-self.K*exp(-self.r*self.T)*normalcum(self.d2)

    def putprice(self):
        return self.K*exp(-self.r*self.T)*normalcum(-self.d2)-self.S*normalcum(-self.d1)

    def digitcall(self):
        return exp(-self.r*self.T)*normalcum(self.d2)

    def digitput(self):
        return exp(-self.r*self.T)*normalcum(-self.d2)

    def ZC(self):
        return exp(-self.r*self.T)

    def putcallparitycheck(self):
        V = BSprice.callprice(self)- \
            BSprice.putprice(self)-BSprice.forwardprice(self)
        if round(V, 4) == 0 :
            print("The put call parity is verified")
            return True
        else:
            print("The put call parity is not verified")
            return False

    def boundarycheckcall(self):
        count=0
        if BSprice.callprice(self) > self.S:
            count=count+1
            print('The higher bound of the call price is not respected')
        if BSprice.callprice(self) < (self.S-self.K*exp(-self.r*self.T)):
            count = count + 1
            print('The lower bound of the call price is not respected')
        if count == 0:
            print('The boundaries are respected for this call')
        return 0


class MCsimulations(object):
    def __init__(self, T, r, d, K, S, vol, N, seed=None):
        self.T = T
        self.r = r
        self.d = d
        self.K = K
        self.S = S
        self.vol = vol
        self.N = N
        self.W1, self.W2 = gaussiangenerator(self.N, None, seed)
        self.ST = self.S*exp((self.r - self.d)*self.T - 0.5*self.vol*self.vol*self.T +
                             self.vol*sqrt(self.T)*np.concatenate((self.W1,
                                                                   self.W2, -self.W1, -self.W2), axis=0))
        self.N = 4*self.N

    def callprice(self):
        V = []

        for i in range(self.N):
            V.append(max(self.ST[i]-self.K, 0))

        V = [exp(-self.r*self.T)*V[i]/self.N for i in range(self.N)]
        variance = var(V)
        sd = std(V)
        return sum(V), variance, sd


    def putprice(self):
        V = []

        for i in range(self.N):
            V.append(max(self.K - self.ST[i], 0))

        V = [exp(-self.r * self.T) * V[i] / self.N for i in range(self.N)]
        variance = var(V)
        sd = std(V)
        return sum(V), variance, sd

    def digitcall(self):
        V = []

        for i in range(self.N):
            V.append(indicatrice(self.ST[i]-self.K > 0))

        V = [exp(-self.r * self.T) * V[i] / self.N for i in range(self.N)]
        variance = var(V)
        sd = std(V)
        return sum(V), variance, sd

    def digitput(self):
        V = []

        for i in range(self.N):
            V.append(indicatrice(self.ST[i]-self.K < 0))

        V = [exp(-self.r * self.T) * V[i] / self.N for i in range(self.N)]
        variance = var(V)
        sd = std(V)
        return sum(V), variance, sd

class MCLEuler(object):
    def __init__(self, T, r, d, K, S, vol, N, M):
        self.T = T
        self.r = r
        self.d = d
        self.K = K
        self.S = S
        self.vol = vol
        self.N = N
        self.M = M
        self.W1, self.W2 = gaussiangenerator(self.N, self.M)
        self.N = 4*self.N
        self.ST = self.S * np.ones(self.N)
        self.dt = self.T/self.M

        for i in range(self.M):
            self.ST = self.ST * (1 + self.dt*self.r + self.vol*sqrt(self.dt) *
                                 np.concatenate((self.W1[:, i], self.W2[:, i], -self.W1[:, i], -self.W2[:, i]), axis=0))


    def callprice(self):
        V = []

        for i in range(self.N):
            V.append(max(self.ST[i]-self.K, 0))

        V = [exp(-self.r*self.T)*V[i]/self.N for i in range(self.N)]
        variance = var(V)
        sd = std(V)
        return sum(V), variance, sd

    def putprice(self):
        V = []

        for i in range(self.N):
            V.append(max(self.K - self.ST[i], 0))

        V = [exp(-self.r * self.T) * V[i] / self.N for i in range(self.N)]
        variance = var(V)
        sd = std(V)
        return sum(V), variance, sd

    def digitcall(self):
        V = []

        for i in range(self.N):
            V.append(indicatrice(self.ST[i]-self.K > 0))

        V = [exp(-self.r * self.T) * V[i] / self.N for i in range(self.N)]
        variance = var(V)
        sd = std(V)
        return sum(V), variance, sd

    def digitput(self):
        V = []

        for i in range(self.N):
            V.append(indicatrice(self.ST[i]-self.K < 0))

        V = [exp(-self.r * self.T) * V[i] / self.N for i in range(self.N)]
        variance = var(V)
        sd = std(V)
        return sum(V), variance, sd

