from numpy import random, sqrt, log, sin, cos, pi

def gaussiangenerator(N, M = None, s = None):
    if s == None:
        random.seed()
    else:
        random.seed(s)

    if M == None:
        u1 = random.rand(N)
        u2 = random.rand(N)
    else:
        u1 = random.rand(N, M)
        u2 = random.rand(N, M)

    return sqrt(-2*log(u1))*cos(2*pi*u2), sqrt(-2*log(u1))*sin(2*pi*u2)

