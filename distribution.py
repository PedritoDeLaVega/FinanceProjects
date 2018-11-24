from math import log, sqrt, pi, exp

def norminv(x):

    moronumbers = [2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637,
                   -8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833,
                   0.3374754822726147, 0.9761690190917186, 0.1607979714918209,
                   0.0276438810333863, 0.038405729373609, 0.0003951896511919,
                   0.0000321767881768, 0.000000288167364, 0.0000003960315187]
    y = x-0.5
    a = 0
    b = 0
    c = 0
    if abs(y) < 0.42:
        r = y*y
        for i in range(4):
            a = a + moronumbers[i]*pow(r,i)
            b = b + moronumbers[i+5]*pow(r,i+1)
        return y*a/(b+1)
    else:
        if y < 0:
            r = x
        else:
            r = 1-x
        s = log(-log(r))
        for i in range(9):
            c = c+moronumbers[i+7]*pow(s,i)
        return c


def normalcum(x):
    if x < 0:
        x = -x
        k = 1/(1+0.2316419*x)
        return (1/sqrt(2*pi))*exp(-0.5*pow(x, 2))*k*(0.319381530+k*(-0.356563782 +
            k * (1.781477937+k*(-1.821255978+1.330274429*k))))
    else:
        k = 1/(1+0.2316419*x)
        return 1 - (1/sqrt(2*pi))*exp(-0.5*pow(x, 2))*k*(0.319381530+k*(-0.356563782 +
            k * (1.781477937+k*(-1.821255978+1.330274429*k))))

def normaldens(x):
    return exp(-0.5*x*x)/sqrt(2*pi)
