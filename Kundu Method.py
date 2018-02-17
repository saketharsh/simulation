
# This Python 3 environment comes with many helpful analytics libraries installed
# It is defined by the kaggle/python docker image:
# https://github.com/kaggle/docker-python

import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import matplotlib.pyplot as plt
import math
import time

np.random.seed(237)  # To reproduce the result later


# Generation from Gamma(alpha, beta) when alpha is integer
def generate(alpha, beta, size):  # size is the number of sample
    sample = np.zeros(size)
    for i in range(size):
        U = np.random.uniform(0, 1, alpha)
        lnU = np.log(U)
        sample[i] = (-1 * np.sum(lnU)) / beta
    return sample

import scipy.stats as stat

#Kundu Method Ravi Waala 
for f in [1.5, 3.2, 100.7]:
    integral = int(f)
    alpha = f - integral
    start = time.time()
    rejection = np.zeros(1000)
    samples_int = generate(integral, 1, 1000)
    samples = np.zeros(1000)
    d = 1.0334 - 0.0766 * math.exp(2.2942 * alpha)
    a = math.pow(2, alpha) * math.pow((1 - math.exp(-1 * d / 2)), alpha)
    b = alpha * math.pow(d, alpha - 1) * math.exp(-1 * d)
    c = a + b
    for i in range(1000):
        reject = 0
        while True:
            U = np.random.uniform(0, 1)
            if(U <= a/(a+b)):
                X = -2 * np.log(1 - (math.pow(c*U, 1/alpha)/2))
            else:
                X = -1 * np.log((c * (1 - U))/(alpha * math.pow(d, alpha - 1)))
            V = np.random.uniform(0, 1)
            if(X<=d):
                if(V <= (math.pow(X, alpha - 1) * math.exp(-1 * X / 2))/(math.pow(2, alpha - 1) * math.pow((1 - math.exp(-1 * X / 2)), alpha - 1))):
                    samples[i] = X
                    break
                else:
                    reject += 1
            else:
                if(V <= math.pow(d/X, 1-alpha)):
                    samples[i] = X
                    break
                else:
                    reject += 1
        rejection[i] = reject
    end = time.time()
    if (f == 100.7):
        x = np.linspace(60, 140, 100)
    else:
        x = np.linspace (0, 10, 100)     
    #print(integral + alpha)
    samples = samples + samples_int
    y1 = stat.gamma.pdf(x, a = f, scale=1)
    #print(integral + alpha)
    plt.plot(x, y1, "y-")
    plt.hist(samples, bins = 20, normed = True)
    plt.show()
    print(np.sum(rejection)/(1000 + np.sum(rejection)))
    print(end-start)