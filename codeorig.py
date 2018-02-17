# This Python 3 environment comes with many helpful analytics libraries installed
# It is defined by the kaggle/python docker image: https://github.com/kaggle/docker-python

import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import matplotlib.pyplot as plt

# Generation from Gamma(alpha, beta) when alpha is integer
def generate(alpha, beta, size): # size is the number of sample
    sample = np.zeros(size)
    for i in range(size):
        U = np.random.uniform(0, 1, alpha)
        lnU = np.log(U)
        sample[i] = (-1 * np.sum(lnU))/beta
    return sample

import scipy.stats as stat

# The shape and scale parameters for Gamma distribution
b = 1
generated = np.empty(1000) # The array to store the output
np.random.seed(237) # To reproduce the result later
# Generating for integral part
time_list = []
import time
for a_1 in [1.5, 3.2, 100.7]:
    start = time.time()
    integral = int(a_1)
    #sample_int = np.zeros(1000)
    rejection = np.zeros(1000)
    #for j in range(1000):
    #    for i in range(integral):
     #       U_1 = np.random.uniform(0, 1)
      #      sample_int[j] = sample_int[j] + np.log(U_1)
    sample_int = generate(integral, 1, 1000) # the sample from
    delta = a_1 - integral # The delta for which it remains to be found
    import math
    #  Ahrens-Dieter acceptance–rejection method
    sample_float = np.zeros(1000)
    e = math.exp(1)
    for i in range(1000):
        reject = 0
        while True:
            U = -1 * np.random.uniform(-1, 0)
            V = -1 * np.random.uniform(-1, 0)
            W = -1 * np.random.uniform(-1, 0)
            if(U <= e/(e+delta)):
                epsilon = math.pow(V, 1/delta)
                eta = W * math.pow(epsilon, delta - 1)
            else:
                epsilon = 1 - np.log(V)
                eta = W * math.exp(-1 * epsilon)
            if eta <= math.pow(epsilon, delta - 1) * math.exp(-1 * epsilon):
                sample_float[i] = epsilon
                break
            else:
                reject+=1
        rejection[i] = reject
    sample = sample_int + sample_float
    # The generation part ends here
    end = time.time()
    time_list.append(end-start)
    if (a_1 == 100.7):
        x = np.linspace(60, 140, 100)
    else:
        x = np.linspace (0, 10, 100)
    y1 = stat.gamma.pdf(x, a= a_1, scale=1)
    plt.plot(x, y1, "y-")
    plt.hist(sample, bins = 20, normed = True)
    plt.show()
    print("Arhens ",a_1)
    print(np.sum(rejection)/ (1000 + np.sum(rejection)))
print(time_list)

# Ratio of Uniform Method
# To be proved in the report that f(x) and x^2f(x) are bounded

for alpha in [1.5, 3.2, 100.7]:
    start = time.time()
    samples = np.zeros(1000)
    rejection = np.zeros(1000)
    for i in range(1000):
        reject = 0
        A = math.sqrt(stat.gamma.pdf(alpha - 1, a = alpha, scale = 1)) # Max value of sqrt(f(x))
        B = (alpha + 1) * math.sqrt(stat.gamma.pdf(alpha + 1, a = alpha, scale = 1)) # Max of x * sqrt(f(x))
        while True:
            U = np.random.uniform(0, A)
            V = np.random.uniform(-1 * B, B)
            if(math.pow(U, 2) <= stat.gamma.pdf(V/U, a = alpha, scale = 1)):
                samples[i] = V/U
                break
            else:
                reject += 1
        rejection[i] = reject
    end = time.time()
    if (alpha == 100.7):
        x = np.linspace(60, 140, 100)
    else:
        x = np.linspace (0, 10, 100)
    y1 = stat.gamma.pdf(x, a= alpha, scale=1)
    plt.plot(x, y1, "y-")
    plt.hist(samples, bins = 20, normed = True)
    plt.show()
    print("Ratio", alpha)
    print(np.sum(rejection)/ (1000 + np.sum(rejection)))
    print(end-start) 
    
 # Marsaglia and Tsang’s Method
def generate_normal(): # Box Muller Transform
    U_1 = np.random.uniform(0, 1)
    U_2 = np.random.uniform(0, 1)
    R = math.sqrt(-2 * np.log(U_1))
    theta = 2 * math.pi * U_2
    return R * math.cos(theta)

for alpha in [1.5, 3.2, 100.7]:
    start = time.time()
    sample = np.zeros(1000)
    rejection = np.zeros(1000)
    for i in range(1000):
        d = alpha - (1/3)
        c = 1/ math.sqrt(9*d)
        while True:
            Z = generate_normal()
            U = np.random.uniform(0, 1)
            V = math.pow(1 + (c * Z), 3)
            if((Z > -1 / c) & (np.log(U) < math.pow(Z, 2)/2 + d - d * V + d * np.log(V))):
                sample[i] = d * V
                break
            else:
                reject += 1
        rejection[i] = reject
    end = time.time()
    if (alpha == 100.7):
        x = np.linspace(60, 140, 100)
    else:
        x = np.linspace (0, 10, 100)
    y1 = stat.gamma.pdf(x, a= alpha, scale=1)
    plt.plot(x, y1, "y-")
    plt.hist(sample, bins = 20, normed = True)
    plt.show()
    print("Tsang", alpha)
    print(np.sum(rejection) / (1000 + np.sum(rejection)))
    print(end-start) 
          
