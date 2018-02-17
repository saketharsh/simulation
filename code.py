import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import matplotlib.pyplot as plt
import math
import time
import scipy.stats as stat

time_list = []



np.random.seed(157) 
shape_paremeter = 1  

##########################################################################################################################################

#  Integral Gamma Generation , that can generate Exp and then sum and send the results 
def expGen(alpha, beta, size): # size is the number of sample
sample = np.zeros(size)
for i in range(size):
	U = np.random.uniform(0, 1, alpha)
	lnU = np.log(U)
	sample[i] = (-1 * np.sum(lnU))/beta
	return sample


def expGen_Normal(): # Box Muller Transform   
	U_1 = np.random.uniform(0, 1)
	U_2 = np.random.uniform(0, 1)
	R = math.sqrt(-2 * np.log(U_1))
	theta = 2 * math.pi * U_2
	return R * math.cos(theta)






#########################################################################################################################################
def Ratio_Uniform (shape_param):
	init_time = time.time()
	accepted_samples = []
	rejected_samples = []
	num_accepted =0
	reject = 0
	i=-1
    A = math.sqrt(stat.gamma.pdf(shape_param - 1, a = shape_param, scale = 1)) # Max value of sqrt(f(x))
    B = (alpha + 1) * math.sqrt(stat.gamma.pdf(shape_param + 1, a = alpha, scale = 1)) # Max of x * sqrt(f(x))
    while num_accepted is not 1000:        
    	U = np.random.uniform(0, A)
    	V = np.random.uniform(-1 * B, B)
    	if(math.pow(U, 2) <= stat.gamma.pdf(V/U, a = shape_param, scale = 1)):
    		accepted_samples.append(V/U);
    	else:
    		reject += 1
    		rejected_samples.append(reject)
    end_time = time.time()
	if (alpha == 100.7):
		x = np.linspace(60, 140, 100)
	else:
		x = np.linspace (0, 10, 100)
    y1 = stat.gamma.pdf(x, a= alpha, scale=1)
    plt.plot(x, y1, "y-")
    plt.hist(accepted_samples, bins = 20, normed = True)
	plt.show()
	print(np.mean(rejected_samples))
	print(end_time - init_time) 

######################################################################################################################################


def MSMethod(shape_param):
	init_time = time.time()
	accepted_samples = []
	rejected_samples= []
	num_accepted =0
	reject =0
	i=-1
	D = shape_param - (1/3)
	C = 1/math.sqrt(9*D)
	while num_accepted is not 1000:
		Z = expGen_Normal()
		U = np.random.uniform(0, 1)
		V = math.pow(1 + (C * Z), 3)
		if((Z > -1 / C) & (np.log(U) < math.pow(Z, 2)/2 + D -  D * V +  D * np.log(V))):
			accepted_samples.append(D*V)
		else:
			reject += 1
			rejected_samples.append(reject)
	end = time.time()
	if (alpha == 100.7):
		x = np.linspace(60, 140, 100)
	else:
		x = np.linspace (0, 10, 100)
	y1 = stat.gamma.pdf(x, a= alpha, scale=1)
	plt.plot(x, y1, "y-")
	plt.hist(sample, bins = 20, normed = True)
	plt.show()
	print(np.mean(rejection))
	print(end-start) 


##############################################################################################################################################







# sample_array = np.empty(1000) # The array to store the output


# for shape_param in [1.5, 3.2, 100.7]:
#     start = time.time()
#     rejection = np.zeros(1000)
#     sample_int = expGen(int(shape_param) , 1, 1000) # the sample from genExp that gives integral part 
#     delta = shape_param - int(shape_param)  # The delta for which it remains to be found




#     #  Ahrens-Dieter acceptance–rejection method
#     sample_float = np.zeros(1000)
#     e = math.exp(1)
#     for i in range(1000):
#         reject = 0
#         while True:
#             U = -1 * np.random.uniform(-1, 0)
#             V = -1 * np.random.uniform(-1, 0)
#             W = -1 * np.random.uniform(-1, 0)
#             if(U <= e/(e+delta)):
#                 epsilon = math.pow(V, 1/delta)
#                 eta = W * math.pow(epsilon, delta - 1)
#             else:
#                 epsilon = 1 - np.log(V)
#                 eta = W * math.exp(-1 * epsilon)
#             if eta <= math.pow(epsilon, delta - 1) * math.exp(-1 * epsilon):
#                 sample_float[i] = epsilon
#                 break
#             else:
#                 reject+=1
#         rejection[i] = reject
#     sample = sample_int + sample_float
#     # The generation part ends here
#     end = time.time()
#     time_list.append(end-start)
#     if (shape_param == 100.7):
#         x = np.linspace(60, 140, 100)
#     else:
#         x = np.linspace (0, 10, 100)
#     y1 = stat.gamma.pdf(x, a= shape_param, scale=1)
#     plt.plot(x, y1, "y-")
#     plt.hist(sample, bins = 20, normed = True)
#     plt.show()
#     print(np.mean(rejection))
# print(time_list)



# Ratio of Uniform Method
for alpha in [1.5, 3.2, 100.7]:
	Ratio_Uniform(alpha)
	print("Ratio Of uniform Method ends here for shape parameter " ++ alpha)
	print("----------------------------------------------------------------")




# Marsaglia and Tsang's Method 
for alpha in [1.5, 3.2, 100.7]:
	MSMethod(alpha)
	print("MSMethod ends here , note all observations " ++ alpha)



