import numpy as np

def Lambda(m):
    nodes = np.linspace(-1, 1, m+1) # m + 1 equally spaced points
    
    x_values = np.linspace(-1, 1, 101) # 101 equally spaced points
    
    max_sum = 0  #initialising max_sum values
    
    for x in x_values:   #calculating L_k for each 101 x values between [-1,1]
        
        sum_l_k = 0          
        for k in range(m + 1):
            L_k = 1
            for j in range(m + 1):
             
                if j !=k:
                    L_k *= (x - nodes[j]) / (nodes[k] - nodes[j])   #The definition of L_k
            sum_l_k += abs(L_k) #calulating the sum of absolute values
            
        max_sum = max(max_sum, sum_l_k)
    return max_sum

def f(x): #Defining this function to compare with Lambda
    return 2**x / x

def LambdaC(m): #Defining the function in exactly the same way but with different nodes
    
    nodes = [np.cos(np.pi * i /m) for i in range(m+1)] #Defining the nodes as required in the question
    
    x_values = np.linspace(-1, 1, 301)
    
    max_sum = 0  #initialising max_sum values
    
    for x in x_values:   #calculating L_k for each 101 x values between [-1,1]
        
        sum_l_k = 0          
        for k in range(m + 1):
            L_k = 1
            for j in range(m + 1):
             
                if j !=k:
                    L_k *= (x - nodes[j]) / (nodes[k] - nodes[j])   #The definition of L_k
            sum_l_k += abs(L_k) #calulating the sum of absolute values
            
        max_sum = max(max_sum, sum_l_k)
    return max_sum

def g(x): #Defining this function in order to compare with LambdaC
    return 2/np.pi * np.log(x)
