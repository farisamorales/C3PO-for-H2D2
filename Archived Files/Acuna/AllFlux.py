"""
JonathanAcuna
"""
import math
import numpy as np
import matplotlib.pyplot as plt
"""CONSTANTS"""
pi = math.pi
e = math.e
h = 6.62607004*10**(-34)
k = 1.38064852*10**(-23)
c = 299792458   
"""CONSTANTS"""
"""OPENS THE DESIRED DATA FILE"""
sample1 = open('C:\Users\jma48203\Desktop\Thesis\Initial data\DDisk_17jun13_allIRflux.txt' , 'r')
s1=sample1.readlines()
s1 = np.array(map(lambda x: x.strip(), s1))
d = np.array(map(lambda x: x.split(), s1)).reshape((-1,3))
d = np.delete(d, [0,1,2]).reshape((-1,3))
d = d.astype(float)
x = []
y = []
z = []
for i in d:
    x.append(i[0])
    y.append(i[1])
    z.append(i[2])
X = np.array(x)
Y = np.array(y)
Z = np.array(z)
"""DESIRED DATA FILE"""
"""WAVE LENGTH RANGE"""
wave = np.arange(1e-6, 500e-6, 1e-6) 
"""WAVE LENGTH RANGE"""
"""PLANK FUNCTION"""
def planck(Wav, T=100, Norm=1):
    a = 2*h*c
    b = h*c/(Wav*k*T)
    Flux = (a/Wav**3)*(1/(e**(b)-1))*10**-19*Norm
    return Flux
"""PLANK FUNCTION END"""
"""STAR VALUES"""
T_star = 10000
"""STAR VALUES"""
"""INITIAL FLUX OUTPUTS"""
Flux_S = planck(wave,T_star)
Flux_C = planck(wave)
Flux_H = planck(wave)
Flux_Tot = Flux_S + Flux_C + Flux_H
"""INITIAL FLUX OUTPUTS"""
"""INITIAL CHI VALUE"""
result = []
for i in x:
    Star = planck(i*1e-6, T_star)
    Cold = planck(i*1e-6)
    Hot = planck(i*1e-6)
    result.append(Star + Cold + Hot)
Result = np.array(result)
Chi = sum(((Result - Y)**2)/(Z**2))
Re_Chi = Chi/(len(x)-4)
"""INITIAL CHI VALUE"""
"""PARAMETERS TO BE OPTIMIZED"""
Norm = 1
T_cold = 100
Norm_cold = 1
T_hot = 100
Norm_hot = 1
"""PARAMETERS TO BE OPTIMIZED"""
"""OPTIMIZATION ROUTINE"""
"""LISTS"""
CHI_LIST_NORM = []#                                                            CHI LIST NORM
NORM_LIST = []#                             PARAMETER   NORM                   NORM LIST
A_NORM_LIST = []#                                                              A LIST NORM

CHI_LIST_NORM_COLD = []#                                                       CHI LIST NORM COLD
NORM_COLD_LIST = []#                        PARAMETER   NORM COLD              NORM COLD LIST
A_NORM_COLD_LIST = []#                                                         A NORM COLD LIST

CHI_LIST_T_COLD = []#                                                          CHI LIST T COLD
T_COLD_LIST = []#                           PARAMETER   T COLD                 T COLD
A_T_COLD_LIST = []#                                                            A T COLD LIST

CHI_LIST_NORM_HOT = []#                                                        CHI LIST NORM HOT
NORM_HOT_LIST = []#                         PARAMETER   NORM HOT               NORM HOT LIST
A_NORM_HOT_LIST = []#                                                          A NORM HOT

CHI_LIST_T_HOT = []#                                                           T HOT CHI LIST
LIST_T_HOT = []#                            PARAMETER   T HOT                  T HOT LIST
A_T_HOT_LIST = []#                                                             A T HOT LIST

"""-----------------------------------"""
"""-----STAR NOMALIZER AND LOCK-------"""
"""-----------------------------------"""
loop = 0
while loop < 100:
    loop += 1
    print loop    
    
    result = []#                                                               NORM
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi_Norm = sum(((Result - Y)**2)/(Z**2))
    Re_Chi_Norm = Chi_Norm/(len(x)-4)
    
    
    UP = 1.1
    DOWN = 0.9
    a = UP
    count = 0                                   
    while count < 10:                           
        count += 1
        Norm *= a
        NORM_LIST.append(Norm)
        result_loop_Norm = []#                                                 RESULT LIST
        for i in x[:10:]:                                       
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop_Norm.append(Star + Cold + Hot)
        X_ = np.array(x[:10:])
        Y_ = np.array(y[:10:])
        Z_ = np.array(z[:10:])
        Result_LOOP_NORM = np.array(result_loop_Norm)
        Chi_Norm = sum(((Result_LOOP_NORM - Y_)**2)/(Z_**2))
        Re_Chi_loop_Norm = Chi_Norm/(len(x)-4)
        
        if a == UP and Re_Chi_loop_Norm < Re_Chi_Norm:                
            a = UP
        elif a == UP and Re_Chi_loop_Norm > Re_Chi_Norm:
            a = DOWN
        elif a == DOWN and Re_Chi_loop_Norm < Re_Chi_Norm:
            a = DOWN
        elif a == DOWN and Re_Chi_loop_Norm > Re_Chi_Norm:
            a = UP
        
        CHI_LIST_NORM.append(Re_Chi_loop_Norm)
        A_NORM_LIST.append(a)
        Re_Chi_Norm = Re_Chi_loop_Norm 



loop = 0
while loop < 100:
    loop += 1
    print loop   
    
    result = []#                                                               NORM PRECISION
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi_Norm = sum(((Result - Y)**2)/(Z**2))
    Re_Chi_Norm = Chi_Norm/(len(x)-4)
    
    
    UP = 1.001
    DOWN = 0.999
    a = UP
    count = 0                                   
    while count < 10:                           
        count += 1
        Norm *= a
        NORM_LIST.append(Norm)
        result_loop_Norm = []#                                                 RESULT LIST
        for i in x[:10:]:                                       
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop_Norm.append(Star + Cold + Hot)
        X_ = np.array(x[:10:])
        Y_ = np.array(y[:10:])
        Z_ = np.array(z[:10:])
        Result_LOOP_NORM = np.array(result_loop_Norm)
        Chi_Norm = sum(((Result_LOOP_NORM - Y_)**2)/(Z_**2))
        Re_Chi_loop_Norm = Chi_Norm/(len(x)-4)
        
        if a == UP and Re_Chi_loop_Norm < Re_Chi_Norm:                
            a = UP
        elif a == UP and Re_Chi_loop_Norm > Re_Chi_Norm:
            a = DOWN
        elif a == DOWN and Re_Chi_loop_Norm < Re_Chi_Norm:
            a = DOWN
        elif a == DOWN and Re_Chi_loop_Norm > Re_Chi_Norm:
            a = UP
        
        CHI_LIST_NORM.append(Re_Chi_loop_Norm)
        A_NORM_LIST.append(a)
        Re_Chi_Norm = Re_Chi_loop_Norm 
"""-----------------------------------"""
"""-----------------------------------"""
"""-----------------------------------"""
loop_Norms = 0
while loop_Norms < 100:
    loop_Norms += 1
    print loop_Norms 

    result = []#                                                               NORM COLD
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi = sum(((Result - Y)**2)/(Z**2))
    Re_Chi = Chi/(len(x)-4)

    UP = 1.1
    DOWN = 0.9
    a = UP
    count = 0                                  
    while count < 10:                            
        count += 1
        Norm_cold *= a
        NORM_COLD_LIST.append(Norm_cold)
        result_loop_Norm_cold = []#                                            RESULT LIST NORM COLD
        for i in x:                                      
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop_Norm_cold.append(Star + Cold + Hot)
        Result = np.array(result_loop_Norm_cold)
        Chi_Norm = sum(((Result - Y)**2)/(Z**2))
        Re_Chi_loop = Chi_Norm/(len(x)-4)

        if a == UP and Re_Chi_loop < Re_Chi:                
            a = UP
        elif a == UP and Re_Chi_loop > Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop < Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop > Re_Chi:
            a = UP
    
        CHI_LIST_NORM_COLD.append(Re_Chi_loop)
        A_NORM_COLD_LIST.append(a)
        Re_Chi = Re_Chi_loop
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------       
    result = []#                                                               NORM HOT
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi = sum(((Result - Y)**2)/(Z**2))
    Re_Chi = Chi/(len(x)-4)

    UP = 1.1
    DOWN = 0.9
    a = UP
    count = 0
    while count < 10:
        count += 1
        Norm_hot *= a
        NORM_HOT_LIST.append(Norm_hot)
        result_loop = []#                                                      RESULT NORM HOT
        for i in x:
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop.append(Star + Cold + Hot)
        Result = np.array(result_loop)
        Chi_Norm = sum(((Result - Y)**2)/(Z**2))
        Re_Chi_loop = Chi_Norm/(len(x)-4)

        if a == UP and Re_Chi_loop < Re_Chi:               
            a = UP
        elif a == UP and Re_Chi_loop > Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop < Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop > Re_Chi:
            a = UP
    
        CHI_LIST_NORM_HOT.append(Re_Chi_loop)
        A_NORM_HOT_LIST.append(a)
        Re_Chi = Re_Chi_loop   
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
T_loop = 0
while T_loop < 100:
    T_loop += 1
    print T_loop        
        
        
    result = []#                                                               T HOT
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi = sum(((Result - Y)**2)/(Z**2))
    Re_Chi = Chi/(len(x)-4)

    UP = 1.1
    DOWN = 0.9
    a = UP
    count = 0
    while count < 10:
        count += 1
        T_cold *= a
        T_COLD_LIST.append(T_cold)
        result_loop = []#                                                      RESULT T HOT
        for i in x:
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop.append(Star + Cold + Hot)
        Result = np.array(result_loop)
        Chi_Norm = sum(((Result - Y)**2)/(Z**2))
        Re_Chi_loop = Chi_Norm/(len(x)-4)

        if a == UP and Re_Chi_loop < Re_Chi:
            a = UP
        elif a == UP and Re_Chi_loop > Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop < Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop > Re_Chi:
            a = UP
        
        A_T_COLD_LIST.append(a)
        CHI_LIST_T_COLD.append(Re_Chi_loop)
        Re_Chi = Re_Chi_loop 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
    result = []#                                                               T HOT
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi = sum(((Result - Y)**2)/(Z**2))
    Re_Chi = Chi/(len(x)-4)

    UP = 1.1
    DOWN = 0.9
    a = UP
    count = 0
    while count < 10:
        count += 1
        T_hot *= a
        LIST_T_HOT.append(T_hot)
        result_loop = []#                                                      RESULT T HOT
        for i in x:
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop.append(Star + Cold + Hot)
        Result = np.array(result_loop)
        Chi_Norm = sum(((Result - Y)**2)/(Z**2))
        Re_Chi_loop = Chi_Norm/(len(x)-4)

        if a == UP and Re_Chi_loop < Re_Chi:
            a = UP
        elif a == UP and Re_Chi_loop > Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop < Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop > Re_Chi:
            a = UP
       
        A_T_HOT_LIST.append(a)
        CHI_LIST_T_HOT.append(Re_Chi_loop)
        Re_Chi = Re_Chi_loop 
        

#                       LOOP 2 FOR Precision
loop_2 = 0
while loop_2 < 200:
    loop_2 += 1
    print loop_2                            

       
    result = []#                                                               NORM COLD
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi = sum(((Result - Y)**2)/(Z**2))
    Re_Chi = Chi/(len(x)-4)

    UP = 1.01
    DOWN = 0.99
    a = UP
    count = 0                                  
    while count < 10:                            
        count += 1
        Norm_cold *= a
        NORM_COLD_LIST.append(Norm_cold)
        result_loop_Norm_cold = []#                                            RESULT LIST NORM COLD
        for i in x:                                      
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop_Norm_cold.append(Star + Cold + Hot)
        Result = np.array(result_loop_Norm_cold)
        Chi_Norm = sum(((Result - Y)**2)/(Z**2))
        Re_Chi_loop = Chi_Norm/(len(x)-4)

        if a == UP and Re_Chi_loop < Re_Chi:                
            a = UP
        elif a == UP and Re_Chi_loop > Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop < Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop > Re_Chi:
            a = UP
    
        CHI_LIST_NORM_COLD.append(Re_Chi_loop)
        A_NORM_COLD_LIST.append(a)
        Re_Chi = Re_Chi_loop
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------       
    result = []#                                                               NORM HOT
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi = sum(((Result - Y)**2)/(Z**2))
    Re_Chi = Chi/(len(x)-4)

    UP = 1.01
    DOWN = 0.99
    a = UP
    count = 0
    while count < 10:
        count += 1
        Norm_hot *= a
        NORM_HOT_LIST.append(Norm_hot)
        result_loop = []#                                                      RESULT NORM HOT
        for i in x:
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop.append(Star + Cold + Hot)
        Result = np.array(result_loop)
        Chi_Norm = sum(((Result - Y)**2)/(Z**2))
        Re_Chi_loop = Chi_Norm/(len(x)-4)

        if a == UP and Re_Chi_loop < Re_Chi:               
            a = UP
        elif a == UP and Re_Chi_loop > Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop < Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop > Re_Chi:
            a = UP
    
        CHI_LIST_NORM_HOT.append(Re_Chi_loop)
        A_NORM_HOT_LIST.append(a)
        Re_Chi = Re_Chi_loop   
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
T_loop_2 = 0
while T_loop_2 < 100:
    T_loop_2 += 1
    print T_loop_2        
        
        
    result = []#                                                               T HOT
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi = sum(((Result - Y)**2)/(Z**2))
    Re_Chi = Chi/(len(x)-4)

    UP = 1.01
    DOWN = 0.99
    a = UP
    count = 0
    while count < 10:
        count += 1
        T_cold *= a
        T_COLD_LIST.append(T_cold)
        result_loop = []#                                                      RESULT T HOT
        for i in x:
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop.append(Star + Cold + Hot)
        Result = np.array(result_loop)
        Chi_Norm = sum(((Result - Y)**2)/(Z**2))
        Re_Chi_loop = Chi_Norm/(len(x)-4)

        if a == UP and Re_Chi_loop < Re_Chi:
            a = UP
        elif a == UP and Re_Chi_loop > Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop < Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop > Re_Chi:
            a = UP
        
        A_T_COLD_LIST.append(a)
        CHI_LIST_T_COLD.append(Re_Chi_loop)
        Re_Chi = Re_Chi_loop 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
    result = []#                                                               T HOT
    for i in x:
        Star = planck(i*1e-6, T_star, Norm)
        Cold = planck(i*1e-6, T_cold, Norm_cold)
        Hot = planck(i*1e-6, T_hot, Norm_hot)
        result.append(Star + Cold + Hot)
    Result = np.array(result)
    Chi = sum(((Result - Y)**2)/(Z**2))
    Re_Chi = Chi/(len(x)-4)

    UP = 1.01
    DOWN = 0.99
    a = UP
    count = 0
    while count < 10:
        count += 1
        T_hot *= a
        LIST_T_HOT.append(T_hot)
        result_loop = []#                                                      RESULT T HOT
        for i in x:
            Star = planck(i*1e-6, T_star, Norm)
            Cold = planck(i*1e-6, T_cold, Norm_cold)
            Hot = planck(i*1e-6, T_hot, Norm_hot)
            result_loop.append(Star + Cold + Hot)
        Result = np.array(result_loop)
        Chi_Norm = sum(((Result - Y)**2)/(Z**2))
        Re_Chi_loop = Chi_Norm/(len(x)-4)

        if a == UP and Re_Chi_loop < Re_Chi:
            a = UP
        elif a == UP and Re_Chi_loop > Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop < Re_Chi:
            a = DOWN
        elif a == DOWN and Re_Chi_loop > Re_Chi:
            a = UP
       
        A_T_HOT_LIST.append(a)
        CHI_LIST_T_HOT.append(Re_Chi_loop)
        Re_Chi = Re_Chi_loop 
         
     
"""------------------------END LOOP 2---------------------------------------"""   
        

        

 
"""OPTIMIZATION ROUTINE"""
"""INITIAL FLUX OUTPUTS"""
Flux_S = planck(wave,T_star, Norm)
Flux_C = planck(wave, T_cold, Norm_cold)
Flux_H = planck(wave, T_hot, Norm_hot)
Flux_Tot = Flux_S + Flux_C + Flux_H
"""INITIAL FLUX OUTPUTS"""


plt.hold(True)
plt.loglog(wave*1e6,Flux_S, 'y-')
plt.loglog(wave*1e6,Flux_C, 'b-')
plt.loglog(wave*1e6,Flux_H, 'r-')
plt.loglog(wave*1e6,Flux_Tot, 'g--')
plt.loglog(x,y, 'k.')
plt.errorbar(x,y,yerr=z,linestyle="none")



plt.title('20')
"""PLOT LIMITS"""
plt.xlim([0.1,520])
plt.ylim([0.0000001,100])
"""PLOT LIMITS"""



plt.show()