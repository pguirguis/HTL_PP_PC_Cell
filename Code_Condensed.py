# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 10:26:02 2021

Full model

@author: Peter M Guirguis
"""


import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import pandas as pd
import math
import matplotlib.pyplot as plt



k_n1 = 11
k_n2 = 5
k_n3 = 9

df = pd.read_excel('Kinetics_data_calculated.xls', sheet_name = 'Summary',index_col=None, header=None)


df = df.to_numpy()
x0 = [100, 0, 0, 0]
sx0 = [0,0,0,0]
t = np.append(0,df[3:9,0])


# given data we want to fit

Cell_data0 = np.vstack((x0,df[3:9,1:5] ))
PP_data0 = np.vstack((x0,df[3:8,5:9] ))
PC_data0 = np.vstack((x0,df[3:9,9:13] ))
Cell_sigma0 = np.vstack((sx0,df[11:17,1:5] ))
PP_sigma0 = np.vstack((sx0,df[11:16,5:9] ))
PC_sigma0 = np.vstack((sx0,df[11:17,9:13] ))
PP_Cell_data0 = np.vstack((x0,df[3:9,13:17] ))
PC_Cell_data0 = np.vstack((x0,df[3:9,17:21] ))
PP_PC_data0 = np.vstack((x0,df[3:9,21:25] ))
PP_PC_Cell_data0 = np.vstack((x0,df[3:9,25:29] ))
PP_Cell_sigma0 = np.vstack((sx0,df[11:17,13:17] ))
PC_Cell_sigma0 = np.vstack((sx0,df[11:17,17:21] ))
PP_PC_sigma0 = np.vstack((sx0,df[11:17,21:25] ))
PP_PC_Cell_sigma0 = np.vstack((sx0,df[11:17,25:29] ))

t = np.reshape(t,-1)

df = pd.read_excel('Ratio_effects_results_10_6.xls', sheet_name = 'PP - Cell',index_col=None, header=None)
df = df.to_numpy()
PP_Cell_percent = df[2:9,0]
PP_Cell_B = df[2:9,1] 
PP_Cell_sigma = df[2:9,2] 
df = pd.read_excel('Ratio_effects_results_10_6.xls', sheet_name = 'PC - Cell',index_col=None, header=None)
df = df.to_numpy()
PC_Cell_percent = df[2:9,0] 
PC_Cell_B = df[2:9,1] 
PC_Cell_sigma = df[2:9,2] 
df = pd.read_excel('Ratio_effects_results_10_6.xls', sheet_name = 'PP-PC',index_col=None, header=None)
df = df.to_numpy()
PC_PP_percent = df[2:9,0] 
PC_PP_B = df[2:9,1] 
PC_PP_sigma = df[2:9,2] 



def coupled(k,a,b,x_S1,x_S2):
    if not(math.isfinite(k)) or x_S1 < 0 or x_S2 < 0:
        r =0
    else:
        r = k*pow(x_S1,a)*pow(x_S2,b)
        if not(math.isfinite(r)) :
              r = 100
    return r

def ode1(x, tspan, k):    
    
    x_S = x[0]
    x_B = x[1]
    x_Aq = x[2]

    
    k1  = k[0]
    k2  = k[1]
    k3  = k[2]
    k4  = k[3]
    k5  = k[4]
    k6  = k[5]
    k7  = k[6]
    k8  = k[7]
    k9  = k[8]
    k10 = k[9]
    k11 = k[10]
    
    r1  = k1*x_S  # S -> B
    r2  = k2*x_B  # B -> S
    r3  = k3*x_S  # S -> Aq
    r4  = k4*x_Aq # Aq-> S
    r5  = k5*x_Aq # Aq-> B
    r6  = k6*x_B  # B -> Aq
    r7  = k7*x_B  # B -> G
    r8  = k8*x_Aq # Aq-> G
    r9  = k9*x_S  # S -> G
    r10 = k10*x_B # B -> S_new
    r11 = k11*x_Aq# Aq-> S_new  
    
    sol = [-r1-r3-r9+r2+r4, r1+r5-r2-r6-r7-r10, r3+r6-r4-r5-r8-r11, r7+r8+r9, r10+r11]
    # [S, B, Aq, G, S_new]
    return sol

def ode2(x, tspan, ks,k):    
        x_S1 = x[0]
        x_S2 = x[1]
        x_B_S1 = x[2]
        x_B_S2 = x[3]
        x_Aq_S1 = x[4];
        x_Aq_S2 = x[5];

        
        k1 = ks[0]
        k2 = ks[1]
        k3 = ks[2]
        k4 = ks[3]
        k5 = ks[4]
        k6 = ks[5]
        k7 = ks[6]
        k8 = ks[7]
        k9 = ks[8]
        k10 = ks[9]
        k11 = ks[10]
        k12 = ks[11]
        k13 = ks[12]
        k14 = ks[13]
        k15 = ks[14]
        k16 = ks[15]
        k17 = ks[16]
        k18 = ks[17]
        k19 = ks[18]
        k20 = ks[19]
        k21 = ks[20]
        k22 = ks[21]
        k23 = k[0]
        a = k[1]
        b = k[2]
        y = k[3]
        z = k[4]

        
        r1  = k1*x_S1     # S1  -> B1
        r2  = k2*x_B_S1   # B1  -> S1
        r3  = k3*x_S1     # S1  -> Aq1
        r4  = k4*x_Aq_S1  # Aq1 -> S1
        r5  = k5*x_Aq_S1  # Aq1 -> B1
        r6  = k6*x_B_S1   # B1  -> Aq1
        r7  = k7*x_B_S1   # B1  -> G1
        r8  = k8*x_Aq_S1  # Aq1 -> G1
        r9  = k9*x_S1     # S1  -> G1
        r10 = k10*x_B_S1  # B1  -> S_new1
        r11 = k11*x_Aq_S1 # Aq1 -> S_new1
        r12 = k12*x_S2    # S2  -> B2
        r13 = k13*x_B_S2  # B2  -> S2
        r14 = k14*x_S2    # S2  -> Aq2
        r15 = k15*x_Aq_S2 # Aq2 -> S2
        r16 = k16*x_Aq_S2 # Aq2 -> B2
        r17 = k17*x_B_S2  # B2  -> Aq2
        r18 = k18*x_B_S2  # B2  -> G2
        r19 = k19*x_Aq_S2 # Aq2 -> G2
        r20 = k20*x_S2    # S2  -> G2
        r21 = k21*x_B_S2  # B2  -> S_new2
        r22 = k22*x_Aq_S2 # Aq2 -> S_new2
        
        r23 = coupled(k23,a,b,x_S1,x_S2) # S1 + S2 -> B12
        
        sol = [-r1+r2-r3+r4-r9-(y*r23),-r12+r13-r14+r15-r20-(z*r23), r1+r5-r2-r6-r7-r10, r12-r13+r16-r17-r18-r21, r3+r6-r4-r5-r8-r11, r14+r17-r15-r16-r19-r22,  r23, r7+r8+r9,r18+r19+r20,r10+r11,r21+r22]
        # [S1, S2, B1, B2, Aq1, Aq2, B12, G1, G2, S_new1, S_new2]
        return sol

def ode3(x, tspan, ks, k):  
    x_S1 = x[0]
    x_S2 = x[1]
    x_B_S1 = x[2]
    x_B_S2 = x[3]
    x_B_S12 = x[4]
    x_Aq_S1 = x[5]
    x_Aq_S2 = x[6]
    x_Aq_S12 = x[7]

    
    k1 = ks[0]
    k2 = ks[1]
    k3 = ks[2]
    k4 = ks[3]
    k5 = ks[4]
    k6 = ks[5]
    k7 = ks[6]
    k8 = ks[7]
    k9 = ks[8]
    k10 = ks[9]
    k11 = ks[10]
    k12 = ks[11]
    k13 = ks[12]
    k14 = ks[13]
    k15 = ks[14]
    k16 = ks[15]
    k17 = ks[16]
    k18 = ks[17]
    k19 = ks[18]
    k20 = ks[19]
    k21 = ks[20]
    k22 = ks[21]
    a = ks[22]
    b = ks[23]
    y = ks[24]
    z = ks[25]
    
    k23 = k[0]
    k24 = k[1]
    k25 = k[2]
    k26 = k[3]
    k27 = k[4]
    k28 = k[5]
    k29 = k[6]
    k30 = k[7]
    k31 = k[8]
    
    
    r1  = k1*x_S1        # S1  -> B1
    r2  = k2*x_B_S1      # B1  -> S1
    r3  = k3*x_S1        # S1  -> Aq1
    r4  = k4*x_Aq_S1     # Aq1 -> S1
    r5  = k5*x_Aq_S1     # Aq1 -> B1
    r6  = k6*x_B_S1      # B1  -> Aq1
    r7  = k7*x_B_S1      # B1  -> G1
    r8  = k8*x_Aq_S1     # Aq1 -> G1
    r9  = k9*x_S1        # S1  -> G1
    r10 = k10*x_B_S1     # B1  -> S_new1
    r11 = k11*x_Aq_S1    # Aq1 -> S_new1
    r12 = k12*x_S2       # S2  -> B2
    r13 = k13*x_B_S2     # B2  -> S2
    r14 = k14*x_S2       # S2  -> Aq2
    r15 = k15*x_Aq_S2    # Aq2 -> S2
    r16 = k16*x_Aq_S2    # Aq2 -> B2
    r17 = k17*x_B_S2     # B2  -> Aq2
    r18 = k18*x_B_S2     # B2  -> G2
    r19 = k19*x_Aq_S2    # Aq2 -> G2
    r20 = k20*x_S2       # S2  -> G2
    r21 = k21*x_B_S2     # B2  -> S_new2
    r22 = k22*x_Aq_S2    # Aq2 -> S_new2
    r25 = k25*x_Aq_S12   # Aq12  -> B12
    r26 = k26*x_B_S12    # B12   -> Aq12
    r27 = k27*x_B_S12    # B12   -> G12
    r28 = k28*x_Aq_S12   # Aq12  -> G12
    r30 = k30*x_B_S12    # B12   -> S_new12
    r31 = k31*x_Aq_S12   # Aq12  -> S_new12

    r24 = coupled(k24,1,1,x_S1,x_S2) # S1+S2 -> Aq12
    r29 = coupled(k29,1,1,x_S1,x_S2) # S1+S2 -> G12
    r23 = coupled(k23,a,b,x_S1,x_S2) # S1 + S2 -> B12
    

    sol = [-r1+r2-r3+r4-r9-(y*r23)-r24-r29, -r12+r13-r14+r15-r20-(z*r23)-r24-r29, r1+r5-r2-r6-r7-r10, r12-r13+r16-r17-r18-r21,r23+r25-r26-r27-r30, r3+r6-r4-r5-r8-r11, r14+r17-r15-r16-r19-r22, r24-r25+r26-r28-r31, r7+r8+r9,r18+r19+r20,r27+r28+r29,r10+r11,r21+r22,r30+r31]
    # S1 , S2, B1, B2, B12, Aq1, Aq2, Aq12, G1, G2, G12, S_new1, S_new2, S_new12  
    return sol

def Fit_Fun1(tspan,data,g,s):
    
    row,col = np.shape(data)
    t0 = 0
    tf = tspan[-1]
    tspan = np.asarray(tspan, dtype = np.float64, order ='C')
    s = np.asarray(s, dtype = np.float64, order ='C')
    data = np.asarray(data, dtype = np.float64, order ='C')

   

    def my_ls_func1(tspan,*teta):
        """definition of function for LS fit
            x gives evaluation points,
            teta is an array of parameters to be varied for fit"""
        # create an alias to f which passes the optional params    
        f2 = lambda t,y: ode1(y, t, teta)
        # calculate ode solution, retuen values for each entry of "x"
    
        r = solve_ivp(f2,(t0,tf),data0, t_eval=tspan,method='Radau')
        #in this case, we only need one of the dependent variable values
        pred = np.transpose(r.y)
        pred[:,0] = np.add(pred[:,0],pred[:,4])
        pred = pred[:,:4]
        pred = np.reshape(pred,-1)
        pred = np.asarray(pred, dtype = np.float64, order ='C')
        return pred
    

    #solve the system - the solution is in variable c
    # l =1 
    data0 = np.append(data[0],0) #inital conditions for ODEs
 
    # o = np.subtract( data, my_ls_func(tspan, guess))
    data = np.reshape(data,-1)
    s = np.reshape(s, -1)
    for q in range(len(s)):
        if s[q] ==0 or not(np.isfinite(s[q])):
            s[q] = 1e-15
    k, var = curve_fit(my_ls_func1,tspan, data,g,s, bounds = (0,1))#get params

    return k , var



PP_kf,PP_k_var = Fit_Fun1(t[:-1],PP_data0,np.random.rand(k_n1)*0.00001,PP_sigma0)
PP_sd = np.sqrt(np.diag(PP_k_var))
PC_kf,PC_k_var = Fit_Fun1(t,PC_data0,np.random.rand(k_n1)*0.00001,PC_sigma0)
PC_sd = np.sqrt(np.diag(PC_k_var))
Cell_kf,Cell_k_var = Fit_Fun1(t,Cell_data0,np.random.rand(k_n1)*0.00001,Cell_sigma0)
Cell_sd = np.sqrt(np.diag(Cell_k_var))


def Fit_Fun2(datain,dataout,g,ks,s):
    
    ks = [item for subl in ks for item in subl]
    tf =30
    
    datain = np.asarray(datain, dtype = np.float64, order ='C')
    dataout = np.asarray(dataout, dtype = np.float64, order ='C')
    s = np.asarray(s, dtype = np.float64, order ='C')
    
    
    
    def my_ls_func2(per,*teta):
        """definition of function for LS fit
            x gives evaluation points,
            teta is an array of parameters to be varied for fit"""
        # create an alias to f which passes the optional params    
        f2 = lambda t,y: ode2(y, t,ks,teta)
        # calculate ode solution, retuen values for each entry of "x"
        pred = np.zeros(len(per))
        for q in range(len(per)):
            data0 = [per[q],100-per[q], 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
            r = solve_ivp(f2,(0,tf),data0, t_eval=[0, tf],method='BDF')
            #in this case, we only need one of the dependent variable values
            y = r.y
        
            pred[q] = y[2,-1] + y[3,-1] + y[6,-1]
            
        return pred
    

    #solve the system - the solution is in variable 
    # l =1 
  


    k, var = curve_fit(my_ls_func2,datain, dataout,g,s, bounds = (-3,3))#get params
    
    # [k, kvg] = leastsq(f_resid,guess)
       
 
    return k, var

PP_Cell_power,PP_Cell_power_var = Fit_Fun2(PP_Cell_percent,PP_Cell_B,np.zeros(k_n2),[PP_kf,Cell_kf],PP_Cell_sigma)
PP_Cell_power_sd = np.sqrt(np.diag(PP_Cell_power_var))
PC_Cell_power,PC_Cell_power_var = Fit_Fun2(PC_Cell_percent,PC_Cell_B,np.zeros(k_n2),[PC_kf,Cell_kf],PC_Cell_sigma)
PC_Cell_power_sd = np.sqrt(np.diag(PC_Cell_power_var))
PC_PP_power,PC_PP_power_var = Fit_Fun2(PC_PP_percent,PC_PP_B,np.zeros(k_n2) ,[PC_kf,PP_kf],PC_PP_sigma)
PC_PP_power_sd = np.sqrt(np.diag(PC_PP_power_var))
PP_Cell_power_sd = PP_Cell_power_sd[1:]
PC_Cell_power_sd = PC_Cell_power_sd[1:]
PP_PC_power_sd = [PC_PP_power_sd[2],PC_PP_power_sd[1],PC_PP_power_sd[4],PC_PP_power_sd[3]]


def Fit_Fun3(tspan,data,ks,g,frac,s):
    row,col = np.shape(data)
    t0 = 0
    tf = tspan[-1]
    tspan = np.asarray(tspan, dtype = np.float64, order ='C')
    data = np.asarray(data, dtype = np.float64, order ='C')
    s = np.asarray(s, dtype = np.float64, order ='C')

    
    def my_ls_func3(tspan,*teta):
        """definition of function for LS fit
            x gives evaluation points,
            teta is an array of parameters to be varied for fit"""
        # create an alias to f which passes the optional params    
        f2 = lambda t,y: ode3(y, t, ks ,teta)
        # calculate ode solution, retuen values for each entry of "x"
    
        r = solve_ivp(f2,(t0,tf),data0, t_eval=tspan,method ='Radau')
        #in this case, we only need one of the dependent variable values
        pred =  np.transpose(r.y)
        pred[:,0] = np.add(np.add(np.add(np.add(pred[:,0],pred[:,1]),pred[:,11]),pred[:,12]),pred[:,13])
        pred[:,1] = np.add(np.add(pred[:,2],pred[:,3]),pred[:,4])
        pred[:,2] = np.add(np.add(pred[:,5],pred[:,6]),pred[:,7])
        pred[:,3] = np.add(np.add(pred[:,8],pred[:,9]),pred[:,10])
        pred = pred[:,:4]
        pred = np.reshape(pred,-1)
        if len(pred) < 28:
            pred = np.ones(28)*0
        return pred

    #solve the system - the solution is in variable c
    # l =1 

    data0 = [100*frac, 100*(1-frac), 0, 0, 0,0,0,0,0,0,0,0,0,0] #inital conditions for ODEs

    data = np.reshape(data,-1)
    s = np.reshape(s,-1)
    for q in range(len(s)):
        if s[q] ==0 or not(np.isfinite(s[q])):
            s[q] = 1e-15
    k, var = curve_fit(my_ls_func3,tspan, data,g,s, bounds = (0,1))#get params


    return k ,var

PP_PC_power_n = [PC_PP_power[2],PC_PP_power[1],PC_PP_power[4],PC_PP_power[3]]
PP_Cell_power_n = PP_Cell_power[1:]
PC_Cell_power_n = PC_Cell_power[1:]
k_PP_PC = np.append(PP_kf,PC_kf)
k_PP_PC = np.append(k_PP_PC,PP_PC_power_n)
k_PP_Cell = np.append(PP_kf,Cell_kf)
k_PP_Cell = np.append(k_PP_Cell,PP_Cell_power_n)
k_PC_Cell = np.append(PC_kf,Cell_kf)
k_PC_Cell = np.append(k_PC_Cell,PC_Cell_power_n)


PP_PC_kf, PP_PC_var = Fit_Fun3(t,PP_PC_data0,k_PP_PC,np.zeros(k_n3),0.5,PP_PC_sigma0)
PP_PC_sd = np.sqrt(np.diag(PP_PC_var))
PC_Cell_kf, PC_Cell_var = Fit_Fun3(t,PC_Cell_data0,k_PC_Cell,np.zeros(k_n3),.8,PC_Cell_sigma0)
PC_Cell_sd = np.sqrt(np.diag(PC_Cell_var))
PP_Cell_kf, PP_Cell_var = Fit_Fun3(t,PP_Cell_data0,k_PP_Cell,np.zeros(k_n3),.5,PP_Cell_sigma0)
PP_Cell_sd = np.sqrt(np.diag(PP_Cell_var))



def Plot_Fun1(t,k):
    t0 = t[0]
    tf = t[-1]

        
    def ode1(x, tspan, k):    
        
        x_S = x[0]
        x_B = x[1]
        x_Aq = x[2]

        
        k1  = k[0]
        k2  = k[1]
        k3  = k[2]
        k4  = k[3]
        k5  = k[4]
        k6  = k[5]
        k7  = k[6]
        k8  = k[7]
        k9  = k[8]
        k10 = k[9]
        k11 = k[10]
        
        r1  = k1*x_S  # S -> B
        r2  = k2*x_B  # B -> S
        r3  = k3*x_S  # S -> Aq
        r4  = k4*x_Aq # Aq-> S
        r5  = k5*x_Aq # Aq-> B
        r6  = k6*x_B  # B -> Aq
        r7  = k7*x_B  # B -> G
        r8  = k8*x_Aq # Aq-> G
        r9  = k9*x_S  # S -> G
        r10 = k10*x_B # B -> S_new
        r11 = k11*x_Aq# Aq-> S_new  
        
        sol = [-r1-r3-r9+r2+r4, r1+r5-r2-r6-r7-r10, r3+r6-r4-r5-r8-r11, r7+r8+r9, r10+r11]
        # [S, B, Aq, G, S_new]
        return sol
    
    
    
    def my_ls_func1(tspan,teta):
        """definition of function for LS fit
            x gives evaluation points,
            teta is an array of parameters to be varied for fit"""
        # create an alias to f which passes the optional params    
        f2 = lambda t,y: ode1(y, t, teta)
        # calculate ode solution, retuen values for each entry of "x"
    
        r = solve_ivp(f2,(t0,tf),data0, t_eval=tspan,method='Radau')
        #in this case, we only need one of the dependent variable values
        pred = np.transpose(r.y)
        pred[:,0] = np.add(pred[:,0],pred[:,4])
        pred = pred[:,:4]
        pred = np.asarray(pred, dtype = np.float64, order ='C')
        return pred
    
    data0 = [100,0,0,0,0] #inital conditions for ODEs
    
    fit = my_ls_func1(t, k)
    
    return fit

tfit = np.linspace(0,60)

PP_fitf = Plot_Fun1(tfit, PP_kf)
PC_fitf = Plot_Fun1(tfit, PC_kf)
Cell_fitf = Plot_Fun1(tfit, Cell_kf)

titles = ['Solids','Biocrude','Aquoeous','Gas']

def plot_data_fit3(t,data,tfit,fit,s,lab):
    fig, axs = plt.subplots(2,2,constrained_layout=True,sharex=True)
    m = 0
    r, c = np.shape(fit)
    tfit = tfit[:r]
    for i in range(2):
        for j in range(2):
            axs[i,j].plot(t, data[:,m], 'ro', label='data')
            axs[i,j].plot(tfit, fit[:,m], 'b-', label='fit')
            axs[i,j].legend(loc='best')
            axs[i,j].set_title(titles[m])
            axs[i,j].errorbar(t,data[:,m], yerr=s[:,m],linestyle="None")
            m= m+1
    fig.text(0.5, 0, 'time (min)', ha='center')
    fig.text(0.04, 0.5, 'wt%', va='center', rotation='vertical')
    fig.text(0.5, 1, lab, ha='center')
    plt.show()


plot_data_fit3(t[:-1],PP_data0,tfit,PP_fitf,PP_sigma0,"PP")
plot_data_fit3(t,PC_data0,tfit,PC_fitf,PC_sigma0,"PC")
plot_data_fit3(t,Cell_data0,tfit,Cell_fitf,Cell_sigma0,"Cell")





def Plot_Fun2(ks):
    
    ks = [item for subl in ks for item in subl]
    tf =30

    def ode2(x, tspan, ks):    
        x_S1 = x[0]
        x_S2 = x[1]
        x_B_S1 = x[2]
        x_B_S2 = x[3]
        x_Aq_S1 = x[4];
        x_Aq_S2 = x[5];

        
        k1 = ks[0]
        k2 = ks[1]
        k3 = ks[2]
        k4 = ks[3]
        k5 = ks[4]
        k6 = ks[5]
        k7 = ks[6]
        k8 = ks[7]
        k9 = ks[8]
        k10 = ks[9]
        k11 = ks[10]
        k12 = ks[11]
        k13 = ks[12]
        k14 = ks[13]
        k15 = ks[14]
        k16 = ks[15]
        k17 = ks[16]
        k18 = ks[17]
        k19 = ks[18]
        k20 = ks[19]
        k21 = ks[20]
        k22 = ks[21]
        k23 = ks[22]
        a = ks[23]
        b = ks[24]
        y = ks[25]
        z = ks[26]

        
        r1  = k1*x_S1     # S1  -> B1
        r2  = k2*x_B_S1   # B1  -> S1
        r3  = k3*x_S1     # S1  -> Aq1
        r4  = k4*x_Aq_S1  # Aq1 -> S1
        r5  = k5*x_Aq_S1  # Aq1 -> B1
        r6  = k6*x_B_S1   # B1  -> Aq1
        r7  = k7*x_B_S1   # B1  -> G1
        r8  = k8*x_Aq_S1  # Aq1 -> G1
        r9  = k9*x_S1     # S1  -> G1
        r10 = k10*x_B_S1  # B1  -> S_new1
        r11 = k11*x_Aq_S1 # Aq1 -> S_new1
        r12 = k12*x_S2    # S2  -> B2
        r13 = k13*x_B_S2  # B2  -> S2
        r14 = k14*x_S2    # S2  -> Aq2
        r15 = k15*x_Aq_S2 # Aq2 -> S2
        r16 = k16*x_Aq_S2 # Aq2 -> B2
        r17 = k17*x_B_S2  # B2  -> Aq2
        r18 = k18*x_B_S2  # B2  -> G2
        r19 = k19*x_Aq_S2 # Aq2 -> G2
        r20 = k20*x_S2    # S2  -> G2
        r21 = k21*x_B_S2  # B2  -> S_new2
        r22 = k22*x_Aq_S2 # Aq2 -> S_new2
        
        r23 = coupled(k23,a,b,x_S1,x_S2) # S1 + S2 -> B12
        
        sol = [-r1+r2-r3+r4-r9-(y*r23),-r12+r13-r14+r15-r20-(z*r23), r1+r5-r2-r6-r7-r10, r12-r13+r16-r17-r18-r21, r3+r6-r4-r5-r8-r11, r14+r17-r15-r16-r19-r22,  r23, r7+r8+r9,r18+r19+r20,r10+r11,r21+r22]
        # [S1, S2, B1, B2, Aq1, Aq2, B12, G1, G2, S_new1, S_new2]
        return sol
    
    
    
    def my_ls_func2(tf,teta,data0):
        """definition of function for LS fit
            x gives evaluation points,
            teta is an array of parameters to be varied for fit"""
        # create an alias to f which passes the optional params    
        f2 = lambda t,y: ode2(y, t,teta)
        # calculate ode solution, retuen values for each entry of "x"
        t0 = 0
    
        r = solve_ivp(f2,(t0,tf),data0, t_eval=[0, tf],method='BDF')
        #in this case, we only need one of the dependent variable values
        return r.y
    
    # [k, kvg] = leastsq(f_resid,guess)
       
 
    fit_datain = np.linspace(0,100,50)
    fit_datain = fit_datain[1:-1]
    
    row = len(fit_datain)
    fit = np.zeros((row,1))
    j = 0
    for i in fit_datain:
        y = my_ls_func2(tf,ks,np.array([i , 100-i, 0,0,0,0,0,0,0,0,0]))
        fit[j] = y[2,-1] + y[3,-1] + y[6,-1]
        j = j+1
    return fit


def plot_data_fit2(t,data,tfit,fit,s,lab):
    plt.plot(t, data, 'ro', label='data')
    plt.plot(tfit, fit, 'b-', label='fit')
    plt.legend(loc='best')
    plt.ylabel('Biocrude wt%')
    plt.xlabel('Wt% of first species')
    plt.errorbar(t,data, yerr=s,linestyle="None")
    plt.title(lab)

    plt.show()

fit_per = np.linspace(0,100,50)
fit_per = fit_per[1:-1]

k_PC_PP = [PC_kf,PP_kf,PC_PP_power]
k_PC_Cell = [PC_kf,Cell_kf,PC_Cell_power]
k_PP_Cell = [PP_kf,Cell_kf,PP_Cell_power]
PC_PP_fit = Plot_Fun2(k_PC_PP)
PC_Cell_fit = Plot_Fun2(k_PC_Cell)
PP_Cell_fit = Plot_Fun2(k_PP_Cell)

plot_data_fit2(PC_PP_percent,PC_PP_B,fit_per,PC_PP_fit,PC_PP_sigma,"PC_PP")
plot_data_fit2(PP_Cell_percent,PP_Cell_B,fit_per,PP_Cell_fit,PP_Cell_sigma,"PP_Cell")
plot_data_fit2(PC_Cell_percent,PC_Cell_B,fit_per,PC_Cell_fit,PC_Cell_sigma,"PC_Cell")

def Plot_Fun3(tspan,ks,frac):
    ks = [item for subl in ks for item in subl]
    t0 = 0
    tf = tspan[-1]
    
    
    def ode3(x, tspan, ks):  
        x_S1 = x[0]
        x_S2 = x[1]
        x_B_S1 = x[2];
        x_B_S2 = x[3];
        x_B_S12 = x[4]
        x_Aq_S1 = x[5];
        x_Aq_S2 = x[6];
        x_Aq_S12 = x[7]


        
        k1 = ks[0]
        k2 = ks[1]
        k3 = ks[2]
        k4 = ks[3]
        k5 = ks[4]
        k6 = ks[5]
        k7 = ks[6]
        k8 = ks[7]
        k9 = ks[8]
        k10 = ks[9]
        k11 = ks[10]
        k12 = ks[11]
        k13 = ks[12]
        k14 = ks[13]
        k15 = ks[14]
        k16 = ks[15]
        k17 = ks[16]
        k18 = ks[17]
        k19 = ks[18]
        k20 = ks[19]
        k21 = ks[20]
        k22 = ks[21]
        a = ks[22]
        b = ks[23]
        y = ks[24]
        z = ks[25]
        
        k23 = ks[26]
        k24 = ks[27]
        k25 = ks[28]
        k26 = ks[29]
        k27 = ks[30]
        k28 = ks[31]
        k29 = ks[32]
        k30 = ks[33]
        k31 = ks[34]
        
        
        r1  = k1*x_S1        # S1  -> B1
        r2  = k2*x_B_S1      # B1  -> S1
        r3  = k3*x_S1        # S1  -> Aq1
        r4  = k4*x_Aq_S1     # Aq1 -> S1
        r5  = k5*x_Aq_S1     # Aq1 -> B1
        r6  = k6*x_B_S1      # B1  -> Aq1
        r7  = k7*x_B_S1      # B1  -> G1
        r8  = k8*x_Aq_S1     # Aq1 -> G1
        r9  = k9*x_S1        # S1  -> G1
        r10 = k10*x_B_S1     # B1  -> S_new1
        r11 = k11*x_Aq_S1    # Aq1 -> S_new1
        r12 = k12*x_S2       # S2  -> B2
        r13 = k13*x_B_S2     # B2  -> S2
        r14 = k14*x_S2       # S2  -> Aq2
        r15 = k15*x_Aq_S2    # Aq2 -> S2
        r16 = k16*x_Aq_S2    # Aq2 -> B2
        r17 = k17*x_B_S2     # B2  -> Aq2
        r18 = k18*x_B_S2     # B2  -> G2
        r19 = k19*x_Aq_S2    # Aq2 -> G2
        r20 = k20*x_S2       # S2  -> G2
        r21 = k21*x_B_S2     # B2  -> S_new2
        r22 = k22*x_Aq_S2    # Aq2 -> S_new2
        r25 = k25*x_Aq_S12   # Aq12  -> B12
        r26 = k26*x_B_S12    # B12   -> Aq12
        r27 = k27*x_B_S12    # B12   -> G12
        r28 = k28*x_Aq_S12   # Aq12  -> G12
        r30 = k30*x_B_S12    # B12   -> S_new12
        r31 = k31*x_Aq_S12   # Aq12  -> S_new12

        r24 = coupled(k24,1,1,x_S1,x_S2) # S1+S2 -> Aq12
        r29 = coupled(k29,1,1,x_S1,x_S2) # S1+S2 -> G12
        r23 = coupled(k23,a,b,x_S1,x_S2) # S1 + S2 -> B12
        

        sol = [-r1+r2-r3+r4-r9-(y*r23)-r24-r29, -r12+r13-r14+r15-r20-(z*r23)-r24-r29, r1+r5-r2-r6-r7-r10, r12-r13+r16-r17-r18-r21, r23+r25-r26-r27-r30, r3+r6-r4-r5-r8-r11, r14+r17-r15-r16-r19-r22, r24-r25+r26-r28-r31, r7+r8+r9,r18+r19+r20,r27+r28+r29,r10+r11,r21+r22,r30+r31]
        # S1 , S2, B1, B2, B12, Aq1, Aq2, Aq12, G1, G2, G12, S_new1, S_new2, S_new12  
        return sol
    
    
    
    def my_ls_func3(tspan,ks):
        """definition of function for LS fit
            x gives evaluation points,
            teta is an array of parameters to be varied for fit"""
        # create an alias to f which passes the optional params    
        f2 = lambda t,y: ode3(y, t ,ks)
        # calculate ode solution, retuen values for each entry of "x"
    
        r = solve_ivp(f2,(t0,tf),data0, t_eval=tspan,method ='Radau')
        #in this case, we only need one of the dependent variable values
        return np.transpose(r.y)
    

    #solve the system - the solution is in variable c
    # l =1 

    data0 = [100*frac, 100*(1-frac), 0, 0, 0,0,0,0,0,0,0,0,0,0] #inital conditions for ODEs
    fit = my_ls_func3(tspan,ks)
    fit[:,0] = np.add(np.add(np.add(np.add(fit[:,0],fit[:,1]),fit[:,11]),fit[:,12]),fit[:,13])
    fit[:,1] = np.add(np.add(fit[:,2],fit[:,3]),fit[:,4])
    fit[:,2] = np.add(np.add(fit[:,5],fit[:,6]),fit[:,7])
    fit[:,3] = np.add(np.add(fit[:,8],fit[:,9]),fit[:,10])
    fit = fit[:,:4]

    return fit
PP_PC_full = [PP_kf,PC_kf,PP_PC_power_n,PP_PC_kf]
PC_Cell_full = [PC_kf,Cell_kf,PC_Cell_power_n,PC_Cell_kf]
PP_Cell_full = [PP_kf,Cell_kf,PP_Cell_power_n,PP_Cell_kf]
PP_PC_fit= Plot_Fun3(tfit,PP_PC_full,.5)
PC_Cell_fit = Plot_Fun3(tfit,PC_Cell_full,.8)
PP_Cell_fit = Plot_Fun3(tfit,PP_Cell_full,.5)


plot_data_fit3(t,PP_PC_data0,tfit,PP_PC_fit,PP_PC_sigma0,"PP PC")
plot_data_fit3(t,PC_Cell_data0,tfit,PC_Cell_fit,PC_Cell_sigma0,"PC Cell")
plot_data_fit3(t,PP_Cell_data0,tfit,PP_Cell_fit,PP_Cell_sigma0,"PP Cell")


def Fit_Fun4(data,ks):
    t0 = t[0]
    tf = t[-1]
        
    
    def ode4(x, tspanc):    
        x_S1 = x[0]
        x_S2 = x[1]
        x_S3 = x[2]
        x_B_S1 = x[3]
        x_B_S2 = x[4]
        x_B_S3 = x[5]
        x_B_S12 = x[6]
        x_B_S13 = x[7]
        x_B_S23 = x[8]
        x_Aq_S1 = x[9]
        x_Aq_S2 = x[10]
        x_Aq_S3 = x[11]
        x_Aq_S12 = x[12]
        x_Aq_S13 = x[13]
        x_Aq_S23 = x[14]


        
        k1 = ks[0]
        k2 = ks[1]
        k3 = ks[2]
        k4 = ks[3]
        k5 = ks[4]
        k6 = ks[5]
        k7 = ks[6]
        k8 = ks[7]
        k9 = ks[8]
        k10 = ks[9]
        k11 = ks[10]
        k12 = ks[11]
        k13 = ks[12]
        k14 = ks[13]
        k15 = ks[14]
        k16 = ks[15]
        k17 = ks[16]
        k18 = ks[17]
        k19 = ks[18]
        k20 = ks[19]
        k21 = ks[20]
        k22 = ks[21]
        k23 = ks[22]
        k24 = ks[23]
        k25 = ks[24]
        k26 = ks[25]
        k27 = ks[26]
        k28 = ks[27]
        k29 = ks[28]
        k30 = ks[29]
        k31 = ks[30]
        k32 = ks[31]
        k33 = ks[32]
        a = ks[33]
        b = ks[34]
        y = ks[35]
        z = ks[36]
        k34 = ks[37]
        k35 = ks[38]
        k36 = ks[39]
        k37 = ks[40]
        k38 = ks[41]
        k39 = ks[42]
        k40 = ks[43]
        k41 = ks[44]
        k42 = ks[45]
        c = ks[46]
        d = ks[47]
        w = ks[48]
        v = ks[49]
        k43 = ks[50]
        k44 = ks[51]
        k45 = ks[52]
        k46 = ks[53]
        k47 = ks[54]
        k48 = ks[55]
        k49 = ks[56]
        k50 = ks[57]
        k51 = ks[58]
        e = ks[58]
        f = ks[59]
        u = ks[60]
        q = ks[61]
        k52 = ks[62]
        k53 = ks[63]
        k54 = ks[64]
        k55 = ks[65]
        k56 = ks[66]
        k57 = ks[67]
        k58 = ks[68]
        k59 = ks[69]
        k60 = ks[70]
        
        
        r1  = k1*x_S1        # S1    -> B1
        r2  = k2*x_B_S1      # B1    -> S1
        r3  = k3*x_S1        # S1    -> Aq1
        r4  = k4*x_Aq_S1     # Aq1   -> S1
        r5  = k5*x_Aq_S1     # Aq1   -> B1
        r6  = k6*x_B_S1      # B1    -> Aq1
        r7  = k7*x_B_S1      # B1    -> G1
        r8  = k8*x_Aq_S1     # Aq1   -> G1
        r9  = k9*x_S1        # S1    -> G1
        r10 = k10*x_B_S1     # B1    -> S_new1
        r11 = k11*x_Aq_S1    # Aq1   -> S_new1
        r12 = k12*x_S2       # S2    -> B2
        r13 = k13*x_B_S2     # B2    -> S2
        r14 = k14*x_S2       # S2    -> Aq2
        r15 = k15*x_Aq_S2    # Aq2   -> S2
        r16 = k16*x_Aq_S2    # Aq2   -> B2
        r17 = k17*x_B_S2     # B2    -> Aq2
        r18 = k18*x_B_S2     # B2    -> G2
        r19 = k19*x_Aq_S2    # Aq2   -> G2
        r20 = k20*x_S2       # S2    -> G2
        r21 = k21*x_B_S2     # B2    -> S_new2
        r22 = k22*x_Aq_S2    # Aq2   -> S_new2
        r23 = k23*x_S3       # S3    -> B3
        r24 = k24*x_B_S3     # B3    -> S3
        r25 = k25*x_S3       # S3    -> Aq3
        r26 = k26*x_Aq_S3    # Aq3   -> S3
        r27 = k27*x_Aq_S3    # Aq2   -> B3
        r28 = k28*x_B_S3     # B3    -> Aq3
        r29 = k29*x_B_S3     # B3    -> G3
        r30 = k30*x_Aq_S3    # Aq3   -> G3
        r31 = k31*x_S3       # S3    -> G3
        r32 = k32*x_B_S3     # B3    -> S_new3
        r33 = k33*x_Aq_S3    # Aq3   -> S_new3
        r36 = k36*x_Aq_S12   # Aq12  -> B12
        r37 = k37*x_B_S12    # B12   -> Aq12
        r38 = k38*x_B_S12    # B12   -> G12
        r39 = k39*x_Aq_S12   # Aq12  -> G12
        r41 = k41*x_B_S12    # B12   -> S_new12
        r42 = k42*x_Aq_S12   # Aq12  -> S_new12
        r45 = k45*x_Aq_S13   # Aq13  -> B13
        r46 = k46*x_B_S13    # B13   -> Aq13
        r47 = k47*x_B_S13    # B13   -> G13
        r48 = k48*x_Aq_S13   # Aq13  -> G13
        r50 = k50*x_B_S13    # B13   -> S_new13
        r51 = k51*x_Aq_S13   # Aq13  -> S_new13
        r54 = k54*x_Aq_S23   # Aq23  -> B23
        r55 = k55*x_B_S23    # B23   -> Aq23
        r56 = k56*x_B_S23    # B23   -> G23
        r57 = k57*x_Aq_S23   # Aq23  -> G23
        r59 = k59*x_B_S23    # B23   -> S_new23
        r60 = k60*x_Aq_S23   # Aq23  -> S_new23
        
        
        r34 = coupled(k34,a,b,x_S1,x_S2) # S1+S2 -> B12
        r35 = coupled(k35,1,1,x_S1,x_S2) # S1+S2 -> Aq12
        r40 = coupled(k40,1,1,x_S1,x_S2) # S1+S2 -> G12
        r43 = coupled(k43,c,d,x_S1,x_S3) # S1+S3 -> B13
        r44 = coupled(k44,1,1,x_S1,x_S3) # S1+S3 -> Aq13
        r49 = coupled(k49,1,1,x_S1,x_S3) # S1+S3 -> G13
        r53 = coupled(k53,1,1,x_S2,x_S3) # S2+S3 -> Aq23
        r58 = coupled(k58,1,1,x_S2,x_S3) # S2+S3 -> G23
        r52 = coupled(k52,e,f,x_S2,x_S3) # S2+S3 -> B23
        

        
        sol = [-r1-r3-r9-(y*r34)-r35-r40-(w*r43)-r44-r49+r2+r4, -r12-r14-r20-(z*r34)-r35-r40-(u*r52)-r53-r58+r13+r15, -r23-r25-r31-(v*r43)-r44-r49-(q*r52)-r53-r58+r24+r26 , r1+r5-r2-r6-r7-r10, r12+r16-r13-r17-r18-r21, r23+r27-r24-r28-r29-r32, r34+r36-r37-r38-r41, r43+r45-r46-r47-r50, r52+r54-r55-r56-r59, r3+r6-r4-r5-r8-r11, r14+r17-r15-r17-r19-r22, r25+r28-r26-r27-r30-r33, r35+r37-r36-r39-r42, r44+r46-r45-r48-r51, r53+r55-r54-r57-r60, r7+r8+r9, r18+r19+r20, r29+r30+r31, r38+r39+r40, r47+r48+r49, r56+r57+r58, r10+r11, r21+r22, r32+r33, r41+r42, r50+r51, r59+r60]
        # S1 , S2,, S3, B1, B2, B3, B12, B13, B23, Aq1, Aq2, Aq3, Aq12, Aq13, Aq23, G1, G2, G3, G12, G13, G23, S_new1, S_new2, S_new3, S_new12, S_new13, S_new23  
        return sol
    
    
    
    def my_ls_func4(tspan,ks):
        """definition of function for LS fit
            x gives evaluation points,
            teta is an array of parameters to be varied for fit"""
        # create an alias to f which passes the optional params    
        f2 = lambda t,y: ode4(y, t)
        # calculate ode solution, retuen values for each entry of "x"
    
        r = solve_ivp(f2,(t0,tf),data0, t_eval=tfit,method ='Radau')
        #in this case, we only need one of the dependent variable values
        return np.transpose(r.y)

    data0 = [100/3, 100/3, 100/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] #inital conditions for ODEs


    
    tfit = np.linspace(0,t[-1])
    fit = my_ls_func4(tfit,ks)
    fit[:,0] = np.add(np.add(np.add(np.add(np.add(np.add(np.add(np.add(fit[:,0],fit[:,1]),fit[:,2]),fit[:,21]),fit[:,22]),fit[:,23]),fit[:,24]),fit[:,25]),fit[:,26])
    fit[:,1] = np.add(np.add(np.add(np.add(np.add(fit[:,3],fit[:,4]),fit[:,5]),fit[:,6]),fit[:,7]),fit[:,8])
    fit[:,2] = np.add(np.add(np.add(np.add(np.add(fit[:,9],fit[:,10]),fit[:,11]),fit[:,12]),fit[:,13]),fit[:,14])
    fit[:,3] = np.add(np.add(np.add(np.add(np.add(fit[:,15],fit[:,16]),fit[:,17]),fit[:,18]),fit[:,19]),fit[:,20])
    fit = fit[:,:4]
    data = np.array(data)
    fit = np.array(fit)
    return [t,data,tfit,fit]








k_PP_PC_Cell = [PP_kf,PC_kf,Cell_kf,PP_PC_power_n,PP_PC_kf,PP_Cell_power_n,PP_Cell_kf,PC_Cell_power_n,PC_Cell_kf]
sd_PP_PC_Cell = [PP_sd,PC_sd,Cell_sd,PP_PC_power_sd,PP_PC_sd,PP_Cell_power_sd,PP_Cell_sd,PC_Cell_power_sd,PC_Cell_sd]
k_PP_PC_Cell = [item for subl in k_PP_PC_Cell for item in subl]
sd_PP_PC_Cell = [item for subl in sd_PP_PC_Cell for item in subl]

[PP_PC_Cell_t,PP_PC_Cell_data,PP_PC_Cell_tfit,PP_PC_Cell_fit] = Fit_Fun4(PP_PC_Cell_data0,k_PP_PC_Cell)

plot_data_fit3(PP_PC_Cell_t,PP_PC_Cell_data,PP_PC_Cell_tfit,PP_PC_Cell_fit,PP_PC_Cell_sigma0,"PP_PC_Cell")



out = pd.ExcelWriter("Output.xlsx")
Sheet = "Full"
k_PP_PC_Cell = pd.DataFrame(data=k_PP_PC_Cell)
sd_PP_PC_Cell = pd.DataFrame(data=sd_PP_PC_Cell)
PP_tf = pd.DataFrame(data=t)
PP_data = pd.DataFrame(data=PP_data0)
PP_sigma0 = pd.DataFrame(data=PP_sigma0)
PP_tfitf = pd.DataFrame(data=tfit)
PP_fitf = pd.DataFrame(data=PP_fitf)
k_PP_PC_Cell.to_excel(out,sheet_name=Sheet,startcol=0,startrow =1,index=False,header=False)
sd_PP_PC_Cell.to_excel(out,sheet_name=Sheet,startcol=1,startrow =1,index=False,header=False)
PP_tf.to_excel(out,sheet_name=Sheet,startcol=2,startrow =1,index=False,header=False)
PP_data.to_excel(out,sheet_name=Sheet,startcol=3,startrow =1,index=False,header=False)
PP_tf.to_excel(out,sheet_name=Sheet,startcol=2,startrow =8,index=False,header=False)
PP_sigma0.to_excel(out,sheet_name=Sheet,startcol=3,startrow =8,index=False,header=False)
PP_tfitf.to_excel(out,sheet_name=Sheet,startcol=7,startrow =1,index=False,header=False)
PP_fitf.to_excel(out,sheet_name=Sheet,startcol=8,startrow =1,index=False,header=False)


i =1
PC_tf = pd.DataFrame(data=t)
PC_data = pd.DataFrame(data=PC_data0)
PC_sigma0 = pd.DataFrame(data=PC_sigma0)
PC_tfitf = pd.DataFrame(data=tfit)
PC_fitf = pd.DataFrame(data=PC_fitf)
PC_tf.to_excel(out,sheet_name=Sheet,startcol=0+(13*i),startrow =1,index=False,header=False)
PC_data.to_excel(out,sheet_name=Sheet,startcol=1+(13*i),startrow =1,index=False,header=False)
PC_tf.to_excel(out,sheet_name=Sheet,startcol=0+(13*i),startrow =8,index=False,header=False)
PC_sigma0.to_excel(out,sheet_name=Sheet,startcol=1+(13*i),startrow =8,index=False,header=False)
PC_tfitf.to_excel(out,sheet_name=Sheet,startcol=5+(13*i),startrow =1,index=False,header=False)
PC_fitf.to_excel(out,sheet_name=Sheet,startcol=6+(13*i),startrow =1,index=False,header=False)

i =2
Cell_tf = pd.DataFrame(data=t)
Cell_data = pd.DataFrame(data=Cell_data0)
Cell_sigma0 = pd.DataFrame(data=Cell_sigma0)
Cell_tfitf = pd.DataFrame(data=tfit)
Cell_fitf = pd.DataFrame(data=Cell_fitf)
Cell_tf.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =1,index=False,header=False)
Cell_data.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =1,index=False,header=False)
Cell_tf.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =8,index=False,header=False)
Cell_sigma0.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =8,index=False,header=False)
Cell_tfitf.to_excel(out,sheet_name=Sheet,startcol=5+(12*i),startrow =1,index=False,header=False)
Cell_fitf.to_excel(out,sheet_name=Sheet,startcol=6+(12*i),startrow =1,index=False,header=False)


i =3
PP_PC_tf = pd.DataFrame(data=t)
PP_PC_data = pd.DataFrame(data=PP_PC_data0)
PP_PC_sigma0 = pd.DataFrame(data=PP_PC_sigma0)
PP_PC_tfitf = pd.DataFrame(data=tfit)
PP_PC_fitf = pd.DataFrame(data=PP_PC_fit)
PP_PC_tf.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =1,index=False,header=False)
PP_PC_data.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =1,index=False,header=False)
PP_PC_tf.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =8,index=False,header=False)
PP_PC_sigma0.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =8,index=False,header=False)
PP_PC_tfitf.to_excel(out,sheet_name=Sheet,startcol=5+(12*i),startrow =1,index=False,header=False)
PP_PC_fitf.to_excel(out,sheet_name=Sheet,startcol=6+(12*i),startrow =1,index=False,header=False)


i =4
PP_Cell_tf = pd.DataFrame(data=t)
PP_Cell_data = pd.DataFrame(data=PP_Cell_data0)
PP_Cell_sigma0 = pd.DataFrame(data=PP_Cell_sigma0)
PP_Cell_tfitf = pd.DataFrame(data=tfit)
PP_Cell_fitf = pd.DataFrame(data=PP_Cell_fit)
PP_Cell_tf.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =1,index=False,header=False)
PP_Cell_data.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =1,index=False,header=False)
PP_Cell_tf.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =8,index=False,header=False)
PP_Cell_sigma0.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =8,index=False,header=False)
PP_Cell_tfitf.to_excel(out,sheet_name=Sheet,startcol=5+(12*i),startrow =1,index=False,header=False)
PP_Cell_fitf.to_excel(out,sheet_name=Sheet,startcol=6+(12*i),startrow =1,index=False,header=False)


i =5
PC_Cell_tf = pd.DataFrame(data=t)
PC_Cell_data = pd.DataFrame(data=PC_Cell_data0)
PC_Cell_sigma0 = pd.DataFrame(data=PC_Cell_sigma0)
PC_Cell_tfitf = pd.DataFrame(data=tfit)
PC_Cell_fitf = pd.DataFrame(data=PC_Cell_fit)
PC_Cell_tf.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =1,index=False,header=False)
PC_Cell_data.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =1,index=False,header=False)
PC_Cell_tf.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =8,index=False,header=False)
PC_Cell_sigma0.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =8,index=False,header=False)
PC_Cell_tfitf.to_excel(out,sheet_name=Sheet,startcol=5+(12*i),startrow =1,index=False,header=False)
PC_Cell_fitf.to_excel(out,sheet_name=Sheet,startcol=6+(12*i),startrow =1,index=False,header=False)


i =6
PP_PC_Cell_data = pd.DataFrame(data=PP_PC_Cell_data)
PP_PC_Cell_sigma0 = pd.DataFrame(data=PP_PC_Cell_sigma0)
PP_PC_Cell_tfit = pd.DataFrame(data=PP_PC_Cell_tfit)
PP_PC_Cell_t = pd.DataFrame(data=PP_PC_Cell_t)
PP_PC_Cell_fit = pd.DataFrame(data=PP_PC_Cell_fit)
PP_PC_Cell_t.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =1,index=False,header=False)
PP_PC_Cell_data.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =1,index=False,header=False)
PP_PC_Cell_t.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =8,index=False,header=False)
PP_PC_Cell_sigma0.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =8,index=False,header=False)
PP_PC_Cell_tfit.to_excel(out,sheet_name=Sheet,startcol=5+(12*i),startrow =1,index=False,header=False)
PP_PC_Cell_fit.to_excel(out,sheet_name=Sheet,startcol=6+(12*i),startrow =1,index=False,header=False)


k_PC_PP = [PC_kf,PP_kf,PC_PP_power]
k_PC_Cell = [PC_kf,Cell_kf,PC_Cell_power]
k_PP_Cell = [PP_kf,Cell_kf,PP_Cell_power]
PC_PP_fit = Plot_Fun2(k_PC_PP)
PC_Cell_fit = Plot_Fun2(k_PC_Cell)
PP_Cell_fit = Plot_Fun2(k_PP_Cell)

i = 7
PC_PP_percent = pd.DataFrame(data=PC_PP_percent)
PC_PP_B = pd.DataFrame(data=PC_PP_B)
PC_PP_per = pd.DataFrame(data=fit_per)
PC_PP_fit = pd.DataFrame(data=PC_PP_fit)
PC_PP_percent.to_excel(out,sheet_name=Sheet,startcol=0+(12*i),startrow =1,index=False,header=False)
PC_PP_B.to_excel(out,sheet_name=Sheet,startcol=1+(12*i),startrow =1,index=False,header=False)
PC_PP_per.to_excel(out,sheet_name=Sheet,startcol=2+(12*i),startrow =1,index=False,header=False)
PC_PP_fit.to_excel(out,sheet_name=Sheet,startcol=3+(12*i),startrow =1,index=False,header=False)


PC_Cell_percent = pd.DataFrame(data=PC_Cell_percent)
PC_Cell_B = pd.DataFrame(data=PC_Cell_B)
PC_Cell_per = pd.DataFrame(data=fit_per)
PC_Cell_fit = pd.DataFrame(data=PC_Cell_fit)
PC_Cell_percent.to_excel(out,sheet_name=Sheet,startcol=0+(12*i) + 5 ,startrow =1,index=False,header=False)
PC_Cell_B.to_excel(out,sheet_name=Sheet,startcol=1+(12*i) + 5 ,startrow =1,index=False,header=False)
PC_Cell_per.to_excel(out,sheet_name=Sheet,startcol=2+(12*i) + 5 ,startrow =1,index=False,header=False)
PC_Cell_fit.to_excel(out,sheet_name=Sheet,startcol=3+(12*i) + 5 ,startrow =1,index=False,header=False)


PP_Cell_percent = pd.DataFrame(data=PP_Cell_percent)
PP_Cell_B = pd.DataFrame(data=PP_Cell_B)
PP_Cell_per = pd.DataFrame(data=fit_per)
PP_Cell_fit = pd.DataFrame(data=PP_Cell_fit)
PP_Cell_percent.to_excel(out,sheet_name=Sheet,startcol=0+(12*i) + 10 ,startrow =1,index=False,header=False)
PP_Cell_B.to_excel(out,sheet_name=Sheet,startcol=1+(12*i) + 10 ,startrow =1,index=False,header=False)
PP_Cell_per.to_excel(out,sheet_name=Sheet,startcol=2+(12*i) + 10 ,startrow =1,index=False,header=False)
PP_Cell_fit.to_excel(out,sheet_name=Sheet,startcol=3+(12*i) + 10 ,startrow =1,index=False,header=False)

out.save()
