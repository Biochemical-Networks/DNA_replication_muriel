import numpy as np


def P_Sr_step1BACKWARD(gtt):
    return 1/(4 *(3/4 - np.exp(gtt)/(2*(1 + np.exp(gtt)))))

def P_Sw_step1BACKWARD(gtt):
    return 1/(2 *(1 + np.exp(gtt))*(3/4 - np.exp(gtt)/(2 *(1 + np.exp(gtt)))))

def P_Sr_step1FORWARD(gtt):
    return 1/(2 *(1 + np.exp(-gtt))*(1 - 1/(2 *(1 + np.exp(-gtt))) - np.exp(-gtt)/(2* (1 + np.exp(-gtt)))))

def P_Sw_step1FORWARD(gtt):
    return  ((np.exp(-gtt))/(2 *(1 + np.exp(-gtt)) *(1 - 1/(2* (1 + np.exp(-gtt))) - np.exp(-gtt)/(
   2 *(1 + np.exp(-gtt))))))


L_gtt = []
L_Psw = []
L_Psr = []
L_error_F = []
L_error_B = []
for gtt in range(0,11,1):
    error_F= P_Sw_step1FORWARD(gtt)/(P_Sw_step1FORWARD(gtt) + P_Sr_step1FORWARD(gtt))
    error_B= P_Sw_step1BACKWARD(gtt)/(P_Sw_step1BACKWARD(gtt) + P_Sr_step1BACKWARD(gtt))
    L_gtt.append(gtt)
    L_error_F.append(error_F)
    L_error_B.append(error_B)

print("for g_tt is" , L_gtt)
print("Step 1 has forward discrimination: the error for every g_tt in the limit of g_pol going to infinity from 0 to 10 in steps of 1", L_error_F)
print("Step 1 has backward discrimination: the error for every g_tt in the limit of g_pol going to infinity from 0 to 10 in steps of 1", L_error_B)
