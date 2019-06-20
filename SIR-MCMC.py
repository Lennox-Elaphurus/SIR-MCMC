# 采用走完完整的一次再做评估的方法
import random
from scipy.stats import norm
import math

fig = r'simulation'
START_YEAR = 2008
END_YEAR=2020
dt = 1/52
gamma = 22
global mu
mu = 1/50

s0 = 26 / 400
i0 = 0.002
r0 = 1
reportRate = 0.5
global pop
pop = 5000000  # initial population
s0 = 26 / 400
i0 = 0.002
r0 = 1
global S0
S0 = pop * s0 / (s0 + i0 + r0)
global I0
I0 = pop * i0 / (s0 + i0 + r0)
global R0
R0 = pop * r0 / (s0 + i0 + r0)
global S
S = [S0]
global I
I = [I0]
global R
R = [R0]

def get_beta(index):
    return float(180 * (5 + 2 * math.sin(math.pi * index* dt + 5)))

global BETA
BETA= get_beta(0)
global points
points=[]

def import_data(fileLocation):   # sir_case.csv
# 读取x跟y
    with open(fileLocation,"r") as file:
        global points
        # seperate data into items
        pointlist=file.read()
        pointlist=pointlist.replace(",", " ")
        pointlist = pointlist.split()
        temp=0
        flag=True
        for idx in range(2,len(pointlist)):
            if flag:
                temp=float(pointlist[idx])
                flag=False
            else:
                points.append((float(temp),float(pointlist[idx])))
                flag=True


import_data("sir_case.csv")
# for i in range(int((END_YEAR-START_YEAR)/dt)):
#     print(points[i][1])

def estimate(this_gamma): # directly write to global S,I,R
    global S0
    global I0
    global R0
    global S
    global I
    global R
    global BETA
    global pop
    global mu
    # state variable/list
    S = [S0]
    I = [I0]
    R = [R0]
    # update state variables each day
    for i in range(1,int((END_YEAR - START_YEAR)/dt)+1):
        S0 = S[-1]
        I0 = I[-1]
        R0 = R[-1]
        pop = S0 + I0 + R0
        beta=get_beta(i)
        s1 = S0 + (mu * pop - beta * S0 * I0 / pop - mu * S0) * dt
        i1 = I0 + (beta * S0 * I0 / pop - thisGamma * I0 - mu * I0) * dt
        r1 = R0 + (thisGamma * I0 - mu * R0) * dt
        BETA.append(beta)
        S.append(s1)
        I.append(i1)
        R.append(r1)


def get_likelihood(time,sigma):
    global points
    # I1 = I0 + (beta*S0*I0/pop - gamma*I0 - mu*I0)*dt
    if time>=1:
        estimate=points[time-1][1]+
        return float(norm.logpdf(points[time][1], estimate, sigma))
