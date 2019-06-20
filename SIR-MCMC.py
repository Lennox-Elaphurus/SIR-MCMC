# 采用走完完整的一次再做评估的方法
import random
from scipy.stats import norm
import math

MAX_PACE=10000000
fig = r'simulation'
START_YEAR = 2008
END_YEAR=2020
dt = 1/52
gamma =0
lastGamma=20 # 22
Ratio=0
GAMMA=[]

sigma=10
E=0.001
Continue=0
MIN_CONTINUE=100

global mu
mu = 1/50
s0 = 26 / 400
i0 = 0.002
r0 = 1
global reportRate
reportRate=0.5
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
global Infe
Infe = [I0]
global R
R=[R0]
LK=[]


def get_beta(index):
    return float(180 * (5 + 2 * math.sin(math.pi * index* dt + 5)))


global BETA
BETA=[get_beta(0)]
global points
points=[]


def import_data(file_location):   # sir_case.csv
# 读取x跟y
    with open(file_location,"r") as file:
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


def estimate(this_gamma):  # directly write to global S,I,R
    global S0
    global I0
    global R0
    global S
    global Infe
    global R
    global BETA
    global pop
    global mu
    # state variable/list
    # print(type(R))
    S.clear()
    Infe.clear()
    R.clear()
    S.append(S0)
    Infe.append(I0)
    R.append(R0)
    # update state variables each day
    for i in range(1,int((END_YEAR - START_YEAR)/dt)+1):
        S0 = S[-1]
        I0 = Infe[-1]
        R0 = R[-1]
        pop = S0 + I0 + R0
        beta=get_beta(i)
        s1 = S0 + (mu * pop - beta * S0 * I0 / pop - mu * S0) * dt
        i1 = I0 + (beta * S0 * I0 / pop - this_gamma * I0 - mu * I0) * dt
        r1 = R0 + (this_gamma * I0 - mu * R0) * dt
        BETA.append(beta)
        S.append(s1)
        Infe.append(i1)
        R.append(r1)


def get_likelihood(this_sigma):
    global points
    global Infe
    global reportRate
    lk2=0
    for i2 in range(len(points)):
        lk2=lk2+float(norm.logpdf(points[i2][1], Infe[i2]*reportRate, this_sigma))
    return lk2


# MCMC
estimate(lastGamma)
lastLk=get_likelihood(sigma)
# print(lastLk)
for i in range(MAX_PACE):
    gamma=abs(float(norm.rvs(lastGamma, sigma, 1)))
    estimate(gamma)
    lk=get_likelihood(sigma)
    try:
        Ratio = math.exp(lastLk - lk)
    except OverflowError:
        # directly accept
        print("An OverflowError occurred.")
        print("lastLk:", lastLk, "lk:", lk)
        print("lastLk - lk", lastLk - lk)
        if lk<lastLk:
            Ratio=1
        else:
            Ratio=0
    if random.random() < Ratio:
        if len(GAMMA) > 1:
            lastGamma = GAMMA[-1]
            if abs(gamma - lastGamma) < E:
                Continue = Continue + 1
                if Continue % 10 == 0:  # 10 was set by hand
                    sigma = sigma / 5  # 5 was set by hand
                    # need to estimate again when sigma changed
                    estimate(lastGamma)
                    lastLk=get_likelihood(sigma)
                    print("Adjust sigma to ", sigma)
                if Continue > MIN_CONTINUE:
                    print("Stop by sufficiently accurate answer.")
                    print("Last accepted Ratio:", Ratio)
                    break
            else:
                Continue = 0

        lastGamma = gamma
        lastLk=lk
        GAMMA.append(lastGamma)
        print("Accepted Ratio:", Ratio)
        print("Accepted gamma:",lastGamma)

print("Final answer:",lastGamma)