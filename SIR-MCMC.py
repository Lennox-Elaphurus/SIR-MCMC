# 采用走完完整的一次再做评估的方法
import random
from scipy.stats import norm
import math
import matplotlib.pyplot as plt
import pandas

MAX_PACE=100000
fig = r'simulation'
START_YEAR = 2008
END_YEAR=2020
dt = 1/52
gamma = 0
lastGamma =70 # 22
Ratio=0
GAMMA=[]

global reportRate
reportRate=0.5
global sigma
sigma=5
sigma0=sigma
lastSigma=sigma
E=1   # 3我之前不能收敛是因为这个调得太小了
E0=E
lastE=E
Continue=0
MIN_CONTINUE=15

global mu
mu = 1/50
s0 = 26 / 400
i0 = 0.002
r0 = 1

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
LK=[0]
global ignore
ignore=False



def get_beta(index):
    return float(180 * (5 + 2 * math.sin(math.pi * index* dt + 5)))

global H
Htem = (get_beta(0) * S0 * I0 / pop) * dt * reportRate
H0 = float(norm.rvs(Htem,0.1,1)[0])
H=[H0]
global points
points=[]


def import_data(file_location):   # sir_case.csv
    df = pandas.read_csv(file_location)
    Tlist = list(df.loc[:, 'time'])
    Dlist = list(df.loc[:, 'cases'])
    global points
    for i in range(len(Tlist)):
        points.append((float(Tlist[i]),float(Dlist[i])))


def estimate(this_gamma):  # directly write to global S,I,R
    global S0
    global I0
    global R0
    global H0
    global S
    global Infe
    global R
    global BETA
    global pop
    global mu
    global dt
    global H
    # state variable/list
    # print(type(R))
    S.clear()
    Infe.clear()
    R.clear()
    S.append(S0)
    Infe.append(I0)
    R.append(R0)
    H.clear()
    H.append(H0)
    # update state variables each day
    for i in range(1,int((END_YEAR - START_YEAR)/dt)+1):
        s1 = S[-1]
        i1 = Infe[-1]
        r1 = R[-1]
        pop = S0 + I0 + R0
        beta = get_beta(i)
        s2 = s1 + (mu * pop - beta * s1 * i1 / pop - mu * s1) * dt
        i2 = i1 + (beta * s1 * i1 / pop - this_gamma * i1 - mu * i1) * dt
        r2 = r1 + (this_gamma * i1 - mu * r1) * dt
        hTem = (beta * s1 * i1 / pop) * dt * reportRate
        h2 = float(norm.rvs(hTem, 0.1, 1)[0])
        # h2 = float(norm.rvs(i2, 0.1, 1)) # 对预测值随机扰动
        H.append(h2)
        S.append(s2)
        Infe.append(i2)
        R.append(r2)


def get_likelihood(this_sigma):
    global points
    global H
    global reportRate
    global ignore
    lk2=0
    for i2 in range(len(points)):
        try:
            # if points[i2][1] Infe[i2]*reportRate
            # lk2=lk2+float(norm.logpdf(points[i2][1], Infe[i2]*reportRate, 1000))
            lk2 = lk2+abs(points[i2][1]-H[i2]*reportRate)**5
        except OverflowError:
            # print("An OverflowError occurred in abs**5.")
            ignore = True
            break
    lk2=lk2*1e-21  # just to decrease the lk
    # print("lk2",lk2)
    return lk2


def draw():
    global points
    global H
    global gamma
    global sigma
    global reportRate
    global lk
    global fig
    global reportRate
    global E
    global LK
    global GAMMA
    t=[]
    real_i=[]
    for i in range(len(points)):
        t.append(points[i][0])
        real_i.append(points[i][1])
        # H[i]=H[i]*reportRate
    plt.ion()# 绘图或者从磁盘读取图像并进行图像处理操作
    H[0]=0
    # plt.scatter(t,real_i,c='b',marker='.',s= 10,edgecolor='none')
    # plt.scatter(t,H,c='r',marker='.',s= 20,edgecolor='none')
    plt.plot(t,real_i, 'r')
    plt.plot(t,H, 'b')
    plt.legend(['real','estimate'])
    plt.savefig(fig)
    plt.xlabel('Year')
    plt.ylabel('Population')
    plt.suptitle("Gamma="+str(gamma)+" Report Rate="+str(reportRate))
    plt.title("Sigma="+str(sigma)+" E="+str(E))
    plt.show()
    plt.pause(5)
    plt.savefig(fig)
    plt.close()
    plt.scatter(GAMMA,LK,s=1)
    plt.show()
    plt.pause(5)
    plt.savefig("gamma-lk")



import_data("cases_obs_final.csv")
# MCMC
isMCMC=True
estimate(lastGamma)
lastLk=get_likelihood(sigma)
GAMMA.append(lastGamma)
print("First gamma:",lastGamma,"First lk:",lastLk)
# print(lastLk)
for cnt_step in range(MAX_PACE):
    gamma=float(norm.rvs(lastGamma, sigma, 1))
    ignore=False
    estimate(gamma)
    lk=get_likelihood(sigma)
    try:
        # Ratio = math.exp(lastLk - lk)
        Ratio=lastLk/lk
    except OverflowError:
        # directly accept
        # print("An OverflowError occurred in exp.")
        # print("lastLk:", lastLk, "lk:", lk)
        # print("lastLk - lk", lastLk - lk)
        if lk<lastLk:
            Ratio=1
        else:
            Ratio=0
    if ignore is False and random.random() < Ratio:
        if isMCMC is False:
            if Ratio<1:
                continue
        if abs(gamma - lastGamma) < E:
            Continue = Continue + 1
            if Continue % 5 == 0:  # 2 was set by hand
                divide=5
                if sigma < 1:
                    divide=2
                lastE=E
                E=E/divide
                lastSigma=sigma
                sigma = sigma /divide  # 5 was set by hand
                if lastSigma>1:
                    lastSigma=1
                    sigma=1
                # if 1 < sigma <= 2:
                #     lastSigma = sigma
                #     sigma = sigma / 1.1
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
            if sigma!= lastSigma:
                print("Reset sigma to",lastSigma)
            sigma=lastSigma
            E=lastE

        lastGamma = gamma
        lastLk=lk
        GAMMA.append(lastGamma)
        LK.append(lastLk)
        print("Accepted gamma:",lastGamma,end=" ")
        print("Accepted Ratio:", Ratio)
        # print("lk:",lk)
        # if not sigma==sigma0:
        #     print("Sigma=",sigma," E=",E)
        # draw()
    if cnt_step % 1000 == 0:
        print("count step:", cnt_step)
        # isMCMC = False
        # draw()

print("Final answer:",lastGamma)
draw()
