import matplotlib.pyplot as plt
import math
import random
from scipy.stats import norm
import pandas as pd

##
# 程序由三部分组成：1）似然值计算；2）动力学模拟；3）MCMC参数优化
# 1）跟2）由定义函数实现
# 请大家认真阅读三个需要注意的地方
##

## likelihood
def lik(sim,obs):
    lk = 0
    for i in range(len(sim)):
        # 注意1：因为流行病学模型变量变化幅度较大，因此这里没有用单一方差（第三个参数），而是用平均值（第二个）
        # 要不然计算出来的似然值会很大而导致计算错误（比如溢出）；且为了避免出现方差为0（会报错），方差加1
        Tlk = norm.logpdf(obs[i],sim[i],sim[i]+1)
        lk += Tlk
    return lk

##  simulation
def sim(n1):
    #   Parameters
    t0 = 2008
    dt = 1 / 52
    mu = 0.02
    rr = 0.5

    # initial population
    pop = 5000000
    s0 = 26 / 400
    i0 = 0.002
    r0 = 1
    S0 = pop * s0 / (s0 + i0 + r0)
    I0 = pop * i0 / (s0 + i0 + r0)
    R0 = pop * r0 / (s0 + i0 + r0)
    H0 = 0

    # state variable/list
    S = [S0]
    I = [I0]
    R = [R0]
    T = [t0]
    H = [H0]
    Beta = [180 * (5 + 2 * math.sin(math.pi * t0 + 5))]

    # update state variables each day
    for i in range((2020 - 2008) * 52):
        S0 = S[-1]
        I0 = I[-1]
        R0 = R[-1]
        T0 = T[-1]
        pop = S0 + I0 + R0
        beta = 180 * (5 + 2 * math.sin(math.pi * (T0 + dt) + 5))
        S1 = S0 + (mu * pop - beta * S0 * I0 / pop - mu * S0) * dt
        I1 = I0 + (beta * S0 * I0 / pop - n1 * I0 - mu * I0) * dt
        R1 = R0 + (n1 * I0 - mu * R0) * dt
        T1 = T0 + dt

        # 注意2：报告病例并不等于I，也不等于I的增加量，而只是新增病例，这里是(beta * S0 * I0 / pop) * dt
        # 另外，我们以周为单位模拟，如果模型以天为单位模拟，需要把一个周每天的新增病例加起来形成一周的新增病例
        # 然后再跟观测的一周为单位的数据比较；这里采取最简单的方式

        # 如果不想在病例报告这里引入随机性，用下面代码就行
        # H1 = (beta * S0 * I0 / pop) * dt * rr

        # 在病例报告中引入随机性的代码如下；这里为了减少随机性，方差取了0.1
        Htem = (beta * S0 * I0 / pop) * dt * rr
        H1 = int(norm.rvs(Htem,1,1)[0])

        # 因为采用基于均值/方差的正态分布采样，可能出现负数，真实情况不可能出现，所以特殊处理为0
        if H1 <= 0 :
            H1 = 0
        Beta.append(beta)
        S.append(S1)
        I.append(I1)
        R.append(R1)
        T.append(T1)
        H.append(H1)
    return H

## 请指向正确的数据源
fn = r'cases_obs_final.csv'
fig = r'gammaNew'

## observed data
df = pd.read_csv(fn)
Tlist = list(df.loc[:,'time'])
Dlist = list(df.loc[:,'cases'])

## MCMC
# 注意3：如果你从偏离22比较大的数值开始，程序基本不能跑下来，如果想探索这些参数，就需要基于一步模拟评估的方式
# 另外，为了更快找到最优参数，还可以考虑平行MCMC，即用并行程序一次性同时运行多个MCMC链
# 但作为练习，gamma可以从21开始，仅仅体会一下目前系统的复杂性，以后如果用到由基础了就可以很好解决更复杂情况

N = 5000
gamma = 21
Gamma = [gamma]
sim0 = sim(Gamma[-1]) #模拟
lk0 = lik(sim0,Dlist) #计算似然值
LK = [lk0]
Continue=0
E=0.0001

for i in range(N):
    # 参数随机游走；我采用了0.01的方差
    gammai = norm.rvs(Gamma[-1],0.01,1)[0]

    # simulation
    simi = sim(gammai)

    # likelihood calculation
    lki = lik(simi,Dlist)

    # odd ratio: R > 1 无条件接受； R < 1，以R概率接受
    R = math.exp(lki-LK[-1])
    # R = math.exp(LK[-1]-lki)
    if random.random() < R:
        if abs(gamma-gammai)<E:
            Continue=Continue+1
        else:
            Continue =0
        gamma = gammai
        Gamma.append(gamma)
        LK.append(lki)
        if Continue>=10:
            break
        print(i,gamma,lki)
    else:
        print(i,gammai,lki,"Not accepted:",R)

print("Last accepted",Gamma[-1])
fn = r'cases_obs_final.csv'
## observed data
df = pd.read_csv(fn)
H=sim(Gamma[-1])
Dlist = list(df.loc[:,'cases'])
# 作图，画出拟合情况
plt.subplot(211)
plt.plot(Tlist,Dlist, 'r')
plt.plot(Tlist,H, 'b')
plt.legend(['real','estimate'])
plt.xlabel('Year')
plt.ylabel('Population')
plt.title("Gamma="+str(Gamma[-1]))
# 作图，画出Gamma的优化过程
plt.subplot(212)
plt.scatter(LK,Gamma,c='r',marker='.',s= 20,edgecolor='none')
# plt.plot(LK,Gamma,'r')
plt.savefig(fig)
plt.show()