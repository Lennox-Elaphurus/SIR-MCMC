import matplotlib.pyplot as plt
import math
import numpy
rate=1/50
S_0=26/400
I_0=0.0002
R_0=1
popsize=5E06
gamma=22
N=S_0+I_0+R_0
paceOfTime=1/365
START_YEAR=2008
END_YEAR=2020
S=[]
S.append(S_0*popsize)
I=[]
I.append(I_0*popsize)
R=[]
R.append(R_0*popsize)
X=[]
X.append(START_YEAR)
def beta(t):
    return 180*(5+2*math.sin(math.pi*(t)+5))


time=1
# run nian? the effect is too tiny
while time*paceOfTime < END_YEAR-START_YEAR+1: # time =t*365
    N=S[-1]+I[-1]+R[-1]
    S.append(S[-1]+(rate*N-beta(time*paceOfTime)*S[-1]*I[-1]/N-rate*S[-1])*paceOfTime)
    I.append(I[-1]+(beta(time*paceOfTime)*S[-1]*I[-1]/N-gamma*I[-1]-rate*I[-1])*paceOfTime)
    R.append(R[-1]+(gamma*I[-1]-rate*R[-1])*paceOfTime)
    X.append(time*paceOfTime+START_YEAR)
    time=time+1

print(len(S))
# print(len(numpy.arange(START_YEAR,END_YEAR+paceOfTime)))
plt.scatter(X,S,c='b',marker=',',s= 2,edgecolor='none')
plt.scatter(X,I,c='c',marker=',',s= 2,edgecolor='none')
plt.scatter(X,R,c='g',marker=',',s= 2,edgecolor='none')
plt.legend(['S','I','R'])
plt.xlabel('Year')
plt.ylabel('Population')
plt.show()