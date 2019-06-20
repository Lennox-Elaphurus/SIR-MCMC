import matplotlib.pyplot as plt
import math

fig = r'simulation'

# parameters
t0 = 2008
dt = 1/365
mu = 1/50
gamma = 22
reportRate=0.5

# initial population
pop = 5000000
s0 = 26/400
i0 = 0.002
r0 = 1
S0 = pop*s0/(s0+i0+r0)
I0 = pop*i0/(s0+i0+r0)
R0 = pop*r0/(s0+i0+r0)

# state variable/list
S = [S0]
I = [I0]
R = [R0]
T = [t0]
Beta = [180*(5+2*math.sin(math.pi*t0+5))]

# update state variables each day
for i in range((2020-2008)*365):
    S0=S[-1]
    I0=I[-1]
    R0=R[-1]
    T0=T[-1]
    pop = S0+I0+R0
    beta = 180*(5+2*math.sin(math.pi*(T0+dt)+5))
    S1 = S0 + (mu*pop - beta*S0*I0/pop - mu*S0)*dt
    I1 = I0 + (beta*S0*I0/pop - gamma*I0 - mu*I0)*dt
    R1 = R0 + (gamma*I0 - mu*R0)*dt
    T1 = T0 + dt
    Beta.append(beta)
    S.append(S1)
    I.append(I1)
    R.append(R1)
    T.append(T1)

plt.plot(T,S,'b--')
plt.plot(T,I,'g--')
plt.plot(T,R,'r--')
plt.text(2008,4000000,"S",fontsize=24,color='b')
plt.text(2008,3000000,"I",fontsize=24,color='g')
plt.text(2008,2000000,"R",fontsize=24,color='r')
plt.savefig(fig)

plt.show()