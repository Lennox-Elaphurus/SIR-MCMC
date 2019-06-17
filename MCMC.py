import random
from scipy.stats import norm
import math
MAX_PACE = 1000000000000000 #MCMC步数
E=0.0000000001
lk1=0
lk2=0
sigma=10
Continue=0
MIN_CONTINUE=100
A = []
B=[] #参数更新
points=[]
# ilk = norm.logpdf(y,ax+b,1)
LK=[] #记录likelihood
# 读取x跟y
with open(r"test-data-xy.csv","r") as file:
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
        # print(points)
    print(points)
    # a = (points[0][1]-points[1][1])/(points[0][0]-points[1][0])  # 初始值
    # b = points[0][1]-a*points[0][0]  # 初始值
    a=0
    b=0
    # print("a\tb")
    a2=0
    b2=0
for i in range(MAX_PACE):
   # 参数随机游走

    a2 = float(norm.rvs(a,sigma,1))
    # print("a:",a,"a2:", a2)
    # 以上一步的参数值“1”作为真值
    # 通过自己设定sigma来设定步长
    b2 = float(norm.rvs(b,sigma,1))
    # print("b:",b,"b2:", b2)
   # 计算likelihood
    lk1=1
    lk2=1
    for i2 in range(len(points)):
        lk1=lk1+float(norm.logpdf(points[i2][1],a*points[i2][0]+b,sigma))
        lk2=lk2+float(norm.logpdf(points[i2][1],a2*points[i2][0]+b2,sigma))
        # lk1=lk1+norm.pdf(points[i2][1],a*points[i2][0]+b,sigma)
        # lk2 = lk2 + norm.pdf(points[i2][1], a2 * points[i2][0] + b2, sigma)

        # 评估与选择
   #  print("lk2",lk2)
   #  print("lk1",lk1)
   #  R = lk2/lk1
    R=math.exp(lk2-lk1)
    # print("R:",R)
    if random.random() < R:
        if len(A)> 1:
            lastA = A[len(A) - 1]
            lastB = B[len(B) - 1]
            if abs(a-lastA)+abs(b-lastB) < E:
                Continue=Continue+1
                if Continue%10 == 0:# 10 was set by hand
                    sigma=sigma/5   # 5 was set by hand
                    print("Adjust sigma to ",sigma)
                if Continue>MIN_CONTINUE:
                    print("Stop by sufficiently accurate answer.")
                    print("Last accepted R:", R)
                    break
            else:
                Continue=0

        a=a2
        b=b2

        A.append(a)
        B.append(b)
        print("Accepted R:",R)
        print(a,b)

print("Final answer:")
print(a2,b2)