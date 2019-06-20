# 采用每步根据真实值推出下一个预测值的方式
import random
from scipy.stats import norm
import math

fig = r'simulation'
START_YEAR = 2008
END_YEAR=2020
dt = 1/52
mu = 1/50
gamma = 22
reportRate=0.5
pop = 5000000   # initial population

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
                points.append(float(temp),float(pointlist[idx]))
                flag=True


import_data("sir_case.csv")
for i in range((END_YEAR-START_YEAR)/dt):
    print(points[i].second())

# def get_likelihood(time,sigma):
    # return float(norm.logpdf(points[i2][1], a2 * points[i2][0] + b2, sigma)
