# SIR-MCMC

## 特点
1. 通过计算欧式距离的5次方*（5是我随便选的，只要大于1，本质上应该是一样的）*来衡量预测值与真实值的差距，从而避免使用允许范围较小的math.exp()
2. 通过采用允许误差较大的预测方式，实现了将gamma从100收敛到22左右
3. 对步长sigma和判定收敛的E做了逐步收敛，从而实现大范围的拟合

## 细节说明
1. 我手动调出合理的报告值为0.23
2. 由于在求5次方后数据太大，我直接乘以1e-21
3. 数据导入我直接采用了老师的做法
4. 我对预测的感染人数在整条曲线估算完后做了随机扰动
5. 使用异常处理解决超出范围的问题
6. 在一轮判定收敛失败后会将sigma和E调整为上一次缩小前的值
7. 由于开始的步长较大，容易走到120-150来回游荡，可能需要较多步数才能走到22附近*（运气好往下走的话就不会）*
