#Voronin Oleksiy
#Zalikova 7203
#variant 3
import random
import math
import matplotlib.pyplot as plt

def signal(n, omega, min_value=0, max_value=1):
    A = [min_value + (max_value - min_value) * random.random() for _ in range(n)]
    phi = [min_value + (max_value - min_value) * random.random() for _ in range(n)]

    def f(t):
        x = 0
        for i in range(n):
            x += A[i]*math.sin(omega/n*t*i + phi[i])
        return x
    return f
def get_m_D(x):
    m = sum(x)/len(x)
    return m, sum([(i - m) ** 2 for i in x]) / (len(x) - 1)
def get_m(x):
    return sum(x)/len(x)
def get_D(x, m=None):
    if m is None:
        m = get_m(x)
    return sum([(i - m) ** 2 for i in x]) / (len(x) - 1)
def get_R(x_gen, y_gen, N, tau=0):
    x = [x_gen(i) for i in range(N)]
    y = [y_gen(i+tau) for i in range(N)]
    m_x, D_x = get_m_D(x)
    m_y, D_y = get_m_D(y)

    R = 0
    for i in range(N//2-1):
        R += (x[i] - m_x)*(y[i+tau] - m_y)/(N//2-1)
    R /= (D_x*D_y)**(1/2)
    return R
def get_F(x):
    N = len(x)
    FR = []
    Fi = []
    for p in range(N):
        FR.append(0)
        Fi.append(0)
        for k in range(N):
            FR[p] += x[k]*math.cos(-2*math.pi*p*k/N)
            Fi[p] += x[k]*math.sin(-2*math.pi*p*k/N)
    return FR, Fi
def get_F_optimized(x):
    N = len(x)
    w = []
    for i in range(N):
        if i < N//4:
            w.append((math.cos(-2*math.pi*i/N),
                      math.sin(-2*math.pi*i/N)))
        elif i < N//2:
            w.append((w[i-N//4][1],
                      -w[i-N//4][0]))
        else:
            w.append((-w[i-N//2][0], -w[i-N//2][1]))
    FR = []
    Fi = []
    for i in range(N):
        FR.append(sum([w[(i*j)%N][0]*x[j] for j in range(N)]))
        Fi.append(sum([w[(i*j)%N][1]*x[j] for j in range(N)]))
    return FR, Fi
def maketable(x,y,i):
    x.append(round(y[i*round(N/10)],2))
n = 8
omega = 1100
N = 256

range_min = 0
range_max = 1

x_gen = signal(n, omega, range_min, range_max)
x = [x_gen(i) for i in range(N)]
(FR, Fi) = get_F(x)
(FR1, Fi1) = get_F_optimized(x)

F = [FR[i] + Fi[i] for i in range(N)]

Ntable = []
FRtable = []
Fitable = []
FR1table = []
Fi1table = []
Ftable = []
xtable = []
rows = ["FR", "Fi", "FR1", "Fi1", "x", "F"]
barcord = []
barheight = []
for i in range(round(N/25)):
    Ntable.append(i*round(N/10))
    maketable(FRtable, FR, i)
    maketable(Fitable, Fi1, i)
    maketable(FR1table, FR1, i)
    maketable(Fi1table, Fi1, i)
    maketable(xtable, x, i)
    maketable(Ftable, F, i)
    barcord.append(Ntable[i])
    barheight.append(FRtable[i]+FR1table[i]+Fitable[i]+Fi1table[i]+Ftable[i]+x[i])
for row in range(len(Ntable)):
    plt.bar(barcord, barheight, width=20)
tabletext = [FRtable, Fitable, FR1table, Fi1table, xtable, Ftable]
the_table = plt.table(cellText=tabletext,rowLabels=rows,colLabels=Ntable)
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.show()
