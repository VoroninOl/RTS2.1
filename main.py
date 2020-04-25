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

fig = plt.figure()
fig2 = plt.figure()

ax_1 = fig.add_subplot(3, 1, 1)
ax_2 = fig.add_subplot(3, 1, 2)
ax_3 = fig.add_subplot(3, 1, 3)

ax_4 = fig2.add_subplot(3, 1, 1)
ax_5 = fig2.add_subplot(3, 1, 2)
ax_6 = fig2.add_subplot(3, 1, 3)


ax_1.plot(range(N), FR)
ax_2.plot(range(N), Fi)
ax_3.plot(range(N), FR1)
ax_4.plot(range(N), Fi1)
ax_5.plot(range(N), x)
ax_6.plot(range(N), F)

ax_1.set(title='FR')
ax_2.set(title='Fi')
ax_3.set(title='FR_opt')
ax_4.set(title='Fi_opt')
ax_5.set(title='x')
ax_6.set(title='F')

plt.show()
