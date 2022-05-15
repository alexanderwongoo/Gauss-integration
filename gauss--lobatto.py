import numpy as np
import math
import matplotlib.pyplot as plt
from sympy import diff
from sympy import symbols
def function(x,n):#定义积分函数
    return x**n

def factorial(n): #计算阶乘
    ns = 1
    for i in range(1, n + 1):
        ns = ns * i
    return ns

def func(n,x):#构造legendre多项式
    l=math.floor(n/2)+1
    sum=0
    for i in range(l):
        sum+=x**(n-2*i)*((-1)**i*factorial(2*n-2*i))/(factorial(i)*factorial(n-i)*factorial(n-2*i))
    temp=(1/2)**n*sum
    return temp

def diff_x(n): #计算n-1legendre项式的导数
    x = symbols("x")
    temp=diff(func(n-1,x),x)
    return temp

def search_xi(n): #找到插值节点
    t = diff_x(n)
    t = t.as_poly()
    xi=np.roots(t.all_coeffs())
    return xi

def gauss_lerangde(a,b,n):
    omiga_i=[2/(n*(n-1))]#边界权重单独计算
    t = diff_x(n)
    t = t.as_poly()
    xi=np.roots(t.all_coeffs())
    xi=np.sort(xi)
    sum=0
    for i in range(len(xi)):
        temp=2/(n*(n-1)*(func(n-1,xi[i])**2))#中间节点权重
        omiga_i.append(temp)
    omiga_i.append(2/(n*(n-1)))
    for i in range(1,n-1):
        sum+=omiga_i[i]*function(xi[i-1],n)
    sum+=omiga_i[0]*function(a,n)+omiga_i[n-1]*function(b,n)
    return sum
print(gauss_lerangde(-1,1,6))
# #记勒让德多项式为pn(x)
plt.rcParams["font.sans-serif"] = ["SimHei"]
plt.rcParams['axes.unicode_minus'] = False
t = np.linspace(-1, 1, 3000)
ax = plt.figure(figsize=(9,5)).add_subplot(1,1,1)
ax.plot(t, func(5,t), label="勒让德多项式曲线")
plt.legend()
plt.show()
