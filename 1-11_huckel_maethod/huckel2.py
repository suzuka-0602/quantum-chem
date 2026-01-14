import numpy as np
from scipy.optimize import root_scalar

#-----設定-----
n=7 #結合に参加する原子数
combination=[ #結合の仕方をインデックスで管理
    [0,4],
    [4,1],
    [1,5],
    [5,2],
    [2,6],
    [6,3],
    [3,4]
]

a=0 #coulomb積分
b=1 #共鳴積分


#-----行列生成-----
def matrix(e):
    matrix=np.zeros((n,n))
    for i in range(n):
        matrix[i][i]=a-e
    for comb in combination:
        index1, index2=comb
        matrix[index1][index2]=b
        matrix[index2][index1]=b
    return matrix

#永年方程式を解く
def det(e):
    return np.linalg.det(matrix(e))
result=root_scalar(det, bracket=(1, 10))
print(result.root)