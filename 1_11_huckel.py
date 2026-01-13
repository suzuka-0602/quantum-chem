#行列式の計算を、LU分解を用いた方法として実装
#永年方程式を二分法を用いて解く

#ここでは自分で実装したけど、行列式の計算と永年方程式の解の計算はscipyでやるのが楽

import numpy as np

#-----設定-----
n=100 #結合に参加する原子数
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

# def ethylene(e):
#     return [
#         [a-e, b],
#         [b, a-e]
#     ]

# def butadiene(e):
#     return [
#         [a-e, b,   0,   0],
#         [b,   a-e, b,   0],
#         [0,   b,   a-e, b],
#         [0,   0,   b,   a-e]
#     ]

# def toluene(e):
#     return [
#         [a-e, 0, 0, 0, b, 0, 0],
#         [0, a-e, 0, 0, 1, 1, 0],
#         [0, 0, a-e, 0, 0, b, b],
#         [0, 0, 0, a-e, b, 0, b],
#         [b, b, 0, b, a-e, 0, 0],
#         [0, b, b, 0, 0, a-e, 0],
#         [0, 0, b, b, 0, 0, a-e]
#     ]


#-----LU分解-----
def LUD(matrix): #分解したい行列を引数、分解したL,Uを返す関数
    n=len(matrix) #n×n行列
    L=np.zeros((n, n)) #分解した下三角行列用
    U=np.eye(n) #分解した上三角行列用

    A=matrix
    for i in range(n):
        #l00=a00
        l00=A[0][0]

        #l1=a1
        l1=np.zeros((n-i-1, 1))
        for j in range(n-i-1):
            l1[j]=A[j+1][0]

        #u1^T=α1^T/l00
        u1=np.zeros((1,n-i-1))
        for j in range(n-i-1):
            u1[0][j]=A[0][j+1]/l00

        #A1=l1@u1^T+L1@U1　 → 　AをA1-l1@u1^Tに置き換え
        A1=np.zeros((n-i-1, n-i-1)) #A1はAの一番目の行と列を取り除いたもの
        l1u1=l1@u1
        for j in range(n-i-1):
            for k in range(n-i-1):
                A1[j][k]=A[j+1][k+1]-l1u1[j][k]
        A=A1


        #l1をLのi列目に代入
        L[i][i]=l00
        for j in range(len(l1)):
            L[j+i+1][i]=l1[j][0]
        #u1をUのi行目に代入 
        for j in range(len(l1)):
            U[i][j+i+1]=u1[0][j]

    return L, U

#-----行列式計算-----
def det(matrix):
    L, U=LUD(matrix)
    total_det=1.0
    for i in range(len(matrix)):
        total_det*=L[i][i]
    return total_det


#-----二分法-----
region=[-10.0, 10.0] #探索範囲
step=0.1
eps=1e-4 #許容誤差

def judge(u:float, v:float):
    return det(matrix(u))*det(matrix(v))

def bisection(left, right): #(left, right)内の零点を探索
    dx=right-left
    while dx > eps: #dxが誤差範囲に収まるまで
        x=(right+left)/2
        if judge(x, left)<0: #(x, left)内に零点があるならその範囲を探索
            right=x
        else: #(x, left)外にあるなら範囲外(right, x)を探索
            left=x
        dx=abs(right-left)
    else:
        print('{:10.6f} hartree'.format(x))

def main():
    x0=region[0]
    while x0 <= region[1]-step: #探索範囲内をstepずつシフト
        if judge(x0, x0+step) < 0: #step内に零点があれば
            bisection(x0, x0+step) #二分法で零点を探す
        x0+=step

main()