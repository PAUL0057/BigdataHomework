#coding:utf-8

#加载包
import time
import numpy as np
import numpy.linalg as alg
from sklearn.model_selection import KFold


##生成协方差矩阵
def generate_sigma(rho = 0.9,dim_num = 1000):
    sigma = np.eye(dim_num)
    for i in np.arange(dim_num):
        for j in np.arange(i,dim_num):
            sigma[i,j] = rho **(j-i)
    sigma += sigma.T
    sigma -= np.diag(np.ones(dim_num))
    return sigma
    
#前代
def mforwardsolve(L,b):
    m,n = L.shape
    a = b.copy()
    if m != n:
        print("Wrong dimensions of matrix L!")
        return
    else:
        x = np.zeros(m)
        for i in range(m):
            x[i] = a[i]/L[i,i]
            if i< (m-1):
                a[(i+1):] = a[(i+1):] - x[i] * L[(i+1):,i]
        return x
        
#回代法
def mbacksolve(L,b):
    m,n = L.shape
    a = b.copy()
    if m != n or m!= len(b):
        print("Wrong dimensions of matrix L!")
        return
    else:
        x = np.zeros(m)
        for i in range(m-1,-1,-1):
            x[i] = a[i]/L[i,i]
            if i>0:
                a[:i] = a[:i] - x[i]*L[:i,i] 
        return x
        
        

class Ridge_regressor():
    '''
    岭回归
    '''
    def __init__(self,x_train,x_test,y_train,y_test,lambd):
    ##初始化设置
        self.X_train = x_train
        self.y_train = y_train
        self.X_test = x_test
        self.y_test = y_test
        self.lambd = lambd
    
    def cal_beta(self):
    ##计算beta系数
        n,p = self.X_train.shape
        V = np.dot(self.X_train.T,self.X_train) + self.lambd * np.eye(p)
        U = np.dot(self.X_train.T,self.y_train)
        L = alg.cholesky(V)
        cache = mforwardsolve(L,U)
        self.beta = mbacksolve(L.T,cache)
        
        self.train_error = (np.sum(self.y_train**2) - np.sum(cache**2)) / (n-p)
        
    def train_mse(self):
    ##计算训练集MSE
        return self.train_error
        
    def test_mse(self):
    ##计算测试集MSE
        n,p = self.X_train.shape
        pred = np.dot(self.X_test,self.beta)
        self.test_error = np.sum((self.y_test - pred)**2) / (n-p)
        return self.test_error

def kcv_RidgeReg(lambd,x,y,k=10):
    ##10折交叉验证
    n,p = x.shape
    kf = KFold(k)
    train_error_cache = []
    test_error_cache = []
    for train_index, test_index in kf.split(x):
        reg = Ridge_regressor(x[train_index,:],x[test_index,:],y[train_index],y[test_index],lambd)
        reg.cal_beta()
        train_error_cache.append(reg.train_mse())
        test_error_cache.append(reg.test_mse())
    train_error = sum(train_error_cache) / len(train_error_cache)
    test_error = sum(test_error_cache) / len(test_error_cache)
    return train_error,test_error


def lattice(x,y):
    
    lat = np.arange(0,1,0.01)
    f = []
    for i in lat:
        f0 = kcv_RidgeReg(i,x,y)[1]
        f.append(f0)
    fin=min(f)
    return lat[f.index(fin)],fin


def tripple(x,y,x0=0,x1=1,L=100):
    
    for i in range(L):
        x2=x0+(x1-x0)/3
        x3=x1-(x1-x0)/3
        f2=kcv_RidgeReg(x2,x,y)[1]
        f3=kcv_RidgeReg(x3,x,y)[1]
        if f2<f3:
            x1=x3
            if (x1-x0)<0.01:
                return x2,f2
        elif f2>f3:
            x0=x2
            if (x1-x0)<0.01:
                return x3,f3
    x_fin=(x2+x3)/2
    f_fin=(f2+f3)/2
    return x_fin,f_fin

#%%
if __name__ == '__main__':
    '''
    主程序
    '''
    start = time.time()
    n = 100000 #样本数
    dim_num = 1000#变量数
    rho = 0.9#rho值
    mu =np.zeros(dim_num)
    sigma = generate_sigma(rho,dim_num)
    beta = np.tile(np.hstack([1,-1]),500)
    x = np.random.multivariate_normal(mu,sigma,size = n)
    y = np.dot(x,beta) + np.random.randn(n)
    
    lambd1,mse1 = lattice(x,y)
    lambd2,mse2 = tripple(x,y)
    
    print("隔点法，最优lambda:{:.4f},对应的MSE:{:.4f}".format(lambd1,mse1))
    print("割线法，最优lambda:{:.4f},对应的MSE:{:.4f}".format(lambd2,mse2))
    
    '''
    最后的运行结果是：
    隔点法，最优lambda:0.5500,对应的MSE:0.1019
    割线法，最优lambda:0.5400,对应的MSE:0.1009
    '''