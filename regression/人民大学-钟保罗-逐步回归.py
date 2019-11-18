import numpy as np
from numpy.linalg import *
import time


##生成协方差矩阵
def generate_sigma(rho = 0.9,dim_num = 1000):
    sigma = np.eye(dim_num)
    for i in np.arange(dim_num):
        for j in np.arange(i,dim_num):
            sigma[i,j] = rho **(j-i)
    sigma += sigma.T
    sigma -= np.diag(np.ones(dim_num))
    return sigma

#前代法
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
        
        

#gives分解
def gives(mx,lmx):
    mc = mx[0]/lmx
    ms = mx[1]/lmx
    matrix1 = np.array([mc,ms,-ms,mc]).reshape(2,2)
    
#gives矩阵分解（去掉某个变量
def mgives(L,k):
    #L = LA
    #k = 0    
    p = L.shape[1]
    if k>=p :
        return("wrong input of k!")
    Lk = np.delete(L,k,axis = 0)
    mk = k
    while(mk < p-1):
        mx = Lk[mk,mk:(mk+2)].copy()
        lmx = norm(mx)
        Lk[mk,mk:(mk+2)] = [lmx,0]
        if mk < p-2:
            l1 = (Lk[(mk+1):(p-1),mk]*mx[0]/lmx + Lk[(mk+1):(p-1),mk+1]*mx[1]/lmx).copy()
            l2 = (-Lk[(mk+1):(p-1),mk]*mx[1]/lmx + Lk[(mk+1):(p-1),mk+1]*mx[0]/lmx).copy()
            Lk[(mk+1):(p-1),mk] = l1
            Lk[(mk+1):(p-1),mk+1] = l2
        mk = mk+1
    new_Lk = np.delete(Lk,p-1,axis = 1)
    return new_Lk
    
    
#向前更新（增加某个变量
def forupdate(L,xxk,xkxk):
    lk = mforwardsolve(L,xxk)
    lkk = np.sqrt(xkxk -np.sum(lk * lk))
    m = L.shape[0]
    new_L = np.zeros((m+1,m+1))
    new_L[:m,:m] = L
    new_L[m,:m] = lk
    new_L[m,m] = lkk
    return new_L
    

#逐步回归
def steplm(X,Y,method = "AIC"):
    n,p= X.shape
    x = X.copy()
    y = Y.copy() 
    xn = np.array(range(p))
    xtx = np.dot(x.T,x)
    xty = np.dot(x.T,y)
    yty = np.sum(y*y)
    
    L = cholesky(xtx)
    tb = mforwardsolve(L,xty)
    b = mbacksolve(L.T,tb)
    
    RSS = yty - np.sum(tb*tb)
    if method == "BIC":
        eval_f = n*np.log(RSS/n) + p*np.log(n)
    elif method == "CP":
        eval_f = (n-p-1)-n + 2*p   
    else:
        eval_f = n*np.log(RSS/n) + 2*p
    
    #A 索引
    A = np.array(range(p))
    B_1 = np.array(range(p))
    LA = L.copy()
    MEVAL = eval_f
    meval = eval_f
    flag = np.tile(True,p)
    #beta
    hbb = b.copy()
    i = 0
    while(True):
        if len(A)< p:
        ##把删掉的变量放在最后
            B = B_1[~flag]
            AB = np.hstack((A,B))
        else:
            AB = A.copy()
            
        #储存beta用
        bm = np.zeros((p,p))
        #储存得分
        evalm = np.tile(0,p)
        
        ff = 0
        for k in AB:
            if flag[k]:
                #删除第ff个变量
                Lk= mgives(LA,ff)
                #更改变量索引
                tA = np.delete(A,ff)
                xtyk = xty[tA]
                ff = ff+1
            else:
                #加进去变量
                xxk = xtx[A,k]
                xkxk = xtx[k,k]
                #加第ff个变量进去
                Lk = forupdate(LA,xxk,xkxk)
                #更改变量索引
                tA = np.hstack((A,k))
                xtyk = xty[tA]
                
            ###上面的ifelse是为了生成xtyk矩阵
            tbk = mforwardsolve(Lk,xtyk)
            bk = mbacksolve(Lk.T,tbk)
            #求解bk
            
            ## 储存beta在bm里
            bm[k,tA] = bk
            RSSk = yty-np.sum(tbk*tbk) 
            if method == "BIC":
                evalk = n*np.log(RSSk/n)+len(tA)*np.log(n)
            elif method == "CP":
                evalk = (n-p-1)*RSSk/RSS-n+2.0*len(tA)
            else:
                evalk = n*np.log(RSSk/n)+2.0*len(tA)
                
            evalm[k] = evalk
            
            if evalk < meval:
            #存储最小得分的k变量
                mink = k
                mtA = tA.copy()
                meval = evalk
                mLA = Lk.copy()
                hb = bm[k,:]
        #i += 1
        #print(meval)
        #print(i)
        if meval >= MEVAL:
            break
        else:
        #排除最小的
            flag[mink] = ~flag[mink]
            A = mtA.copy()
            MEVAL = meval
            hbb = hb.copy()
            LA = mLA.copy()
    
    re_beta = hbb[flag]
    re = xn[flag]
    return(re_beta,re)
    
    
 

if __name__=="__main__":
    start = time.time()
    n = 100000
    dim_num = 1000
    rho = 0.9
    mu = np.zeros(dim_num)
    sigma = generate_sigma(rho,dim_num = dim_num)
    beta0 = np.array([1,0,-1,0])
    beta0 = np.repeat(beta0,[1,99,1,99])
    beta = np.tile(beta0,5)
    beta_T = ((beta!=0)+ 0).reshape((1,dim_num))[0]
    for met in ["AIC","BIC","CP"]:
        met_acc = []
        for i in range(100):
            time1 = time.time()
            x = np.random.multivariate_normal(mean = mu,cov = sigma,size = n)            
            y = (np.dot(x,beta) + np.random.rand(n))
            re_beta,re = steplm(x,y,method=met)
            print("{0} No.{1} Time used:".format(met,str(i+1)), time.time()-time1)
            print(":\n",re)
            beta_t = np.zeros(dim_num)
            beta_t[re] = 1
            met_acc.append( (beta_t==beta_T).all() )
        print(met_acc)
        print(np.sum(met_acc)/float(len(met_acc)))
    elapsed = (time.time() - start)
    print("Total Time used:",elapsed)
    

'''
最终结论：在大数据集1000维的数据下需要时间较长，在测试100维*10000数据时5.01微秒
'''