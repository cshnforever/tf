import numpy as np
import tf
import tfftn
import tfsys

def sys2tfsys(s):

    a=s.A
    b=s.B
    c=s.C
    d=s.D

    a1=tfftn.tfalg_sys_sI_A(a)
    inva=tfftn.tfalg_sys_inverse(a1)
    b1=tf.Tfalgsys()
    b1.setMat(b)
    b2=b1.tfalg2tf()
    c1=tf.Tfalgsys()
    c1.setMat(c)
    c2=c1.tfalg2tf()
    d1=tf.Tfalgsys()
    d1.setMat(d)
    d2=d1.tfalg2tf()

    t1=tfftn.tfsys_mul(c2,inva,True)
    t2=tfftn.tfsys_mul(t1,b2,True)
    t=tfftn.tfsys_sum(t2,d2)

    return t


def sys_poleplace(s,p_ls,step=10000,num=3):

    cost=10000*np.ones((num,1))
    kk=np.zeros((num,s.len_u,s.len_x))

    for iter in range(num):

        k=np.random.randn(s.len_u,s.len_x)

        dk=np.ones(k.shape)

        for ii in range(1,step):
            g=np.zeros(k.shape)

            s1=tfftn.tfalg_sys_sI_A(s.A-np.matmul(s.B,k))
            char_eq1=tfftn.tfalg_sys_char_eq(s1,True)
            cost1=0
            for p in p_ls:
                x=char_eq1.calc(p)
                cost1 += np.sqrt(x[0]*x[0]+x[1]*x[1])
            cost1 /= len(p_ls)
            
            for i in range(s.len_u):
                for j in range(s.len_x):
                    k[i][j] = k[i][j] + dk[i][j]
                    s2=tfftn.tfalg_sys_sI_A(s.A-np.matmul(s.B,k))
                    char_eq2=tfftn.tfalg_sys_char_eq(s2,True)
                    cost2=0
                    for p in p_ls:
                        x=char_eq2.calc(p)
                        cost2 += np.sqrt(x[0]*x[0]+x[1]*x[1])
                    cost2 /= len(p_ls)

                    k[i][j] = k[i][j] - dk[i][j]
                    g[i][j]=(cost2-cost1)/(dk[i][j]+1e-6)

            abs_g=np.linalg.norm(g)
            if cost1<1e-2 or abs_g<1e-3:
                break
            g=g/abs_g
            dk=-(g)/(ii+1)*2
            k=k+dk
        
        cost[iter]=cost1
        kk[iter]=k
    
    min_cost=10000
    min_ind=0

    for iter in range(num):
        if cost[iter]<min_cost:
            min_ind=iter
            min_cost=cost[iter]

    return kk[min_ind]