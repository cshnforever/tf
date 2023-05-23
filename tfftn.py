import numpy as np
import itertools
from tf import *

def sort(a):
    a1=list(a)
    count=0
    for i in range(len(a)-1):
        for j in range(len(a)-i-1):
            if a1[j]>a1[j+1]:
                temp=a1[j]
                a1[j]=a1[j+1]
                a1[j+1]=temp
                count+=1
    return count

########## tfalg ############

def tfalg_mul(tf1,tf2,coef_mode=False):
    s=Tfalg()
    a=tf1.coef
    b=tf2.coef
    if a==[0] or b==[0]:
        s.setRoots([],0)
        return s
    if coef_mode==False:
        s.setRoots(tf1.getRoots()+tf2.getRoots(),tf1.n_coef*tf2.n_coef)
    else:
        coef=np.zeros((len(a)+len(b)-1))
        for i in range(len(a)):
            for j in range(len(b)):
                coef[i+j] += (a[i]*b[j])
        s.setCoef(list(coef))
    return s

def tfalg_mul_ls(ls,coef_mode=False):
    s=Tfalg()
    if coef_mode==False:
        r=[]
        n=1
        for tf in ls:
            if tf.coef==[0]:
                s.setCoef([0])
                return s
            n = n * tf.n_coef
            r += tf.getRoots()
        s.setRoots(r,n)
    else:
        s.setCoef([1])
        for tf in ls:
            if tf.coef==[0]:
                s.setCoef([0])
                return s
            s=tfalg_mul(s,tf)

    return s

def tfalg_sum(tf1,tf2,coef_mode=False):
    tfalg=Tfalg()
    a=tf1.getCoef()
    b=tf2.getCoef()
    n=max(len(a),len(b))
    a1=[0]*(n-len(a))+a
    b1=[0]*(n-len(b))+b

    coef=list(np.array(a1)+np.array(b1))
    for i in range(len(coef)):
        if coef[i]!=0:
            coef1=coef[i:]
            break
        if i==len(coef)-1:
            coef1=[0]
    tfalg.setCoef(coef1)

    if coef_mode==False:
        tfalg.setRoots(tfalg.find_roots(),coef1[0])

    return tfalg

def tfalg_minus(tf1,tf2,coef_mode=False):
    tfalg=Tfalg()
    a=tf1.getCoef()
    b=tf2.getCoef()
    n=max(len(a),len(b))
    a1=[0]*(n-len(a))+a
    b1=[0]*(n-len(b))+b
    coef=list(np.array(a1)-np.array(b1))
    for i in range(len(coef)):
        if coef[i]!=0:
            coef1=coef[i:]
            break
        if i==len(coef)-1:
            coef1=[0]
    tfalg.setCoef(coef1)
    if coef_mode==False:
        tfalg.setRoots(tfalg.find_roots(),coef1[0])
    
    return tfalg

########## tf #################3


def tf_mul(tf1,tf2,coef_mode=False):
    s=Tf()
    if coef_mode==False:
        s1=tfalg_mul(tf1.num,tf2.num)
        s2=tfalg_mul(tf1.den,tf2.den)
        s.setZeros(s1.find_roots(),s1.coef[0])
        s.setPoles(s2.find_roots(),s2.coef[0])
    else:
        s1=tfalg_mul(tf1.num,tf2.num,coef_mode=True)
        s2=tfalg_mul(tf1.den,tf2.den,coef_mode=True)
        s.num=s1
        s.den=s2
    return s

def tf_sum(tf1,tf2,coef_mode=False):
    s=Tf()
    if coef_mode==False:
        s11=tfalg_mul(tf1.num,tf2.den)
        s12=tfalg_mul(tf1.den,tf2.num)
        s1=tfalg_sum(s11,s12)
        s2=tfalg_mul(tf1.den,tf2.den)
        s.setZeros(s1.find_roots(),s1.coef[0])
        s.setPoles(s2.find_roots(),s2.coef[0])
    else:
        s11=tfalg_mul(tf1.num,tf2.den,coef_mode=True)
        s12=tfalg_mul(tf1.den,tf2.num,coef_mode=True)
        s1=tfalg_sum(s11,s12,coef_mode=True)
        s2=tfalg_mul(tf1.den,tf2.den,coef_mode=True)
        s.num=s1
        s.den=s2
    return s


############### tfalg_sys ###################

def tfalg_sys_sI_A(A):
    a=Tfalgsys()
    if A.shape[0]!=A.shape[1]:
            print("Tfalgsys(sI-A): Not square Matrix")
            return
    n=A.shape[0]
    a.setSize(n,n)

    for i in range(n):
        for j in range(n):
            if i==j:
                a.tfalg_sys[i][j].setCoef([1,-A[i][j]])
            else:
                a.tfalg_sys[i][j].setCoef([-A[i][j]])
    return a

def tfalg_sys_char_eq(a,coef_mode=False):
    n=len(a.tfalg_sys)
    iter_ls=list(itertools.permutations(list(range(n)), n))
    char_eq=Tfalg()
    char_eq.setCoef([0])
    for it in iter_ls:
        b=Tfalg()
        for i in range(len(it)):
            b=tfalg_mul(b,a.tfalg_sys[i][it[i]],True)
        if sort(it)%2==0:
            char_eq=tfalg_sum(char_eq,b,True)
        else:
            char_eq=tfalg_minus(char_eq,b,True)

    if coef_mode==False:
        char_eq.setRoots(char_eq.find_roots(),char_eq.coef[0])

    return char_eq

def tfalg_sys_inverse(a):
    n=len(a.tfalg_sys)
    adj=Tfalgsys()
    adj.setSize(n,n)

    for i in range(n):
        for j in range(n):
            ls1=list(range(i))+list(range(i+1,n))
            ls2=list(range(j))+list(range(j+1,n))
            small=a.crop(ls1,ls2)
            if (i+j)%2==0:
                adj.tfalg_sys[j][i]=tfalg_sys_char_eq(small)
            else:
                aa=Tfalg()
                aa.setRoots([],-1)
                adj.tfalg_sys[j][i]=tfalg_mul(tfalg_sys_char_eq(small),aa)

    det=tfalg_sys_char_eq(a,False)
    roots=det.getRoots()

    inv=Tfsys()
    inv.setSize(n,n)
    for i in range(n):
        for j in range(n):
            inv.tf_sys[i][j].setZeros(adj.tfalg_sys[i][j].find_roots(),adj.tfalg_sys[i][j].coef[0])
            inv.tf_sys[i][j].setPoles(roots,det.coef[0])
    return inv

############### tf_sys

def tfsys_mul(tf1,tf2,same_roots=False):

    a=tf1.tf_sys
    b=tf2.tf_sys
    if len(a[0])!=len(b):
        print("tfsys_mul: Not match shape of Matrix")
        return
    m=len(tf1.tf_sys)
    n=len(tf1.tf_sys[0])
    r=len(tf2.tf_sys[0])
    aa=Tfsys()
    aa.setSize(m,r)

    if same_roots==False:
        for i in range(m):
            for k in range(r):
                x=Tf()
                x.setZeros([],0)
                for j in range(n):
                    temp=tf_mul(tf1.tf_sys[i][j],tf2.tf_sys[j][k],True)
                    x=tf_sum(x,temp,True)
                aa.tf_sys[i][k]=x

    else:
        den=tfalg_mul(tf1.tf_sys[0][0].den,tf2.tf_sys[0][0].den,True)
        for i in range(m):
            for k in range(r):
                x=Tfalg()
                x.setRoots([],0)
                for j in range(n):
                    temp=tfalg_mul(tf1.tf_sys[i][j].num,tf2.tf_sys[j][k].num,True)
                    x=tfalg_sum(x,temp,True)
                aa.tf_sys[i][k].num=x
                aa.tf_sys[i][k].den=den

    print(aa.tfsys_print())
    aa.setRoots(same_roots)
    
    return aa

def tfsys_sum(tf1,tf2):

    a=tf1.tf_sys
    b=tf2.tf_sys

    if len(a)!=len(b) or len(a[0])!=len(b[0]):
        print("tfsys_sum: Not match shape of Matrix")
        return
    m=len(tf1.tf_sys)
    n=len(tf1.tf_sys[0])

    aa=Tfsys()
    aa.setSize(m,n)

    for i in range(m):
        for j in range(n):
            aa.tf_sys[i][j]=tf_sum(a[i][j],b[i][j])
    
    return aa