import numpy as np

def comp_sum(a,b):
    return [a[0]+b[0],a[1]+b[1]]

def comp_mul(a,b):
    return [a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]]

def comp_sum_ls(ls):
    s=[0,0]
    for elem in ls:
        elem=list(elem)
        s = comp_sum(s,elem)
    return s

def comp_mul_ls(ls):
    s=[1,0]
    for elem in ls:
        elem=list(elem)
        s = comp_mul(s,elem)
    return s