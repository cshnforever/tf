import numpy as np
import itertools

class Sys:
    def __init__(self,A,B,C=0,D=0):
        if A.shape[0]!=A.shape[1]:
            print("Not match A's Shape")
            return
        if A.shape[1]!=B.shape[0]:
            print("Not match number of variable and input")
            return
        self.A=A
        self.B=B
        self.C=C
        self.D=D

        self.len_x=A.shape[1]
        self.len_u=B.shape[1]

        if type(C)==int:
            self.C=np.eye(self.len_x)
        self.len_y=self.C.shape[0]

        if type(D)==int:
            self.D=np.zeros((self.len_y,self.len_u))

    def is_controllable(self):
        cont_mat=self.B
        temp=self.B
        for i in range(1,self.len_x):
            temp=np.matmul(self.A,temp)
            cont_mat=np.hstack((cont_mat,temp))
        
        rank=np.linalg.matrix_rank(cont_mat)
        print(cont_mat)
        if rank<self.B.shape[0]:
            print("Uncontrollable (Not full rank)")
            return
        else:
            print("Full rank")
        if cont_mat.shape[0]==cont_mat.shape[1]:
            a=np.linalg.det(cont_mat)
            if abs(a)<1e-3:
                print("May be Uncontrollable (Small Determinant)")
            else:
                print("Controllable")
        


    def sys_print(self):
        print("System Resume")
        print("Numbers of Variables: "+str(self.len_x))
        print("Numbers of Inputs: "+str(self.len_u))
        print("Numbers of Outputs: "+str(self.len_y))
        print("A: ")
        print(self.A)
        print("B: ")
        print(self.B)
        print("C: ")
        print(self.C)
        print("D: ")
        print(self.D)