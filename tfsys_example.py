import numpy as np
import matplotlib.pyplot as plt
import shapes
import tf,tfftn,tfsys,tfsysftn

#fig,ax = plt.subplots(1,2,figsize=(10,5))
fig,ax = plt.subplots(3,2,figsize=(15,10))

x_lim=[-3,4]
y_lim=[-3,3]

x0=np.linspace(-100,100,2)
y0=np.zeros(x0.shape)

t=np.linspace(0,10,101)

a=np.array([[-1,1,1],[2,0,3],[0,1,1]])
b=np.array([[1,1],[-1,0],[0,2]])
#c=np.array([[1,0,0]])
s=tfsys.Sys(a,b)
s.is_controllable()
p=[[-1,-3],[-1,3],[-2,0]]

################### Controller K -> Pole Placement
k=tfsysftn.sys_poleplace(s,p)
ss=tfftn.tfalg_sys_sI_A(s.A-np.matmul(s.B,k))
char_eq=tfftn.tfalg_sys_char_eq(ss)
#print(char_eq.getRoots())
print(k)

################## Make Closed-loop System(with Controller K)
#k=np.array([[0.595, 0.0281,-0.5389],[1.0779,1.9544,1.1770]])
cont_s=tfftn.tfalg_sys_sI_A(a-np.matmul(b,k))
char_eq1=tfftn.tfalg_sys_char_eq(cont_s)
#print(char_eq1.getRoots())
cont_tfs=tfftn.tfalg_sys_inverse(cont_s)
#print(cont_tfs.tfsys_print())
x_init=np.array([[1],[1],[1]])
####### Tracking #######
x_ss=np.array([[-1],[0],[2]])
u_ss=np.array([[-4/3],[5/3]])
#####################
tfalg_init=tf.Tfalgsys()
tfalg_init.setMat(x_init)
tf_init=tfalg_init.tfalg2tf()


result=tfftn.tfsys_mul(cont_tfs,tf_init,same_roots=True)
result_sys=result.tf_sys
y=np.zeros((3,101))
for i in range(len(result_sys)):
    for j in range(len(result_sys[0])):
        temp=result_sys[i][j]
        y[i]=temp.tf_response(t)+x_ss[i]
        ax[i][0].plot(x0,y0,'k',alpha=0.3)
        ax[i][0].plot(y0,x0,'k',alpha=0.3)
        ax[i][0].plot(t,y[i])
        ax[i][0].set_title(f'Response of {i+1}th state variable')
        ax[i][0].set_xlim([0,10])
        ax[i][0].set_ylim([-3,3])
        ax[i][0].set_xlabel('time(t)')

u=-np.matmul(k,y)
for i in range(2):
    ax[i][1].plot(x0,y0,'k',alpha=0.3)
    ax[i][1].plot(y0,x0,'k',alpha=0.3)
    ax[i][1].plot(t,u[i]+u_ss[i])
    ax[i][1].set_title(f'{i+1}th input')
    ax[i][1].set_xlim([0,10])
    ax[i][1].set_ylim([-3,3])
    ax[i][1].set_xlabel('time(t)')


#fig.savefig("Pole_placement_track.jpg")

plt.show()
