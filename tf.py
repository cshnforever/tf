import numpy as np
import itertools
import cpx

def tf_mul(tf1,tf2):
    s=Tfalg()
    s.setRoots(tf1.getRoots()+tf2.getRoots(),tf1.n_coef*tf2.n_coef)
    return s

def tf_mul_ls(ls):
    s=Tfalg()
    r=[]
    n=1
    for tf in ls:
        n = n * tf.n_coef
        r += tf.getRoots()
    s.setRoots(r,n)
    return s

def fact(n):
    f=1
    for i in range(n):
        f *= (i+1)
    return f

class Tfalg:
    def __init__(self):
        self.root_dict=dict()
        self.n_coef=1
        self.order=0
        self.coef=[1]

    def setCoef(self,coef):
        self.coef=coef
        self.n_coef=coef[0]
        self.order=len(coef)-1

    def setRoots(self,roots,n):
        ls=[]
        set_roots=set(roots)
        for root in set_roots:
            m=0
            for i in range(len(roots)):
                if root[0]==roots[i][0] and root[1]==roots[i][1]:
                    m+=1
            ls.append([root,m])
        self.root_dict=dict(ls)

        self.n_coef=n
        self.coef=self.calc_coef()
        self.order=len(self.coef)-1

    def getRoots(self):
        roots=[]
        for key,val in self.root_dict.items():
            for i in range(val):
                roots.append(key)
        roots=sorted(roots)
        return roots
    
    def getCoef(self):
        return self.coef

    def calc_coef(self):
        roots=self.getRoots()
        n=len(roots)
        if n==0:
            coef=[self.n_coef]
            return coef
        coef=[]*n
        for i in range(n):
            s=[0,0]
            ls=list(itertools.combinations(roots,i+1))
            for elem_ls in ls:
                s=cpx.comp_sum(s,cpx.comp_mul_ls(elem_ls))
            if i%2==0:
                s[0]=-s[0]
            coef.append(s[0]*self.n_coef)
        coef=[self.n_coef]+coef
        return coef
    
    def calc(self,x):
        n=len(self.coef)-1
        a=0
        b=0
        for i in reversed(range(n+1)):
            l=[x]*i
            
            m=cpx.comp_mul_ls(l)

            a += self.coef[n-i]*m[0]
            b += self.coef[n-i]*m[1]
        return [a,b]

    def find_roots(self,step=100000):
        n=len(self.coef)-1
        ls=[]
        if n==1:
            x=list(np.around(np.array([-self.coef[-1]/self.coef[0],0]),4))
            ls.append(x)
            #print(ls)
            return ls
        elif n==2:
            a=self.coef[0]
            b=self.coef[1]
            c=self.coef[2]
            if b*b-4*a*c>=0:
                x1=list(np.around(np.array([(-b+np.sqrt(b*b-4*a*c))/(2*a),0]),4))
                x2=list(np.around(np.array([(-b-np.sqrt(b*b-4*a*c))/(2*a),0]),4))
                ls.append(x1)
                ls.append(x2)
            else:
                x1=list(np.around(np.array([(-b)/(2*a),np.sqrt(4*a*c-b*b)/(2*a)]),4))
                x2=list(np.around(np.array([(-b)/(2*a),-np.sqrt(4*a*c-b*b)/(2*a)]),4))
                if abs(x1[1])<1e-2:
                    x1[1]=0
                    x2[1]=0
                ls.append(x1)
                ls.append(x2)
            #print(ls)
            return ls
        else:
            aa=1
            
        x=list(np.random.randn(2))
        dx=0.1
        dy=0.1

        for i in range(1,step):
        #print(x)
            p=self.calc(x)
            px=self.calc([x[0]+dx,x[1]])
            py=self.calc([x[0],x[1]+dy])
            
            gx=(2*p[0]*(px[0]-p[0])+2*p[1]*(px[1]-p[1]))/(dx+1e-9)
            gy=(2*p[0]*(py[0]-p[0])+2*p[1]*(py[1]-p[1]))/(dy+1e-9)

            g=np.sqrt(gx*gx+gy*gy)

            dx = -(gx/g)*(1/(i+10))
            dy = -(gy/g)*(1/(i+10))
            d=np.sqrt(dx*dx+dy*dy)
            
            x[0]=x[0]+dx
            x[1]=x[1]+dy
            '''
            if i%1000==0:
                print(x)
                print(p[0]*p[0]+p[1]*p[1])
                #plt.plot(x[0],x[1],'ro')
            '''
            if np.sqrt(p[0]*p[0]+p[1]*p[1])<1e-10:
                break

        x=list(np.around(np.array(x),4))
        if abs(x[1])<1e-2:
            x[1]=0
        
        new_tf=Tfalg()
        if x[1]==0:
            ls.append(x)
            A=np.zeros((n+1,n))
            b=np.array(self.coef)
            for i in range(n):
                A[i:i+2,i]=np.array([1,-x[0]])
            new_coef=np.matmul(np.linalg.pinv(A),b)
            new_tf.setCoef(new_coef)
            
        else:
            ls.append(x)
            ls.append([x[0],-x[1]])
            A=np.zeros((n+1,n-1))
            b=np.array(self.coef)
            for i in range(n-1):
                A[i:i+3,i]=np.array([1,-2*x[0],x[0]*x[0]+x[1]*x[1]])
            new_coef=np.matmul(np.linalg.pinv(A),b)
            new_tf.setCoef(new_coef)

        new_ls=new_tf.find_roots()
        ls = ls + new_ls

        return ls


    def tfalg_print(self):
        n=len(self.coef)-1
        t=""
        for i in range(n+1):
            if i==0:
                t = t + f"{self.coef[i]}s^{n-i}"
            elif i==n:
                if self.coef[i]>=0:
                    t = t + f"+{self.coef[i]}"
                else:
                    t = t + f"{self.coef[i]}"
            else:
                if self.coef[i]>=0:
                    t = t + f"+{self.coef[i]}s^{n-i}"
                else:
                    t = t + f"{self.coef[i]}s^{n-i}"
        return t

class Tf(Tfalg):
    def __init__(self):
        self.num=Tfalg()
        self.den=Tfalg()

    def setZeros(self,zeros,n=1):
        self.num.setRoots(zeros,n)
        self.elim()
    
    def setPoles(self,poles,n=1):
        self.den.setRoots(poles,n)
        self.elim()

    def getZeros(self):
        return self.num.getRoots()
    
    def getPoles(self):
        return self.den.getRoots()

    def elim(self):
        elim_num=[]
        elim_den=[]
        for key1 in self.num.root_dict.keys():
            for key2 in self.den.root_dict.keys():
                if key1==key2:
                    if self.num.root_dict[key1]>self.den.root_dict[key1]:
                        self.num.root_dict[key1]=self.num.root_dict[key1]-self.den.root_dict[key1]
                        elim_den.append(key1)
                    elif self.num.root_dict[key1]<self.den.root_dict[key1]:
                        self.den.root_dict[key1]=self.den.root_dict[key1]-self.num.root_dict[key1]
                        elim_num.append(key1)
                    else:
                        elim_den.append(key1)
                        elim_num.append(key1)
        
        for key in elim_num:
            del self.num.root_dict[key]
        for key in elim_den:
            del self.den.root_dict[key]

        self.num.setRoots(self.num.getRoots(),self.num.n_coef)
        self.den.setRoots(self.den.getRoots(),self.den.n_coef)

    def tf_divide(self):
        tf_ls=[]
        var_num_ls=[]
        
        for key,val in self.den.root_dict.items():
            if key[1]==0:
                for i in range(val):
                    tf=Tf()
                    p=[key]*(i+1)
                    tf.setPoles(p)
                    var_num_ls.append(1)
                    tf_ls.append(tf)
            elif key[1]>0:
                for i in range(val):
                    tf=Tf()
                    p=[key,(key[0],-key[1])]*(i+1)
                    tf.setPoles(p)
                    var_num_ls.append(2)
                    tf_ls.append(tf)
            else:
                continue

        var_num=np.sum(var_num_ls)
        A=np.zeros((self.den.order,var_num))
        
        j=0
        for i in range(len(tf_ls)):
            tf1=Tf()
            tf1.setZeros(self.den.getRoots())
            tf1.setPoles(tf_ls[i].den.getRoots())
            for k in reversed(range(var_num_ls[i])):
                a=[0]*(self.den.order-len(tf1.num.calc_coef()))+tf1.num.calc_coef()+[0]*k
                A[:,j]=a[k:]
                j+=1

        b=np.array([0]*(self.den.order-len(self.num.calc_coef()))+self.num.calc_coef())
        b = b / self.num.n_coef

        c=np.matmul(np.linalg.inv(A),b).reshape(-1)

        j=0
        for i in range(len(tf_ls)):
            tf_num_coef=np.around(np.array(c[j:j+var_num_ls[i]])*self.num.n_coef,4)
            j += var_num_ls[i]
            tf_ls[i].num.setCoef(tf_num_coef)
            tf_ls[i].den.n_coef=self.den.n_coef

        return tf_ls

    def tf_2nd_divide(self):
        den_coef=self.den.calc_coef()
        num_coef=[0,0]
        if len(self.num.coef)==2:
            num_coef=self.num.coef
        else:
            num_coef[1]=self.num.coef
        a=den_coef[0]
        b=den_coef[1]
        c=den_coef[2]
        d=num_coef[0]
        e=num_coef[1]
        
        if c-b*b/(4*a)<=0:
            return [0,0]
        else:
            k=(e-b*d/(2*a))/np.sqrt(c-b*b/(4*a))
            return [d,k]


    def tf_response(self,t):
        tf_ls=self.tf_divide()
        
        y=np.zeros(t.shape)

        for tf in tf_ls:
            roots=tf.den.root_dict
            for key,val in roots.items():
                if key[1]==0:
                    sigma=key[0]
                    y += np.exp(sigma*t)*(t**(val-1))/(fact(val))
                elif key[1]>0:
                    cs=tf.tf_2nd_divide()
                    c=cs[0]
                    s=cs[1]
                    sigma=key[0]
                    w_n=abs(key[1])
                    y += np.exp(sigma*t)*(c*np.cos(w_n*t)+s*np.sin(w_n*t))*(t**(val-1))/(fact(val))
                else:
                    continue
        
        return y
    
    def tf_print(self):
        t=""
        t += self.num.tfalg_print()
        t += " / "
        t += self.den.tfalg_print()
        print(t)
        return