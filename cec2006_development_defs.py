from __future__ import annotations
import numpy as np
from cec2006_unseen_defs import P

class G1(P):
    def __init__(self):
        xl=np.zeros(13); xu=np.array([1,1,1,1,1,1,1,1,1,100,100,100,1],float); ps=[1,1,1,1,1,1,1,1,1,3,3,3,1]
        super().__init__('g1',13,xl,xu,ps)
    def _eval(self,x):
        f=5*np.sum(x[:,:4],axis=1)-5*np.sum(x[:,:4]**2,axis=1)-np.sum(x[:,4:13],axis=1)
        g1=2*x[:,0]+2*x[:,1]+x[:,9]+x[:,10]-10
        g2=2*x[:,0]+2*x[:,2]+x[:,9]+x[:,11]-10
        g3=2*x[:,1]+2*x[:,2]+x[:,10]+x[:,11]-10
        g4=-8*x[:,0]+x[:,9]; g5=-8*x[:,1]+x[:,10]; g6=-8*x[:,2]+x[:,11]
        g7=-2*x[:,3]-x[:,4]+x[:,9]; g8=-2*x[:,5]-x[:,6]+x[:,10]; g9=-2*x[:,7]-x[:,8]+x[:,11]
        return f,np.column_stack([g1,g2,g3,g4,g5,g6,g7,g8,g9]),None

class G3(P):
    def __init__(self):
        n=10; ps=np.full(n,1/np.sqrt(n)); super().__init__('g3',n,np.zeros(n),np.ones(n),ps)
    def _eval(self,x):
        f=-(np.sqrt(self.n_var)**self.n_var)*np.prod(x,axis=1); h=np.sum(x**2,axis=1)-1
        return f,None,h

class G4(P):
    def __init__(self): super().__init__('g4',5,[78,33,27,27,27],[102,45,45,45,45],[78,33,29.9952560256815985,45,36.7758129057882073])
    def _eval(self,x):
        f=5.3578547*x[:,2]**2+0.8356891*x[:,0]*x[:,4]+37.293239*x[:,0]-40792.141
        u=85.334407+0.0056858*x[:,1]*x[:,4]+0.0006262*x[:,0]*x[:,3]-0.0022053*x[:,2]*x[:,4]
        v=80.51249+0.0071317*x[:,1]*x[:,4]+0.0029955*x[:,0]*x[:,1]+0.0021813*x[:,2]**2
        w=9.300961+0.0047026*x[:,2]*x[:,4]+0.0012547*x[:,0]*x[:,2]+0.0019085*x[:,2]*x[:,3]
        return f,np.column_stack([-u,u-92,-v+90,v-110,-w+20,w-25]),None

class G6(P):
    def __init__(self): super().__init__('g6',2,[13,0],[100,100],[14.095,5-np.sqrt(100-(14.095-5)**2)])
    def _eval(self,x):
        f=(x[:,0]-10)**3+(x[:,1]-20)**3
        g1=-((x[:,0]-5)**2)-(x[:,1]-5)**2+100
        g2=(x[:,0]-6)**2+(x[:,1]-5)**2-82.81
        return f,np.column_stack([g1,g2]),None

class G8(P):
    def __init__(self): super().__init__('g8',2,[1e-5,1e-5],[10,10],[1.22797135260752599,4.24537336612274885])
    def _eval(self,x):
        f=-(np.sin(2*np.pi*x[:,0])**3*np.sin(2*np.pi*x[:,1]))/(x[:,0]**3*(x[:,0]+x[:,1]))
        g1=x[:,0]**2-x[:,1]+1; g2=1-x[:,0]+(x[:,1]-4)**2
        return f,np.column_stack([g1,g2]),None

class G9(P):
    def __init__(self):
        ps=[2.33049949323300210,1.95137239646596039,-0.47754041766198602,4.36572612852776931,-0.62448707583702823,1.03813092302119347,1.59422663221959926]
        super().__init__('g9',7,np.full(7,-10.0),np.full(7,10.0),ps)
    def _eval(self,x):
        f=(x[:,0]-10)**2+5*(x[:,1]-12)**2+x[:,2]**4+3*(x[:,3]-11)**2+10*x[:,4]**6+7*x[:,5]**2+x[:,6]**4-4*x[:,5]*x[:,6]-10*x[:,5]-8*x[:,6]
        v1=2*x[:,0]**2; v2=x[:,1]**2
        g1=v1+3*v2**2+x[:,2]+4*x[:,3]**2+5*x[:,4]-127
        g2=7*x[:,0]+3*x[:,1]+10*x[:,2]**2+x[:,3]-x[:,4]-282
        g3=23*x[:,0]+v2+6*x[:,5]**2-8*x[:,6]-196
        g4=2*v1+v2-3*x[:,0]*x[:,1]+2*x[:,2]**2+5*x[:,5]-11*x[:,6]
        return f,np.column_stack([g1,g2,g3,g4]),None

REG={p.name:p for p in [G1(),G3(),G4(),G6(),G8(),G9()]}
def get_problem(name): return REG[name]
