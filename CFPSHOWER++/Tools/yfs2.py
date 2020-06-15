import math as m
import random as r

from vector import Vec4, NT, LT
from particle import Particle
from qcd import AlphaS, NC, TR, CA, CF
#from wdb import Solve

class Kernel:

    def __init__(self,flavs):
        self.flavs = flavs

class Sqq (Kernel):

    def Value(self,pi,pj,q,t):
        return CF*2.*t/(2.*pi*q)*(pi*pj)/(pi*q+pj*q)
    def Estimate(self,y,phi,k02):
        return CF*2.
    def Integral(self,k02):
        return CF*2.*m.log(1./k02)
    def GeneratePoint(self,k02):
        y = .5*m.log(1./k02)*(2.*r.random()-1.)
        phi = 2.*m.pi*r.random()
        return y, phi

class Sgg (Kernel):

    def Value(self,pi,pj,q,t):
        return CA/2.*2.*t/(2.*pi*q)*(pi*pj)/(pi*q+pj*q)
    def Estimate(self,y,phi,k02):
        return CA/2.*2.
    def Integral(self,k02):
        return CA/2.*2.*m.log(1./k02)
    def GeneratePoint(self,k02):
        y = .5*m.log(1./k02)*(2.*r.random()-1.)
        phi = 2.*m.pi*r.random()
        return y, phi

class Shower:

    def __init__(self,alpha,t0,coll):
        self.t0 = t0
        self.alpha = alpha
        self.alphamax = alpha(self.t0)
        self.kernels = [ Sqq([fl,fl,21]) for fl in [-5,-4,-3,-2,-1,1,2,3,4,5] ]
        self.kernels += [ Sgg([21,21,21]) ]
        if coll:
            self.kernels += [ Cqq([fl,fl,21]) for fl in [-5,-4,-3,-2,-1,1,2,3,4,5] ]
            self.kernels += [ Cgg([21,21,21]) ]
            self.kernels += [ Cgq([21,fl,-fl]) for fl in [1,2,3,4,5] ]

    def ESum(self,x,K,Q,qs):
        return m.sqrt(Q.M2()+(x*Q+K).P2())-K[0]- \
            sum([m.sqrt(q[0]+x*x*q[1]) for q in qs])
            
    def MakeKinematics(self,event,t,y,phi,pi,pj):
        nt = NT(pi,pj)
        lt = LT(pi,pj,nt)
        K = m.sqrt(t)*( (m.exp(y)*pi+m.exp(-y)*pj)/(pi+pj).M() \
                        +nt/nt.P()*m.cos(phi)+lt/m.sqrt(abs(lt.M2()))*m.sin(phi) )
        qq = Vec4()
        for p in event[2:]: qq += p.mom
#        u = 1.-(qq.E-m.sqrt(qq.M2()+K.P2())+K[0])/(pi.E+pj.E)
#        if u < 0.: return [],pi,pj,K,0.,0.
        A, B = (pi+pj)*(qq+K)/(pi+pj).M2(), 2.*qq.E*K[0]/(pi+pj).M2()
        if A*A < B: return [],pi,pj,K,0.,0.
        u = 1-(A-m.sqrt(A*A-B))
        sij = (pj+pi).M2()
        # require that k_{T,k}^{ij} < min(k_{T,j}^{ik},k_{T,i}^{jk})
        if u < 0. or abs(y) > .5*m.log(u*u*sij/t): return [],pi,pj,K,0.,0.
        pp = qq-(1.-u)*(pi+pj)
        pm = pp+K
#        print qq, pp, pm, pp.M2(), qq.M2()/pm.M2()
        w = pow(pp.M2()/qq.M2(),1.5)*qq[0]/(pp[0]+K[0])
        w *= pow(u,3.*2.-4.)#/pow(u,2.*(len(qs)-2.))
        #w *= qq.P2()/qq[0]-sum([q.P2()/q[0] for q in qs])
        #w /= qq[0]-(pp*qq)/pp[0]-sum([q[0]-(p*q)/p[0] for p,q in zip(ps,qs)])
        #for p,q in zip(ps,qs): w *= q[0]/p[0]
#        w *= qq.P2()/qq[0]-sum([q.P2()/q[0] for q in qs])
#        w /= qq[0]-(pp*qq)/pp[0]-sum([q[0]-(p*q)/p[0] for p,q in zip(ps,qs)])
        w *= 1./u**2
        Pi = pm.Boost(u*pi)
        Pj = pm.Boost(u*pj)
        K = pm.Boost(K)
        ps = []
        for p in event[2:]:
            if p.mom == pi: ps.append(Pi)
            elif p.mom == pj: ps.append(Pj)
            else: ps.append(pm.Boost(p.mom))
        return ps, Pi, Pj, K, w, u

    def MakeColors(self,flavs,colij,colk):
        self.c += 1
        if flavs[0] != 21:
            if flavs[0] > 0:
                return [ [self.c,0], [colij[0],self.c] ]
            else:
                return [ [0,self.c], [self.c,colij[1]] ]
        else:
            if flavs[1] == 21:
                if colij[0] == colk[1]:
                    if colij[1] == colk[0] and r.random()>0.5:
                        return [ [colij[0],self.c], [self.c,colij[1]] ]
                    return [ [self.c,colij[1]], [colij[0],self.c] ]
                else:
                    return [ [colij[0],self.c], [self.c,colij[1]] ]
            else:
                if flavs[1] > 0:
                    return [ [colij[0],0], [0,colij[1]] ]
                else:
                    return [ [0,colij[1]], [colij[0],0] ]

    def GeneratePoint(self,event):
        while self.t > self.t0:
            t = self.t0
            for split in event[2:]:
                for spect in event[2:]:
                    if spect == split: continue
                    for w in [1,-1]:
                        if not split.ColorConnected(spect,w): continue
                        for sf in self.kernels:
                            if sf.flavs[0] != split.pid: continue
                            m2 = (split.mom+spect.mom).M2()
                            if m2 < self.t0: continue
                            g = self.alphamax/(2.*m.pi)*sf.Integral(self.t0/m2)
                            tt = self.t*m.pow(r.random(),1./g)
                            if tt > t:
                                t = tt
                                s = [ split, spect, sf, m2 ]
            self.t = t
            if t > self.t0:
                y, phi = s[2].GeneratePoint(self.t0/s[3])
                p,pi,pj,q,v,u = self.MakeKinematics(event,t,y,phi,s[0].mom,s[1].mom)
                if v != 0.:
                    w = v*self.alpha(t)/self.alphamax
                    w *= s[2].Value(pi,pj,q,t)/s[2].Estimate(y,phi,self.t0/s[3])
                    if q.E > self.laste: continue
                    if w > r.random():
                        self.laste = q.E
                        cols = self.MakeColors(s[2].flavs,s[0].col,s[1].col)
                        for i in range(0,len(p)):
                            cp = event[2+i]
                            cp.Set(cp.pid,p[i],cp.col)
                        event.append(Particle(s[2].flavs[2],q,cols[1]))
                        s[0].Set(s[2].flavs[1],pi,cols[0])
                        s[1].mom = pj
                        return
    
    def Run(self,event,nem):
        em = 0
        self.c = 2
        self.t = (event[0].mom+event[1].mom).M2()
        self.laste = 1.e37
        while self.t > self.t0:
            if em >= nem: return
            self.GeneratePoint(event)
            em += 1
            
import sys, optparse
from matrix import eetojj
from durham import Analysis
from particle import CheckEvent

parser = optparse.OptionParser()
parser.add_option("-s","--seed",default=123456,dest="seed")
parser.add_option("-e","--events",default=1000,dest="events")
parser.add_option("-f","--file",default="dire",dest="histo")
parser.add_option("-c","--collinear",default=False,action="store_true",dest="coll")
parser.add_option("-n","--nem",default=1000000,dest="nem")
(opts,args) = parser.parse_args()

hardxs = eetojj()
shower = Shower(AlphaS(91.1876,0.118),1.,opts.coll)
jetrat = Analysis()

r.seed(int(opts.seed))
for i in range(int(opts.events)):
    event, weight = hardxs.GeneratePoint()
    shower.Run(event,int(opts.nem))
    if not CheckEvent(event):
        print(event)
    if i % 1000 == 0:
        sys.stdout.write('\rEvent {0}'.format(i))
    sys.stdout.flush()
    jetrat.Analyze(event,weight)
jetrat.Finalize(opts.histo)
print("")
