import Particle

class HardProcess:

    def __init__(self):
        self.InParticles  = []
        self.OutParticles = []


    def AddOutParticle(self, a_Particle):
        self.OutParticles.append(a_Particle)
        print 'OUTPARTICLE {0}'.format(a_Particle)

    def AddInParticle(self, a_Particle):
        self.InParticles.append(a_Particle)

    def __str__(self):
        OutIDs = []
        OutMomenta = []
        InIDs = []
        InMomenta = []
        for i in self.InParticles:
            InIDs.append(i.PDGid) 
            InMomenta.append( [i.px,i.py,i.pz,i.E] )
        for i in self.OutParticles:
            OutIDs.append(i.PDGid)
            OutMomenta.append( [i.px,i.py,i.pz,i.E] )
        return '{0}  --->  {1}'.format(InIDs, OutIDs)
        #return '{0}  --->  {1}\n{2}  ---> {3}'.format(InIDs, OutIDs, InMomenta, OutMomenta)
