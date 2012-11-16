#!/usr/bin/python

class BeamBox:
    def __init__(self):
        self.BeamParticles = []
        self.fillBeamParticles()        
        self.Colliders = []
        self.fillColliders()        
        self.initialiseDefaults()        

    def initialiseDefaults(self):
        self.beam_id      = [0,0]
        self.beam_energy  = [0.,0.]
        self.setCollider(self.Colliders[0][0])

 
    def getCollider(self):
        return self.collider

    def getColliderName(self):
        return self.collider[1]

    def setCollider(self,coltag):
        for collider in self.Colliders:
            if collider[0]==coltag:
                self.collider       = collider
                self.collider_id    = self.collider[0]
                self.collider_tag   = self.collider[1]
                for i in range(0,2):
                    self.beam_id[i]     = self.collider[2+2*i]
                    self.beam_energy[i] = self.collider[3+2*i]

    def setBeam(self,beam_number,beam_id,beam_energy):
        self.beam_id[beam_number]      = beam_id
        self.beam_energy[beam_number]  = beam_energy
        self.collider[2+beam_number*2] = beam_id
        self.collider[3+beam_number*2] = beam_energy

    def getBeamParticles(self):
        return self.BeamParticles

    def getColliders(self):
        return self.Colliders
            
    def write(self,runfile):
        print "Write for ",self.collider[0],self.collider[1]
        tag      = self.collider[1]
        beam_id1 = self.collider[2]
        beam_en1 = self.collider[3]
        beam_id2 = self.collider[4]
        beam_en2 = self.collider[5]
        runfile.write('\n')
        runfile.write('  %%% Beam setup according to: '+tag+'\n\n')
        runfile.write('  BEAM_1 = %8s ; ' % str(beam_id1))
        runfile.write('BEAM_ENERGY_2 = %0.2f ;\n' % (beam_en1))
        runfile.write('  BEAM_2 = %8s ; ' % str(beam_id2))
        runfile.write('BEAM_ENERGY_2 = %0.2f ;\n' % (beam_en2))

    def printDefaults(self):
        print "Default settings in the beam box:"
        print "Collider = ",self.default_collider_id

    def printCollider(self):
        print "Collider: ",self.collider[0],self.collider[1]
        print "  Beam 1: ",self.collider[2],self.collider[3]
        print "  Beam 2: ",self.collider[4],self.collider[5]

    def fillColliders(self):
        self.Colliders.append([0,  "None",0,0.,0,0.])
        self.Colliders.append([123,"LHC pp (8000 GeV)",
                               2212,4000.,2212,4000.])
        self.Colliders.append([124,"LHC pp (13000 GeV)",
                               2212,6500.,2212,6500.])
        self.Colliders.append([125,"LHC pp (14000 GeV)",
                               2212,7000.,2212,7000.])
        self.Colliders.append([120,"LHC pp (900 GeV)",
                               2212,450.,2212,450.])
        self.Colliders.append([121,"LHC pp (2360 GeV)",
                               2212,1180.,2212,1180.])
        self.Colliders.append([122,"LHC pp (7000 GeV)",
                               2212,3500.,2212,3500.])
        self.Colliders.append([126,"LHC pp (33000 GeV)",
                               2212,16500.,2212,16500.])
        self.Colliders.append([1,  "LEP e+e- (91 GeV)",
                               11,45.6,-11,45.6])
        self.Colliders.append([2,  "ILC e+e- (500 GeV)",
                               11,250.,-11,250.])
        self.Colliders.append([110,"Tevatron ppbar Run 1 (1800 GeV)",
                               2212,900.,-2212,900.])
        self.Colliders.append([111,"Tevatron ppbar Run 2 (1960 GeV)",
                               2212,980.,-2212,980.])
        self.Colliders.append([-1, "User defined",0,0.,0,0.])

    def fillBeamParticles(self):
        self.BeamParticles.append([0,"None"])
        self.BeamParticles.append([2212,"proton"])
        self.BeamParticles.append([-2212,"anti-proton"])
        self.BeamParticles.append([11,"electron"])
        self.BeamParticles.append([-11,"positron"])
