#This class stores only the particle information required for the reweighting procedure

import Flavour
import Vec4

class Particle:
    
    def __init__(self, a_InputLine):

        PropList = a_InputLine.strip('\n').split()
        if len(PropList) < 7:
            raise RuntimeError('Failed to read particle info from line "{0}", which contains less than 7 numbers'.format(a_InputLine.strip('\n')))
        try:
            self.Barcode  = int(PropList[1])
            self.PDGid    = int(PropList[2])
            anti          = True if self.PDGid < 0 else False
            self.Flavour  = Flavour.Flavour( abs(self.PDGid), anti )
            self.px       = float(PropList[3])
            self.py       = float(PropList[4])
            self.pz       = float(PropList[5])
            self.E        = float(PropList[6])
            self.FourMom  = Vec4.Vec4D(self.E, self.px, self.py, self.pz)
            self.VertIn  = int(PropList[11])

        except ValueError as Error:
            raise RuntimeError('Failed to read particle info from line "{0}"'.format(a_InputLine.strip('\n')))

    def __str__(self):
        return 'Particle with barcode {0}, PDGid {1}, px {2:.2e}, py {3:.2e}, pz {4:.2e}, energy {5:.2e} and VertIn {6}'.format( self.Barcode, self.PDGid, self.px, self.py, self.pz, self.E, self.VertIn)

    def __repr__(self):
        return self.__str__()
 
# barcode pdg_id px py pz E generated_mass status theta phi 

