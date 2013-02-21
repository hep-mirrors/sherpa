import Particle

class Vertex:

    def __init__(self, a_InputLine):

        PropList = a_InputLine.strip('\n').split()
        if (len(PropList)) < 10:
            raise RuntimeError('Incomplete vertex data line "{0}". Expected 10 entries but found {1}'.format(a_InputLine.strip('\n'), len(PropList)))

        try:
            self.Identifier   = int(PropList[1])
            self.ID           = int(PropList[2])
            self.x            = float(PropList[3])
            self.y            = float(PropList[4])
            self.z            = float(PropList[5])
            self.t            = float(PropList[6])
            self.NOrphansIn   = int(PropList[7])
            self.NPartOut     = int(PropList[8])
            self.WeightsSize  = int(PropList[9])
            self.Weights      = []
            self.Particles = []
            for i in range(10,10+self.WeightsSize):
                self.Weights.append(float(PropList[i]))
        except ValueError as Error:
            raise RuntimeError('Failed to read vertex data from line "{0}"'.format(a_InputLine.strip('\n')))
        if (len(PropList)) != (10+self.WeightsSize):
            raise RuntimeError('Inconsistent vertex data line "{0}". Expected {1} entries but found {2}'.format(a_InputLine.strip('\n'), 10+self.WeightsSize, len(PropList)))

    def AddParticle(self, a_Part):
        if (self.NPartOut + self.NOrphansIn) <= len(self.Particles):
            raise RuntimeError('Cannot add further outgoing particles to this vertex with NPartOut {0}, NOrphansIn {1}, and already {2} particles attached.'.format(self.NPartOut, self.NOrphansIn, len(self.Particles)))
        self.Particles.append(a_Part)
                
    def __str__(self):
        return 'Vertex with identifier {0}, ID {1}, x {2:.2e}, y {3:.2e}, z {4:.2e}, t {5:.2e}, NOrphansIn {6}, NPartOut {7}, WeightsSize {8}, Weights {9}'.format(self.Identifier, self.ID, self.x, self.y, self.z, self.t, self.NOrphansIn, self.NPartOut, self.WeightsSize, self.Weights)
