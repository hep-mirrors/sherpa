import Particle
import Vertex
import HardProcess

class Event:

    def __init__(self, a_InputLine):
        self.EventNum     = 0
        self.NumVert      = 0
        self.RandSize     = 0
        self.WeightsSize  = 0
        self.Weights      = []
        self.SigProcVertIdentifier = 0        
        self.SigProcVert  = Vertex.Vertex("V 0 0 0 0 0 0 0 0 0")
        self.HardProcess  = HardProcess.HardProcess()
        
        PropList = a_InputLine.strip('\n').split()
        print '\nReading event metadata from line "{0}"'.format(a_InputLine.strip('\n'))
        try:
            self.EventNum    = int(PropList[1])
            # It looks like SigProcVert-1 is the identifier of the hard process vertex 
            # in the Event. This is only an assumption and has not yet been verified!
            self.SigProcVertIdentifier = int(PropList[7])
            self.NumVert     = int(PropList[8])
            self.RandSize    = int(PropList[11])
            self.WeightsSize = int(PropList[int(12+self.RandSize)])
            for i in range(13+self.RandSize,13+self.RandSize+self.WeightsSize):
                self.Weights.append(float(PropList[i]))
        except ValueError as Error:
            raise RuntimeError('Failed to read event metadata from line "{0}"'.format(a_InputLine.strip('\n')))

        print 'Created instance of "Event" with event number {0}, signal process identifier {1}, number of vertices {2}, and weights {3}'.format(self.EventNum, self.SigProcVertIdentifier, self.NumVert, self.Weights)
  
