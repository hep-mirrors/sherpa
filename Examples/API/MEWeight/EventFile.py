import Event
import Vertex
import Particle
import HardProcess

class EventFile:

    def __init__(self, a_FileName):
        self.Path          = a_FileName
        self.EventPositions = []
        self.NumEvents     = 0
        self.Infile        = open(self.Path, 'r')
        self.HardProcess   = HardProcess.HardProcess()
        counter = 0
        # count events and store line numbers (EventPositions)
        for Line in self.Infile:
            counter += 1
            if Line.startswith('E '):
                # looks like an Event Line, store the line number
                self.NumEvents+=1
                self.EventPositions.append(counter)
        print 'Found {0} Events in File {1}'.format(self.NumEvents, self.Path)
        

    # set the file iterator to the a_EventNr'th Event line such that a readline afterwards gives the corresponding event line
    def SetIterator(self, a_EventNr):
        if a_EventNr<0 or not a_EventNr%1==0:
            raise ValueError('Argument of "SetIterator" method must be of integer type')        
        if a_EventNr > self.NumEvents:
            raise ValueError('File contains only {0} events, iterator cannot be set to event number {1}'.format(self.NumEvents, a_EventNr))
        # determine the linenumber of the determina_EventNr'th Event
        self.Infile.seek(0)
        for i in range(1, self.EventPositions[a_EventNr-1]):
            self.Infile.readline()
            
    def GetEvent(self, a_EventNr):
        if a_EventNr<0 or not a_EventNr%1==0:
            raise ValueError('Argument of "GetEvent" method must be of integer type')        
        if a_EventNr > self.NumEvents:
            raise ValueError('File contains only {0} events, cant extract event number {1}'.format(self.NumEvents, a_EventNr))
        self.SetIterator(a_EventNr)

        # now read event metadata line and set up an Event return value
        EventInputLine = self.Infile.readline().strip('\n')
        r_Event = Event.Event(EventInputLine)

        # Read weight info line. None of this information is needed so we only check if the next line starts with 'N'.
        # That the subsequent line must be a weight info line starting with an "N " is only an assumption that could not 
        # yet be verified. The same applies to the units line, the cross section line and the PDF info line.
        NInputLine = self.Infile.readline().strip('\n')
        if not NInputLine.startswith('N '):
            raise RuntimeError('Expected weight information starting with "N " but read line "{0}"'.format(NInputLine.strip))
        UnitsInputLine = self.Infile.readline().strip('\n')
        if not UnitsInputLine.startswith('U '):
            raise RuntimeError('Expected units information line starting with "U " but read "{0}"'.format(UnitsInputLine))
        XSInputLine = self.Infile.readline().strip('\n')
        if not XSInputLine.startswith('C '):
            raise RuntimeError('Expected cross section information line starting with "C " but read "{0}"'.format(XSInputLine))
        PDFInputLine = self.Infile.readline().strip('\n')
        if not PDFInputLine.startswith('F '):
            raise RuntimeError('Expected PDF information line starting with "F " but read "{0}"'.format(PDFInputLine))

        # set up a hard process. While parsing particle lines, particles may be added to this hard process instance
        # before the hard process vertex is identified
        r_HardProcess = HardProcess.HardProcess()
        SigProcVertFound = False
        # now read vertex lines
        for VertCount in range( r_Event.NumVert ):
            Line = self.Infile.readline().strip('\n')
            # Vertex lines must begin with a "V "
            if not Line.startswith('V '):
                raise RuntimeError('Expected vertex line starting with "V " but read "{0}"'.format(Line))
            r_Vertex = Vertex.Vertex(Line)
            print 'Created {0}'.format(r_Vertex)

            if r_Vertex.Identifier == r_Event.SigProcVertIdentifier:
                if SigProcVertFound == True:
                    raise RuntimeError('More than one signal process vertex found in event {0}'.format(a_EventNr))
                r_Event.SigProcVert = r_Vertex
                SigProcVertFound = True
            # read particle lines
            for PartCount in range(r_Vertex.NOrphansIn+r_Vertex.NPartOut):
                Line = self.Infile.readline().strip('\n')
                # particle lines must start with a "P "
                if not Line.startswith('P '):
                    raise RuntimeError('Expected particle line starting with "P " but read "{0}"'.format(Line))
                r_Part = Particle.Particle(Line)
                # if the particle an incoming one of the hard process, add it to the hard process
                if r_Part.VertIn == r_Event.SigProcVertIdentifier:
                    r_HardProcess.InParticles.append(r_Part)
                r_Vertex.AddParticle(r_Part)
                print 'Added {0}'.format(r_Part)
        # if no signal process vertex is found in the event, something might be wrong
        if not SigProcVertFound:
            raise RuntimeError('No signal process vertex found for event number {0}'.format(r_Event.EventNum))
        # finally fill all necessary member variables of the hard process and return the event
        # (the outgoing particles of the hard process are the particles belongin to the signal process vertex)
        r_HardProcess.OutParticles.extend(r_Event.SigProcVert.Particles)
        r_Event.HardProcess = r_HardProcess
        return r_Event        
