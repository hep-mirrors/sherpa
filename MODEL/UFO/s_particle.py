#from parameters import ZERO
from s_parameter import s_parameter
from ufo_exception import ufo_exception

class s_particle:
    
    def __init__(self,ufo_particle):
        self.ufo_particle = ufo_particle

    def kf_code(self):
        return self.ufo_particle.pdg_code

    def mass(self):
        return self.ufo_particle.mass

    def width(self):
        pass
        
    def charge_times_three(self):
        # el. charge is not an integer for ufo particles
        ret = int(3*(self.ufo_particle.charge))
        assert isinstance(ret, int)
        return ret

    def mass(self):
        return s_parameter(self.ufo_particle.mass)

    def width(self):
        return s_parameter(self.ufo_particle.width)

    def icharge(self):
        pass

    # sherpa wants a zero for singlets
    def color(self):
        return 0 if self.ufo_particle.color==1 else self.ufo_particle.color

    def spin_times_two(self):
        ret = self.ufo_particle.spin-1
        assert isinstance(ret, int)
        return ret 

    # UFO doesn't comply with any convention here, 
    # have to assume case-insensitively that presence of a 'goldstone'
    # attribute signals a goldstone except for the case where such an 
    # attribute is zero
    def is_goldstone(self):
        gold_attr = [att for att in self.ufo_particle.__dict__ if att.lower().startswith('gold')]
        if len(gold_attr) == 0:
            return False
        falses = [0, 0.0, False]
        return any([ self.ufo_particle.__dict__[attr] not in falses for attr in gold_attr ])

    def majorana(self):
        # ufo_particle.spin is (2s+1) with spin s
        if self.spin_times_two() != 1:
            return 0
        if self.ufo_particle.name == self.ufo_particle.antiname:
            return 1
        return 0

    def self_conjugate(self):
        # for sherpa compatibility: 1 if majorana
        #                           0 if not self conjugate
        #                          -1 if self conj but not majorana
        if self.ufo_particle.name != self.ufo_particle.antiname:
            return 0
        if self.majorana():
            return 1
        return -1
        
    def name(self):
        return self.ufo_particle.name

    def antiname(self):
        return self.ufo_particle.antiname

    def texname(self):
        return self.ufo_particle.texname

    def antitexname(self):
        return self.ufo_particle.antitexname

    def texname(self):
        # replace backslahes by escaped backslashes
        return self.ufo_particle.texname.replace(r"\ ".rstrip(), r"\\")
