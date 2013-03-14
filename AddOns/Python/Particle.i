//%module Particle
%{
#include <ATOOLS/Math/Vector.H>
#include <ATOOLS/Phys/Flavour.H>
#include <ATOOLS/Phys/Flow.H>
#include <ATOOLS/Org/MyStrStream.H>
#include <iostream>
#include <list>
%}

namespace ATOOLS {

  class Blob;
   
  class Particle {
  public:

    Particle();
    Particle(const Particle & );
    Particle(int, Flavour=Flavour(kf_none), Vec4D=Vec4D(0.,0.,0.,0.),
	     char='a');
    ~Particle();

    double            ProperTime();
    double            LifeTime();
    Vec3D             Distance(double = -1.);
    const Vec4D&      Momentum() const;
    double            E() const;
    double            FinalMass() const;
    %extend {
      int Stat() const
      {
	return self->Status();
      }
    };
    char              Info() const;
    Vec4D             XProd() const;
    %extend {
      bool HasProdBlob() const
      {
	return self->ProductionBlob();
      }
      const Blob &ProdBlob() const
      {
	return *self->ProductionBlob();
      }
    };
    Vec4D             XDec() const;
    %extend {
      bool HasDecBlob() const
      {
	return self->DecayBlob();
      }
      const Blob &DecBlob() const
      {
	return *self->DecayBlob();
      }
    };
    double            Time() const ;
    Flavour           Flav() const ;
    const Flavour   & RefFlav() const;
    unsigned int      GetFlow( unsigned int ) const ;
    int               Number() const;
    int               Beam() const;

    %extend {
      std::string __str__() {
	MyStrStream conv;
	conv<<*self;
	return conv.str();
      };
    };

  };

}
