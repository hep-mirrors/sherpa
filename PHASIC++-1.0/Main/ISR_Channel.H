#ifndef ISR_Channel_h
#define ISR_Channel_h

#include "Single_Channel.H"


/*!
  At the moment we assume explicitly a symmetric y-range (see deltay, given
  from outside). For assymmtreic collision we'll have to work this out via
  assymmetric deltay's, i.e. 2 deltay's !!!
  Note further that we played dirty in the backward peaks !
  Maybe we'll need an extra PeakedTralala method set.
*/


namespace PHASIC {
  class SimplePoleUniform : public Single_Channel {
  protected :
    double sprimeexp;
    double deltay[2];
  public :
    SimplePoleUniform(double _sprimeexp,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class SimplePoleCentral : public Single_Channel {
  protected :
    double sprimeexp;
    double deltay[2];
  public :
    SimplePoleCentral(double _sprimeexp,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class SimplePoleForward : public Single_Channel {
  protected :
    double norm,sprimeexp,yexp;
    double deltay[2];
  public :
    SimplePoleForward(double _sprimeexp,double _yexp,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class SimplePoleBackward : public Single_Channel {
  protected :
    double norm,sprimeexp,yexp;
    double deltay[2];
  public :
    SimplePoleBackward(double _sprimeexp,double _yexp,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };






  class ResonanceUniform : public Single_Channel {
  protected :
    double mass,width;
    double deltay[2];
  public :
    ResonanceUniform(double _mass,double _width,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class ResonanceCentral : public Single_Channel {
  protected :
    double mass,width;
    double deltay[2];
  public :
    ResonanceCentral(double _mass,double _width,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class ResonanceForward : public Single_Channel {
  protected :
    double norm,mass,width,yexp;
    double deltay[2];
  public :
    ResonanceForward(double _mass,double _width,double _yexp,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };


  class ResonanceBackward : public Single_Channel {
  protected :
    double norm,mass,width,yexp;
    double deltay[2];
  public :
    ResonanceBackward(double _mass,double _width,double _yexp,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };



  class ThresholdUniform : public Single_Channel {
  protected :
    double mass;
    double deltay[2];
  public :
    ThresholdUniform(double _mass,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class ThresholdCentral : public Single_Channel {
  protected :
    double mass;
    double deltay[2];
  public :
    ThresholdCentral(double _mass,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class ThresholdForward : public Single_Channel {
  protected :
    double norm,mass,yexp;
    double deltay[2];
  public :
    ThresholdForward(double _mass,double _yexp,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };


  class ThresholdBackward : public Single_Channel {
  protected :
    double norm,mass,yexp;
    double deltay[2];
  public :
    ThresholdBackward(double _mass,double _yexp,double _deltay1, double _deltay2);
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };



  class LLUniform : public Single_Channel {
  protected :
    double beta,factor;
    double deltay[2];
  public :
    LLUniform(double _beta,double _deltay1, double _deltay2) {
      name     = std::string("LLUniform");
      beta     = _beta;
      deltay[0] = _deltay1;
      deltay[1] = _deltay2;
      factor   = 1.00000001;
      ms = rans = 0;
    };
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class LLCentral : public Single_Channel {
  protected :
    double beta,factor;
    double deltay[2];
  public :
    LLCentral(double _beta,double _deltay1, double _deltay2) {
      name     = std::string("LLCentral");
      beta     = _beta;
      deltay[0] = _deltay1;
      deltay[1] = _deltay2;
      factor   = 1.00000001;
      ms = rans = 0;
    };
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };



  class LLForward : public Single_Channel {
  protected :
    double norm,yexp;
    double deltay[2];
    double beta,factor;
  public :
    LLForward(double _beta,double _yexp,double _deltay1, double _deltay2) {
      name     = std::string("LLForward");
      deltay[0] = _deltay1;
      deltay[1] = _deltay2;
      beta     = _beta;
      yexp     = _yexp;
      factor   = 1.00000001;
      ms = rans = 0;
    };
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class LLBackward : public Single_Channel {
  protected :
    double norm,yexp;
    double deltay[2];
    double beta,factor;
  public :
    LLBackward(double _beta,double _yexp,double _deltay1, double _deltay2) {
      name     = std::string("LLBackward");
      deltay[0] = _deltay1;
      deltay[1] = _deltay2;
      beta     = _beta;
      yexp     = _yexp;
      factor   = 1.00000001;
      ms = rans = 0;
    };
    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };



  class LBSComptonPeakUniform : public Single_Channel {
  protected :
    double pole,sexp;
    double deltay[2];
  public :
    LBSComptonPeakUniform(double _pole,double _sexp,double _deltay1, double _deltay2) { 
      name      = std::string("LBSComptonPeakUniform");
      pole      = _pole;
      sexp      = _sexp;
      deltay[0] = _deltay1;
      deltay[1] = _deltay2;
      ms = rans = 0;
    }; 

    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class LBSComptonPeakCentral : public Single_Channel {
  protected :
    double pole,sexp;
    double deltay[2];
  public :
    LBSComptonPeakCentral(double _pole,double _sexp,double _deltay1, double _deltay2) { 
      name      = std::string("LBSComptonPeakCentral");
      pole      = _pole;
      sexp      = _sexp;
      deltay[0] = _deltay1;
      deltay[1] = _deltay2;
      ms = rans = 0;
    }; 

    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class LBSComptonPeakForward : public Single_Channel {
  protected :
    double pole,sexp,yexp;
    double deltay[2];
  public :
    LBSComptonPeakForward(double _pole,double _sexp,double _deltay1, double _deltay2,double _yexp) { 
      name      = std::string("LBSComptonPeakForward");
      pole      = _pole;
      sexp      = _sexp;
      deltay[0] = _deltay1;
      deltay[1] = _deltay2;
      yexp     = _yexp;
      ms = rans = 0;
    }; 

    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };

  class LBSComptonPeakBackward : public Single_Channel {
  protected :
    double pole,sexp,yexp;
    double deltay[2];
  public :
    LBSComptonPeakBackward(double _pole,double _sexp,double _deltay1, double _deltay2,double _yexp) { 
      name      = std::string("LBSComptonPeakBackward");
      pole      = _pole;
      sexp      = _sexp;
      deltay[0] = _deltay1;
      deltay[1] = _deltay2;
      yexp     = _yexp;
      ms = rans = 0;
    }; 

    void GenerateWeight(double,double,int); 
    void GeneratePoint(double &,double &,int,double *);
  };


}


#endif
