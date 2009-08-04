#ifndef METOOLS_Main_Three_Particle_Amplitudes_H
#define METOOLS_Main_Three_Particle_Amplitudes_H

#include "METOOLS/Main/Partial_Amplitude_Base.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/Main/XYZFuncs.H"

namespace METOOLS {

  class SSS : public Partial_Amplitude_Base {
  private:
  public:
    SSS(ATOOLS::Flavour *fl,int *i,bool *out);
    ~SSS() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

  class SFF : public Partial_Amplitude_Base {
  private:
    Complex            m_cL, m_cR;
    int                m_bar, m_nonbar;
    XYZFunc          * p_xyz;
  public:
    SFF(ATOOLS::Flavour *fl,int* i,bool *out,const Complex cl,const Complex cR);
    ~SFF() { if (p_xyz) delete p_xyz; }

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

  class SFF_FPI : public Partial_Amplitude_Base {
  private:
    Complex            m_c;
    int                m_bar, m_nonbar;
    XYZFunc          * p_xyz;
  public:
    SFF_FPI(ATOOLS::Flavour *fl,int *i,bool *out,const Complex c=Complex(1.0,0.0));
    ~SFF_FPI() { if (p_xyz) delete p_xyz; }

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

  class SSV : public Partial_Amplitude_Base {
  private:
  public:
    SSV(ATOOLS::Flavour *fl,int *i,bool *out);
    ~SSV() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

  class SVV : public Partial_Amplitude_Base {
  private:
  public:
    SVV(ATOOLS::Flavour *fl,int *i,bool *out);
    ~SVV() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

  class VFF : public Partial_Amplitude_Base {
  private:
    Complex            m_cL, m_cR;
    int                m_bar, m_nonbar;
    XYZFunc          * p_xyz;
  public:
    VFF(ATOOLS::Flavour *fl,int *i,bool *out,const Complex cl,const Complex cR);
    ~VFF() { if (p_xyz) delete p_xyz; }

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

  class VVV : public Partial_Amplitude_Base {
  private:
  public:
    VVV(ATOOLS::Flavour *fl,int *i,bool *out);
    ~VVV() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

  class TSS : public Partial_Amplitude_Base {
  private:
  public:
    TSS(ATOOLS::Flavour *fl,int *i,bool *out);
    ~TSS() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

  class TVS : public Partial_Amplitude_Base {
  private:
  public:
    TVS(ATOOLS::Flavour *fl,int *i,bool *out);
    ~TVS() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

  class TVV : public Partial_Amplitude_Base {
  private:
  public:
    TVV(ATOOLS::Flavour *fl,int *i,bool *out);
    ~TVV() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };
}

#endif
