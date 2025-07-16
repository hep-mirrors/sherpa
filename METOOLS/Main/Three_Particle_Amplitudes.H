#ifndef METOOLS_Main_Three_Particle_Amplitudes_H
#define METOOLS_Main_Three_Particle_Amplitudes_H

#include "METOOLS/Main/Partial_Amplitude_Base.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/Main/XYZFuncs.H"

namespace METOOLS {

  class SSS : public Partial_Amplitude_Base {
  private:
  public:
    SSS(const ATOOLS::Flavour_Vector& flavs,
        const std::vector<int>& i, const std::vector<bool>& out);
    ~SSS() {}

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

  class SFF : public Partial_Amplitude_Base {
  private:
    Complex            m_cL, m_cR;
    int                m_bar, m_nonbar;
    XYZFunc          * p_xyz;
  public:
    SFF(const ATOOLS::Flavour_Vector& flavs,
        const std::vector<int>& i, const std::vector<bool>& out,
        const Complex cl,const Complex cR);
    ~SFF() { if (p_xyz) delete p_xyz; }

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

  class SFF_FPI : public Partial_Amplitude_Base {
  private:
    Complex            m_c;
    int                m_bar, m_nonbar;
    XYZFunc          * p_xyz;
  public:
    SFF_FPI(const ATOOLS::Flavour_Vector& flavs,
            const std::vector<int>& i, const std::vector<bool>& out,
            const Complex c=Complex(1.0,0.0));
    ~SFF_FPI() { if (p_xyz) delete p_xyz; }

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

  class SSV : public Partial_Amplitude_Base {
  private:
  public:
    SSV(const ATOOLS::Flavour_Vector& flavs,
        const std::vector<int>& i, const std::vector<bool>& out);
    ~SSV() {}

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

  class SVV : public Partial_Amplitude_Base {
  private:
  public:
    SVV(const ATOOLS::Flavour_Vector& flavs,
        const std::vector<int>& i, const std::vector<bool>& out);
    ~SVV() {}

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

  class VFF : public Partial_Amplitude_Base {
  private:
    Complex            m_cL, m_cR;
    int                m_bar, m_nonbar;
    XYZFunc          * p_xyz;
  public:
    VFF(const ATOOLS::Flavour_Vector& flavs,
        const std::vector<int>& i, const std::vector<bool>& out,
        const Complex cl,const Complex cR);
    ~VFF() { if (p_xyz) delete p_xyz; }

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

  class VVV : public Partial_Amplitude_Base {
  private:
  public:
    VVV(const ATOOLS::Flavour_Vector& flavs,
        const std::vector<int>& i, const std::vector<bool>& out);
    ~VVV() {}

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

  class TSS : public Partial_Amplitude_Base {
  private:
  public:
    TSS(const ATOOLS::Flavour_Vector& flavs,
        const std::vector<int>& i, const std::vector<bool>& out);
    ~TSS() {}

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

  class TVS : public Partial_Amplitude_Base {
  private:
  public:
    TVS(const ATOOLS::Flavour_Vector& flavs,
        const std::vector<int>& i, const std::vector<bool>& out);
    ~TVS() {}

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

  class TVV : public Partial_Amplitude_Base {
  private:
  public:
    TVV(const ATOOLS::Flavour_Vector& flavs,
        const std::vector<int>& i, const std::vector<bool>& out);
    ~TVV() {}

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };
}

#endif
