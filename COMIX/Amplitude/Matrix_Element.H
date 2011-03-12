#ifndef COMIX_Amplitude_Matrix_Element_H
#define COMIX_Amplitude_Matrix_Element_H

#include "COMIX/Amplitude/Amplitude.H"

namespace PHASIC {

  class Color_Integrator;
  class Helicity_Integrator;

}

namespace COMIX {

  class Matrix_Element {
  private:

    size_t m_nin, m_nout;

    Amplitude m_ampl;

    PHASIC::Color_Integrator    *p_colint;
    PHASIC::Helicity_Integrator *p_helint;

    double m_sf, m_fsf;

    int m_mode, m_tests;

  public:

    // constructor
    Matrix_Element();
    
    // destructor
    ~Matrix_Element();

    // member functions
    bool Initialize(const size_t &nin,const size_t &nout,
		    const std::vector<ATOOLS::Flavour> &flavs,
		    const double &isf,const double &fsf,
		    Model *const model,MODEL::Coupling_Map *const cpls,
		    const size_t &oew=99,const size_t &oqcd=99,
		    const size_t &maxoew=99,const size_t &maxoqcd=99,
		    const size_t &minntc=0);
    bool Map(const Matrix_Element *me,ATOOLS::Flavour_Map &fmap);

    double Differential(const std::vector<ATOOLS::Vec4D> &momenta,
			PHASIC::Color_Integrator *const colint,
			PHASIC::Helicity_Integrator *const helint);
    double Differential(const std::vector<ATOOLS::Vec4D> &momenta,
			const Int_Vector &ci,const Int_Vector &cj,
			const bool set=false);

    bool GaugeTest(std::vector<ATOOLS::Vec4D> momenta);

    void PrintGraphs(const std::string &gpath) const;

    void PrintStatistics(std::ostream &str,const int mode=0) const;

    // inline functions
    inline void SetMode(const int &mode)    { m_mode=mode;   }
    inline void SetTests(const int &tests)  { m_tests=tests; }

    inline void SetColorIntegrator(PHASIC::Color_Integrator *const cint)
    { p_colint=cint; }
    inline void SetHelicityIntegrator(PHASIC::Helicity_Integrator *const hint)
    { p_helint=hint; }

    inline PHASIC::Color_Integrator *ColorIntegrator() const
    { return p_colint; }
    inline PHASIC::Helicity_Integrator *HelicityIntegrator() const
    { return p_helint; }

    inline double FSSymmetryFactor() const { return m_fsf;      }
    inline double ISSymmetryFactor() const { return m_sf/m_fsf; }

    inline size_t MaxOrderEW() const  { return m_ampl.MaxOrderEW();  }
    inline size_t MaxOrderQCD() const { return m_ampl.MaxOrderQCD(); }

    inline Amplitude *GetAmplitude() { return &m_ampl; }

    inline void SetDecayInfos(const ATOOLS::DecayInfo_Vector &id) 
    { m_ampl.SetDecayInfos(id); }
    
  };// end of class Matrix_Element

}// end of namespace COMIX

#endif

