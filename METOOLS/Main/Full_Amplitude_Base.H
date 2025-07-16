#ifndef METOOLS_Main_Full_Amplitude_Base_H
#define METOOLS_Main_Full_Amplitude_Base_H

#include "METOOLS/Main/Spin_Structure.H"

namespace METOOLS {
  class Partial_Amplitude_Base;

  class Full_Amplitude_Base : public Spin_Structure<std::vector<Complex> > {
  protected:
    ATOOLS::Flavour       * p_flavs;
    const ATOOLS::CMatrix * p_colormatrix;
    
    std::vector<Partial_Amplitude_Base * > m_amplitudes;

  public:
    Full_Amplitude_Base(ATOOLS::Flavour* flavs,size_t size);
    ~Full_Amplitude_Base();
    inline void AddPartialAmplitude(Partial_Amplitude_Base * pa) {
      m_amplitudes.push_back(pa);
    }
    size_t NumberOfAmplitudes() const { return m_amplitudes.size(); }

    double SummedSquared(const ATOOLS::Vec4D *,bool anti=false);
    //double SummedSquared(const ATOOLS::Vec4D *,int hels);

    size_t ColorSize() const;
    const  ATOOLS::CMatrix* GetColorMatrix() const;
    void   SetColorMatrix(const ATOOLS::CMatrix* colormatrix);
  };

  inline const ATOOLS::CMatrix * Full_Amplitude_Base::GetColorMatrix() const { 
    return p_colormatrix; 
  }
  inline void Full_Amplitude_Base::SetColorMatrix(const ATOOLS::CMatrix* cm) { 
    p_colormatrix=cm; 
  }
  inline size_t Full_Amplitude_Base::ColorSize() const { return 0; }
}

#endif
