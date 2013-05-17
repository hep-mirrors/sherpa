#ifndef AMEGIC_Amplitude_FullAmplitude_External_H
#define AMEGIC_Amplitude_FullAmplitude_External_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Permutation.H"

namespace AMEGIC {

  class Basic_Sfuncs;
  class Helicity;

  class FullAmplitude_External {
  protected:

    PHASIC::Tree_ME2_Base *p_calc;
    ATOOLS::Flavour_Vector m_fls;

    MODEL::Interaction_Model_Base *p_model;

    size_t m_nin, m_oqcd, m_oew, m_emit, m_spect;
    double m_A;
    Complex m_eip;

    std::vector<PHASIC::Tree_ME2_Base::Map_Info> m_amap;

    std::vector<std::vector<size_t> > m_hmap;
    std::vector<std::vector<int> >    m_pmap;

    std::vector<Complex> m_amps;

    std::vector<std::vector<std::vector<std::vector<Complex> > > > m_colfs;

    void GetPermutation(const std::vector<size_t> &ids,
			std::vector<size_t> &cid,
			ATOOLS::Flavour_Vector &cfl,
			int &nsub,int &psub,int &swap) const;

    void BuildHelicityMap(Helicity *const hel);
    void BuildColorMatrix(const size_t &ci,const size_t &cj);
    void BuildColorMatrix();

  public:
	
    // constructor
    FullAmplitude_External(const PHASIC::Process_Info &pi,
			   MODEL::Model_Base *model,
			   MODEL::Coupling_Map *const cpls,Helicity *hel,
			   const size_t &emit,const size_t &spect); 

    //destructor
    ~FullAmplitude_External(); 

    double Calc(const ATOOLS::Vec4D *p);

    double MSquare(const size_t &h) const;
    double MSquare(const size_t &h,const size_t &i,const size_t &j);
    
    void SetSqMatrix(const double &A,const ATOOLS::Vec4D &pijt,
		     const ATOOLS::Vec4D &eps1);

    inline int OrderQCD() const { return m_oqcd; }
    inline int OrderEW() const  { return m_oew;  }

    inline bool Status() const { return p_calc!=NULL; }

    inline PHASIC::Tree_ME2_Base *Calc() const { return p_calc; }

  };  

}// end of namespace AMEGIC

#endif
