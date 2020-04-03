#ifndef EXTRA_XS_Main_Single_Process_H
#define EXTRA_XS_Main_Single_Process_H

#include "PHASIC++/Process/Single_Process.H"

namespace PHASIC { class Virtual_ME2_Base; }

namespace EXTRAXS {

  class ME2_Base;

  class Single_Process: public PHASIC::Single_Process {
  private:
    ME2_Base                              * p_born_me2;
    PHASIC::Virtual_ME2_Base              * p_virtual_me2;
    ATOOLS::nlo_type::code                  m_nlotype;
    std::set<std::pair<size_t,size_t> >     m_ccombs;
    std::map<size_t,ATOOLS::Flavour_Vector> m_cfls;
    
    bool m_localFS;
  protected:
    double m_sprimemin, m_sprimemax;
  public:

    // constructors
    Single_Process();

    // destructor
    ~Single_Process();

    // member functions
    bool Initialize();
    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);

    double Partonic(const ATOOLS::Vec4D_Vector& momenta,int mode);

    bool Combinable(const size_t &idi, const size_t &idj);

    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);
    bool FillFinalState(const ATOOLS::Vec4D_Vector &p);
    bool FillResonances(ATOOLS::Flavour_Vector & flavs);

    
    // inline functions
    inline ME2_Base *GetME() const { return p_born_me2; }
    //size_t NOut() const;
    //const ATOOLS::Flavour_Vector   & Flavours() const;
    //const ATOOLS::Vec4D_Vector     & Momenta() const;
    std::vector<std::vector<int> > * Colours() const;
    const bool   HasInternalScale() const;
    const double InternalScale() const;
    inline const double SPrimeMin() const { return m_sprimemin; }
    inline const double SPrimeMax() const { return m_sprimemax; }
  };

}// end of namespace EXTRAXS

#endif



