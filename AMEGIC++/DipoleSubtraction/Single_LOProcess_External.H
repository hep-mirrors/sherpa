#ifndef AMEGIC_DipoleSubtraction_Single_LOProcess_External_H
#define AMEGIC_DipoleSubtraction_Single_LOProcess_External_H

#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"

namespace AMEGIC {

  class Single_LOProcess_External : public Single_LOProcess {
  protected:

    PHASIC::Tree_ME2_Base *p_me2;

    size_t m_pn;

  public:

    Single_LOProcess_External(const PHASIC::Process_Info &pi,
			      BEAM::Beam_Spectra_Handler *const beam,
			      PDF::ISR_Handler *const isr);
    ~Single_LOProcess_External();

  public:

    void                AddChannels(std::list<std::string>*);
    bool                NewLibs() {return m_newlib;}
    bool                FoundMappingFile(std::string &,std::string &);
    virtual int         InitAmplitude(MODEL::Model_Base *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    void                WriteMomFlavs(std::ofstream& to) { p_BS->WriteMomFlavs(to); }
    double              Result()                       { return m_iresult; } 

    double              Norm()  {return m_Norm;}

    inline size_t       GetEmit()  { return m_emit; }
    inline size_t       GetSpect() { return m_spect; }
    int                 PerformTests();
    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    void             Minimize();
    int  Type() { return m_type; }
    inline const std::vector<size_t>& PartonList() const { return m_partonlist; }

    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    int            NumberOfDiagrams();
    double         Partonic(const ATOOLS::Vec4D_Vector&,const int);
    virtual double operator()(const ATOOLS::Vec4D_Vector &,const ATOOLS::Vec4D *,
			      std::vector<double>*,std::vector<ATOOLS::Vec4D>*,const int);
    virtual void   Calc_AllXS(const ATOOLS::Vec4D_Vector &,
			      const ATOOLS::Vec4D *,std::vector<std::vector<double> > &,const int);
  };
}



#endif

