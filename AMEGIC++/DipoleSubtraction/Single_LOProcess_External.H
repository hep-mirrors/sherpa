#ifndef AMEGIC_DipoleSubtraction_Single_LOProcess_External_H
#define AMEGIC_DipoleSubtraction_Single_LOProcess_External_H

#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "AMEGIC++/Amplitude/Zfunctions/Mom.H"

namespace AMEGIC {
  class FullAmplitude_External;
  class Single_LOProcess_External : public Single_LOProcess {
  private:
    bool                    m_emitgluon;

    std::vector<ATOOLS::Vec4D>* p_epol;

    FullAmplitude_External    * p_MHVamp;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

   Single_LOProcess_External(const PHASIC::Process_Info &pi,
                        BEAM::Beam_Spectra_Handler *const beam,
                        PDF::ISR_Handler *const isr);
    ~Single_LOProcess_External();


    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  private :
    int         Tests(std::vector<double>* =NULL);
  public:
    int         InitAmplitude(MODEL::Model_Base *,Topology *,
			      std::vector<Process_Base *> &,
			      std::vector<Process_Base *> &);
    int         InitAmplitude(MODEL::Model_Base *,Topology *,
			      std::vector<Process_Base *> &,
			      std::vector<Process_Base *> &,
			      std::vector<ATOOLS::Vec4D>*,std::vector<double>*);

    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double  operator()(const ATOOLS::Vec4D_Vector &,const ATOOLS::Vec4D *,
		       std::vector<double>*,std::vector<ATOOLS::Vec4D>*,const int);
    void    Calc_AllXS(const ATOOLS::Vec4D_Vector &,const ATOOLS::Vec4D *,
		       std::vector<std::vector<double> > &dsij,const int);
    int NumberOfDiagrams();
  };
}



#endif

