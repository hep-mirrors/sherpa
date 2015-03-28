#ifndef AMEGIC_DipoleSubtraction_Single_LOProcess_MHV_H
#define AMEGIC_DipoleSubtraction_Single_LOProcess_MHV_H

#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "AMEGIC++/Amplitude/Zfunctions/Mom.H"

namespace AMEGIC {
  class FullAmplitude_MHV_Base;
  class Single_LOProcess_MHV : public Single_LOProcess {
  private:
    bool                    m_emitgluon;
    bool                    m_ownamps;

    std::vector<ATOOLS::Vec4D>* p_epol;

    FullAmplitude_MHV_Base    * p_MHVamp;
    MomentumList              * p_momlist;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

   Single_LOProcess_MHV(const PHASIC::Process_Info &pi,
                        BEAM::Beam_Spectra_Handler *const beam,
                        PDF::ISR_Handler *const isr);
    ~Single_LOProcess_MHV();


    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  private :
    int         Tests(std::vector<double>* =NULL);
  public:
    int         InitAmplitude(Amegic_Model *,Topology *,
			      std::vector<Process_Base *> &,
			      std::vector<Process_Base *> &,int);
    int         InitAmplitude(Amegic_Model *,Topology *,
			      std::vector<Process_Base *> &,
			      std::vector<Process_Base *> &,
			      std::vector<ATOOLS::Vec4D>*,std::vector<double>*);

    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  private:
    Complex CalculateHelicityPhase(const ATOOLS::Vec4D *);
  public:
    double  operator()(const ATOOLS::Vec4D_Vector &,const ATOOLS::Vec4D *,
		       std::vector<double>*,std::vector<ATOOLS::Vec4D>*,const int);
    void    Calc_AllXS(const ATOOLS::Vec4D_Vector &,
		       const ATOOLS::Vec4D *,std::vector<std::vector<double> > &,const int);
  };
}



#endif

