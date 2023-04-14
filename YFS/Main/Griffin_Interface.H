#ifndef Griffin_Interface_H
#define Griffin_Interface_H


#include "YFS/Griffin/SMval.h"
#include "YFS/Griffin/EWPOZ2.h"
#include "YFS/Griffin/xscnnlo.h"

#include "PHASIC++/Process/Process_Info.H"
#include "YFS/Main/YFS_Base.H"



namespace Griffin {

class Griffin_Interface : public YFS::YFS_Base {
public :

  // SMval *p_griffin;
  Griffin_Interface();
  ~Griffin_Interface();

  SMval griffin_input;  // convert masses from PDG values to complex pole scheme

  // FA_SMNNLO *FAi;
  // FA_SMNNLO *FAf;
  // SW_SMNNLO *SWi;
  // SW_SMNNLO *SWf;

  void PrintLogo(std::ostream &s);
  void RegisterDefaults();
  bool Initialize(const PHASIC::Process_Info& pi);


  // TODO: Clean up
  double EvaluateLoop(const ATOOLS::Vec4D_Vector& momenta);
  void EvaluateBorn(const ATOOLS::Vec4D_Vector& momenta);

  int m_mode, m_inital, m_final;


};


}

#endif
