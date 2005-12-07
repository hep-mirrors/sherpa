#ifndef Veto_Statistics_H
#define Veto_Statistics_H

#define USING__Veto_Info

#ifdef USING__Veto_Info

namespace APACIC {

  struct svc {
    
    enum code {
      no_veto  = 1,
      kt_veto  = 2,
      ang_veto = 4,
      jet_veto = 8,
      lj_veto  = 16,
      mlm_veto = 32,
      size     = 64
    };

  };// end of struct ssc 

}// end of namespace APACIC

#endif

#endif
