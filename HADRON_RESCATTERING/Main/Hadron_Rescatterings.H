#ifndef HADRON_RESCATTERING_Main_Hadron_Rescatterings_H
#define HADRON_RESCATTERING_Main_Hadron_Rescatterings_H

#include "HADRON_RESCATTERING/XSecs/BaryonBaryon.H"
#include "HADRON_RESCATTERING/XSecs/MesonMeson.H"
#include "HADRON_RESCATTERING/XSecs/BaryonMeson.H"
#include "HADRON_RESCATTERING/XSecs/HR_Parameters.H"

#include "ATOOLS/Phys/Blob_List.H"

namespace HADRON_RESCATTERING {
  class Hadron_Rescatterings {
  private:
    bool                m_on;
    long unsigned int   m_nfails;

    double              m_mm2mb;
    BaryonBaryon *      p_BB;
    MesonMeson *        p_MM;
    BaryonMeson *       p_BM;

 public:
    Hadron_Rescatterings(const bool & on);
    ~Hadron_Rescatterings();
    
    ATOOLS::Blob * operator()(ATOOLS::Particle * A,
			      ATOOLS::Particle * B,
			      const double & dist2);

    void Initialize();
    void Reset();
  };
}

#endif
