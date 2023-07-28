#ifndef SHRIMPS_Glauber_Glauber_H
#define SHRIMPS_Glauber_Glauber_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Cross_Sections/Sigma_D.H"
#include "SHRiMPS/Cross_Sections/Sigma_DD.H"
#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_Total.H"
#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Phys/Particle.H"
#include <vector>
#include <map>
#include <cmath>

namespace SHRIMPS {
  class Glauber{
    public:
      int m_type_elas, m_type_abs, m_type_SD1, m_type_SD2, m_type_DD;
      int m_numNucleons1, m_numNucleons2;
      double m_radius1, m_radius2;
      double m_rho0, m_w, m_a;
      std::vector<ATOOLS::Vec4D> m_nucleus1;
      std::vector<ATOOLS::Vec4D> m_nucleus2;
      ATOOLS::Vec4D m_nucleus_1_position;
      ATOOLS::Vec4D m_nucleus_2_position;
      std::vector<std::pair<std::pair<int,int>,int>> m_list_of_interactions;
      bool m_pp, m_pA, m_AA;
      int m_projState;

      Glauber(Cross_Sections * xsecs, int N1 = 1, int N2 = 208);
      ~Glauber();
      ATOOLS::Vec4D distributeNucleon(double radius, ATOOLS::Vec4D nuc_pos = ATOOLS::Vec4D(0.,0.,0.,0.));
      double WoodsSaxon(double r, double Radius);
      void DoCollision(Cross_Sections * xsecs, std::vector<ATOOLS::Vec4D> pos_N1,std::vector<ATOOLS::Vec4D> pos_N2);
      void SaveXSs(Cross_Sections * xsecs,double bmin=0., double bmax=20., int num=100);
      void SaveNucleonPositions();
  };
}
#endif
