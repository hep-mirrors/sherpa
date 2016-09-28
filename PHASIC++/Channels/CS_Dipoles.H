#ifndef PHASIC_Main_CS_Dipoles_H
#define PHASIC_Main_CS_Dipoles_H

#include "PHASIC++/Channels/CS_Dipole.H"

namespace PHASIC {

  class FF_Dipole: public CS_Dipole {
  private:
    double m_yexp, m_zexp, m_mi, m_mj, m_mk;
    double m_mi2, m_mj2, m_mij2, m_mk2;
    bool   m_massive;
  public:

    FF_Dipole(ATOOLS::NLO_subevt *const sub,
	      Phase_Space_Handler *const psh,const bool bmcw=0);

    ~FF_Dipole();

    ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,
     Cut_Data *const cuts,const double *rns);
    double GenerateWeight(const ATOOLS::Vec4D_Vector &p,
			  Cut_Data *const cuts);
    bool ValidPoint(const ATOOLS::Vec4D_Vector& p);

  };// end of class FF_Dipole

  class FI_Dipole: public CS_Dipole {
  private:
    double m_xmin, m_xexp, m_zexp;
    double m_mi, m_mj, m_mi2, m_mj2, m_mij2;
    bool   m_massive;
  public:

    FI_Dipole(ATOOLS::NLO_subevt *const sub,
	      Phase_Space_Handler *const psh,const bool bmcw=0);

    ~FI_Dipole();

    ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,
     Cut_Data *const cuts,const double *rns);
    double GenerateWeight(const ATOOLS::Vec4D_Vector &p,
			  Cut_Data *const cuts);
    bool ValidPoint(const ATOOLS::Vec4D_Vector& p);

  };// end of class FI_Dipole

  class IF_Dipole: public CS_Dipole {
  private:
    double m_xmin,m_xexp, m_uexp, m_mk2;
  public:

    IF_Dipole(ATOOLS::NLO_subevt *const sub,
	      Phase_Space_Handler *const psh,const bool bmcw=0);

    ~IF_Dipole();

    ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,
     Cut_Data *const cuts,const double *rns);
    double GenerateWeight(const ATOOLS::Vec4D_Vector &p,
			  Cut_Data *const cuts);
    bool ValidPoint(const ATOOLS::Vec4D_Vector& p);

  };// end of class IF_Dipole

  class II_Dipole: public CS_Dipole {
  private:
    double m_xmin, m_xexp, m_vexp;
  public:

    II_Dipole(ATOOLS::NLO_subevt *const sub,
	      Phase_Space_Handler *const psh,const bool bmcw=0);

    ~II_Dipole();

    ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,
     Cut_Data *const cuts,const double *rns);
    double GenerateWeight(const ATOOLS::Vec4D_Vector &p,
			  Cut_Data *const cuts);
    bool ValidPoint(const ATOOLS::Vec4D_Vector& p);

    void Calculate(const ATOOLS::Vec4D &pi, const ATOOLS::Vec4D &pj,
                   const ATOOLS::Vec4D &pk,
                   const ATOOLS::Vec4D_Vector& kj,
                   double &x, double &v, double &phi,
                   ATOOLS::Vec4D &pijt, ATOOLS::Vec4D &pkt,
                   ATOOLS::Vec4D_Vector& kjt);
    void Construct(ATOOLS::Vec4D &pi, ATOOLS::Vec4D &pj, ATOOLS::Vec4D &pk,
                   ATOOLS::Vec4D_Vector& kj,
                   const double &x, const double &v, const double &phi,
                   const ATOOLS::Vec4D &pijt, const ATOOLS::Vec4D &pkt,
                   const ATOOLS::Vec4D_Vector& kjt);

  };// end of class II_Dipole

}// end of namespace PHASIC

#endif
