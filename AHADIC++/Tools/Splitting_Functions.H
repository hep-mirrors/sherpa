#ifndef AHADIC_Tools_Splitting_Functions_H
#define AHADIC_Tools_Splitting_Functions_H

#include "AHADIC++/Tools/Strong_Coupling.H"
#include "AHADIC++/Tools/Cluster.H"
#include "AHADIC++/Tools/Dipole.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>

namespace AHADIC {
  struct ZForm {
    enum code {
      fragmentation    = 30,
      splitkernel_frag = 3,
      splitkernel_flat = 2,
      splitkernel      = 1,
      flat             = 0
    };
  };
  ZForm::code DefineZForm(const int & zform);




  class Splitting_Functions {
  private:
    ZForm::code m_zform;
    int         m_masstreatment;
    double      m_alpha,m_sigma,m_kappa;

    double m_m1_2,m_m2_2,m_m3_2,m_Q2,m_y,m_viji,m_vijk,m_zp,m_zm,m_pipj;

    double SelectZFromG2QQSplittingFunction(const double & zmin,const double & zmax);
    double SelectZFromQ2QGSplittingFunction(const double & zmin,const double & zmax);
    double SelectZFromQ2QGFragmentationFunction(const double & zmin,const double & zmax);

    double WeightForG2QQSplittingFunction(const double & z,const double & kt2);
    double WeightForQ2QGSplittingFunction(const double & z,const double & kt2);
    double WeightForQ2QGFragmentationFunction(const double & z,const double & kt2);

    double J(const double & z,const bool & glusplit);

    double IntegratedFromQ2QGSplittingFunction(const double & zmin,const double & zmax);
    double lambda(const double & x,const double & y,const double & z);
  public:
    Splitting_Functions(const ZForm::code & zform,const int & masstreatment);
    ~Splitting_Functions() {};
    double SelectZ(const double & zmin,const double & zmax,
		   const bool & glusplit,const bool & leading=true);
    double Weight(const double & kt2,const double & z,
		  const bool & glusplit,const bool & leading=true);
    double Integrated(const double & zmin,const double & zmax,const bool & glusplit);
    void   SetScales(const double & m1_2,const double & m2_2,const double & m3_2,
		     const double & Q2,const double & y,
		     const double & viji,const double & vijk,
		     const double & zp,const double & zm,const double & pipj) {
      m_m1_2 = m1_2; m_m2_2 = m2_2; m_m3_2 = m3_2;
      m_Q2 = Q2; m_y = y; m_viji = viji; m_vijk = vijk; 
      m_zp = zp; m_zm = zm; m_pipj = pipj;
    }

  };
}

#endif
