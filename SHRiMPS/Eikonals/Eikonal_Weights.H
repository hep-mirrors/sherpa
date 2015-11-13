#ifndef SHRIMPS_Eikonals_Eikonal_Weights_H
#define SHRIMPS_Eikonals_Eikonal_Weights_H

#include "SHRiMPS/Eikonals/Eikonal_Contributor.H"

namespace SHRIMPS {
  class Eikonal_Weights {
  private:
    Eikonal_Contributor * p_Omegaik, * p_Omegaki;
    double m_originalY, m_Ymax, m_bmax;
    double m_lambda, m_Delta;
    absorption::code m_absorp;
    double m_singletwt;

    double DeltaOmega(const double & b1,const double & b2,
		      const double & y1,const double & y2,
		      const double & sup=1.,const int & nbeam=1);
  public:
    Eikonal_Weights(const Eikonal_Parameters & params);
    ~Eikonal_Weights() {}

    void   SetSingleOmegaTerms(Eikonal_Contributor * omegaik, 
			       Eikonal_Contributor * omegaki);
    double MaximalEmissionProbability(const double & b1,const double & b2);
    double EmissionWeight(const double & b1,const double & b2,const double & y,
			  const double & sup=1.);
    double SingletWeight(const double & b1,const double & b2,
			 const double & y1,const double & y2,
			 const double & sup=1.,const int & nbeam=1);
    double OctetWeight(const double & b1,const double & b2,
		       const double & y1,const double & y2,
		       const double & sup=1.,const int & nbeam=1);
    double RescatterProbability(const double & b1,const double & b2,
				const double & y1,const double & y2,
				const double & sup=1.,const int & nbeam=0);
    double EffectiveIntercept(double b1,double b2,const double & y);
    //double Sum(const double & b1,const double & b2,const double & y);
  };
}

#endif
