#ifndef SHRIMPS_Eikonals_Eikonal_Weights_H
#define SHRIMPS_Eikonals_Eikonal_Weights_H

#include "SHRiMPS/Eikonals/Rapidity_Density.H"
#include "SHRiMPS/Event_Generation/Ladder.H"

namespace SHRIMPS {
  class Eikonal_Weights {
  private:
    double m_lambda, m_Delta;
    double m_originalY, m_Ymax, m_bmax;
    double m_singletwt;

    Rapidity_Density m_density;
    
    double m_b1, m_b2;
  public:
    Eikonal_Weights();
    ~Eikonal_Weights() {}

    void SetEikonal(Omega_ik * eikonal);
    void SetImpactParameters(const double & b1, const double & b2);
    void AddRapidities(Ladder * ladder,const double & ymin,const double & ymax);
    colour_type::code PropColour(const double & y1,const double & y2);
    double WeightSingletOverOctet(const double & y1,const double & y2);

    void Test(const std::string & dirname);


    // double EmissionWeight(const double & b1,const double & b2,const double & y,
    // 			  const double & sup=1.);
    // double SingletWeight(const double & b1,const double & b2,
    // 			 const double & y1,const double & y2,
    // 			 const double & sup=1.,const int & nbeam=1);
    // double OctetWeight(const double & b1,const double & b2,
    // 		       const double & y1,const double & y2,
    // 		       const double & sup=1.,const int & nbeam=1);
    // double RescatterProbability(const double & b1,const double & b2,
    // 				const double & y1,const double & y2,
    // 				const double & sup=1.,const int & nbeam=0);
    // double EffectiveIntercept(double b1,double b2,const double & y);
    //double Sum(const double & b1,const double & b2,const double & y);
  };
}

#endif
