#ifndef HADRON_RESCATTERING_XSECS_HERAPARAM_H
#define HADRON_RESCATTERING_XSECS_HERAPARAM_H

#include <vector>
#include <map>


namespace HADRON_RESCATTERING 
{
  struct HeraParameters {
    enum code {
      pPiMinus  = 1,
      pPiPlus   = 2,
      nPiMinus  = 3,
      nPiPlus   = 4,
      pKMinus   = 5,
      pKBarZero = 6,
      nKMinus   = 7,
      nKBarZero = 8
    };
  };


class HeraParam
{
    double m1,m2;
    std::map<HADRON_RESCATTERING::HeraParameters::code, std::vector<double> > m_params;

    void InitParams();
  public:
    double xs_tot(HADRON_RESCATTERING::HeraParameters::code tag,const double & s) ;

    HeraParam();
    ~HeraParam();
};

}

#endif