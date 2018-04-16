#ifndef AMISIC_Tools_MI_Parameters_H
#define AMISIC_Tools_MI_Parameters_H

#include "ATOOLS/Org/Default_Reader.H"

#include <map>
#include <string>

namespace AMISIC {
  struct overlap_form {
    enum code {
      Single_Gaussian = 1,
      Double_Gaussian = 2,
      unknown         = -1
    };
  };
  
  class MI_Parameters {
  private:
    overlap_form::code m_overlapform;
    std::map<std::string,double> m_parameters;

    void   ReadParameters(ATOOLS::Default_Reader *const defaultreader);
    double CalculatePT(const double & pt);
  public:
    MI_Parameters();
    ~MI_Parameters() {}

    bool Init(ATOOLS::Default_Reader *const defaultreader);
    double operator()(std::string keyword);

    const overlap_form::code & GetOverlapForm() { return m_overlapform; } 
  };
  
  extern MI_Parameters * mipars;
}
#endif
