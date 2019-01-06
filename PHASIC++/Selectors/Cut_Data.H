#ifndef PHASIC_Selectors_Cut_Data_H
#define PHASIC_Selectors_Cut_Data_H

#include "ATOOLS/Phys/Flavour.H"
#include <map>

namespace PHASIC {

  class Cut_Data {
    const ATOOLS::Flavour * fl;
    double ** cosmin_save;
    double ** cosmax_save;
    double ** scut_save;
    double *  energymin_save;
    std::map<size_t,double> m_smin_map;
    char GetIndexID(int id);
    double Getscut(std::vector<int>,std::vector<int>,int,int,int);
    void Reset(bool);
  public:
    double ** cosmin;
    double ** cosmax;
    double ** scut;
    double *  energymin;
    double *  etmin;
    int       nin, ncut;
    double    smin, minsijfac;

    Cut_Data();
    ~Cut_Data();
    void Init(int,const ATOOLS::Flavour_Vector &fl);
    void Update(double,double);

    void Complete();
    double Smin() const { return smin; }
    double GetscutAmegic(std::string);
    double Getscut(size_t);
    void Setscut(size_t,double);
  };    


  std::ostream & operator<<(std::ostream &, Cut_Data &);
}

#endif







