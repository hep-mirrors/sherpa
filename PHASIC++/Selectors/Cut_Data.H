#ifndef PHASIC_Selectors_Cut_Data_H
#define PHASIC_Selectors_Cut_Data_H

#include "ATOOLS/Phys/Flavour.H"
#include <map>

namespace PHASIC {

  class Cut_Data {
    ATOOLS::Flavour * fl;
    double ** cosmax_save;
    double ** scut_save;
    double *  energymin_save;
    std::map<std::string,double> m_smin_map;
    char GetIndexID(int id);
    double Getscut(std::vector<char>,std::vector<char>,int,int,int);
  public:
    double ** cosmin;
    double ** cosmax;
    double ** scut;
    double *  energymin;
    double *  energymax;
    double *  etmin;
    int       ncut;
    double    smin;

    Cut_Data();
    ~Cut_Data();
    void Init(int,ATOOLS::Flavour *);
    //void Init(Cut_Data *,Flavour *);
    void Update(double,double);

    void Complete();
    void Reset(bool);
    double Smin() const { return smin; }
    double Getscut(std::string);
    void Setscut(std::string,double);
  };    


  std::ostream & operator<<(std::ostream &, Cut_Data &);
}

#endif







