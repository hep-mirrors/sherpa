#ifndef Cut_Data_h
#define Cut_Data_h

#include "Flavour.H"
#include <map>

namespace ATOOLS {

  class Cut_Data {
    Flavour * fl;
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
    void Init(int,Flavour *);
    //void Init(Cut_Data *,Flavour *);
    void Update(double,double);

    void Complete();
    void Reset(bool);
    const double Smin() const { return smin; }
    double Getscut(std::string);
    void Setscut(std::string,double);
  };    


  std::ostream & operator<<(std::ostream &, Cut_Data &);
}

#endif







