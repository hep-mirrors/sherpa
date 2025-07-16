#ifndef AMEGIC_Amplitude_Amplitude_Output_H
#define AMEGIC_Amplitude_Amplitude_Output_H

#include <fstream>
#include <string> 
#include "AMEGIC++/Main/Point.H"
#include "AMEGIC++/Main/Topology.H"

namespace AMEGIC {

  class Amplitude_Output {
    int counter;
    int maincounter;
    int subcounter;
    bool super_amplitude;
    std::ofstream  pios;
    int m_graphmode;
    int nin, nout, nmed;   // total
    int ic, oc, mc;        // counters
    std::vector<std::string> ins;    // incomming points 
    std::vector<std::string> outs;   // outgoin points
    std::vector<std::string> meds;   // intermediate points

    std::vector<std::string> captions;

    Topology * top;
    Point    * ampl;

    std::string Int2String(const int);
    void LegCount(Point * mo);
    void DrawLine(std::string from, Point * d, bool flip);
    int InclInComming(Point * mo);
    void WriteHeader(const std::string &name);
  public:
    Amplitude_Output(std::string pid, Topology * _top,
                     std::string gpath, int gmode=0);
    void WriteOut(Point *);
    ~Amplitude_Output();
    void BeginSuperAmplitude();
    void EndSuperAmplitude();
  };

  inline void Amplitude_Output::BeginSuperAmplitude() { super_amplitude=true; subcounter=0;}
  inline void Amplitude_Output::EndSuperAmplitude() { super_amplitude=false; if (m_graphmode==1) ++maincounter;}
}

#endif
