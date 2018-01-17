#ifndef PHASIC__EWSudakov__Comix_Interface_H
#define PHASIC__EWSudakov__Comix_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "COMIX/Main/Single_Process.H"

namespace PHASIC {

  class Comix_Interface {
  private:
    typedef std::vector<std::vector<double> > Matrix;
    typedef std::vector<std::vector<Complex> > Complex_Matrix;

    Process_Base *p_proc;
    ATOOLS::Cluster_Amplitude *p_ampl;
    NLOTypeStringProcessMap_Map m_apmap;
    
    void InitializeProcesses();
    void Construct_Ampli(Complex ampl1, Complex ampl2);
    Process_Base * GetComixProc(Cluster_Amplitude * const ampl);
  
  public:
    Comix_Interface(Process_Base *const proc, ATOOLS::Cluster_Amplitude *ampl);
    ~Comix_Interface();

    Complex amplitude(ATOOLS::Cluster_Amplitude *const ampl);
    double  Differential(ATOOLS::Cluster_Amplitude *const ampl);
    double  Born();
  };
} 
#endif
