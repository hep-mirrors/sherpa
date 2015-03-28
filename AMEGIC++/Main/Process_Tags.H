#ifndef AMEGIC_Main_Process_Tags_H
#define AMEGIC_Main_Process_Tags_H

#include "ATOOLS/Phys/Flavour.H"
#include "AMEGIC++/Main/Pol_Info.H"
#include "AMEGIC++/Main/Tools.H"
#include "PHASIC++/Process/Process_Info.H"
#include <vector>
#include <map>

namespace AMEGIC {
  class Point;
  typedef std::map<int, int> FMMap;

  class Process_Tags {
  public:
    ATOOLS::Flavour* p_fl;
    Pol_Info* p_pl;
    std::vector<std::vector<Process_Tags*> > m_sublist;
    int m_maxqcdjets,m_zerowidth;
  protected:
    Process_Tags* GetSubProcess(int n,int& dn);
    int  GetDPOffset(int&,int&);
    Point* MergePointList(Point** plist,Point* np,int &nd, int nin, int &ep);
   public:
    Process_Tags(ATOOLS::Flavour*,Pol_Info*);
    Process_Tags(Process_Tags*);
    ~Process_Tags();
    ATOOLS::Flavour* Flav() {return p_fl;}
    Pol_Info* Pol()         {return p_pl;}
    void AddSubList(int,ATOOLS::Flavour*,Pol_Info*);
    void ResetSubList(int,ATOOLS::Flavour*,Pol_Info*);
    int Nout();
    int TotalNout();
    void GetFlavList(ATOOLS::Flavour*,int n=0);
    size_t GetStableFlavList(ATOOLS::Flavour* fl, int n=0);
    void GetPolList(Pol_Info*);
    int  GetTotalFlavList(ATOOLS::Flavour*,int n=-1);
    size_t GetOnshellFlavList(ATOOLS::Flavour_Vector &fl, std::vector<Process_Tags*> &decaylist, bool first=1);
    int  OSDecays();
    void GetTotalPolList(Pol_Info*);
    Process_Tags* FindDM(std::string);
    bool CheckCompleteness();
    void Expand();
    int  NProcs();
    int  Nmax(int nin);
    int  Ndecays();
    int  GetDPOffset(int);
    Process_Tags* GetSubProcess(int n);
    Process_Tags* GetDecay(int &n);
    void Reshuffle(Process_Tags* =0);
    std::string GenerateName();
    void MergePointList(Point**,Point*,int);

    void GetOSConditions(std::vector<std::pair<std::string, double> >& osc,int &cnt);
    void Print();
    void FullPrint();

    std::string PNID() const;
  };

  class Check_External_Flavours {
  public:
    bool ValidProcess(int _nin,ATOOLS::Flavour * _in,int _nout,ATOOLS::Flavour * _out);
    bool PureGluonic(int _nin,ATOOLS::Flavour * _in,int _nout,ATOOLS::Flavour * _out);
////////////////////////////////////// MHV //////////////////////////////////////////
    bool MHVCalculable(const PHASIC::Process_Info&);
/////////////////////////////////////////////////////////////////////////////////////
  };
  extern Check_External_Flavours CF;                                                  
}



#endif
