#ifndef Process_Info_h
#define Process_Info_h

#include "Flavour.H"
#include "Pol_Info.H"
#include <vector>

namespace AMEGIC {
  class Point;

  class Process_Info {
  public:
    ATOOLS::Flavour* p_fl;
    Pol_Info* p_pl;
    std::vector<std::vector<Process_Info*> > m_sublist;
  protected:
    Process_Info* GetSubProcess(int n,int& dn);
    int  GetDPOffset(int&,int&);
    Point* MergePointList(Point** plist,Point* np,int &nd, int nin, int &ep);
   public:
    Process_Info(ATOOLS::Flavour*,Pol_Info*);
    ~Process_Info();
    ATOOLS::Flavour* Flav() {return p_fl;}
    Pol_Info* Pol()         {return p_pl;}
    void AddSubList(int,ATOOLS::Flavour*,Pol_Info*);
    void ResetSubList(int,ATOOLS::Flavour*,Pol_Info*);
    int Nout();
    int TotalNout();
    void GetFlavList(ATOOLS::Flavour*,int n=0);
    void GetPolList(Pol_Info*);
    int  GetTotalFlavList(ATOOLS::Flavour*,int n=-1);
    void GetTotalPolList(Pol_Info*);
    Process_Info* FindDM(char);
    bool CheckCompleteness();
    void Expand();
    int  NProcs();
    int  Nmax(int nin);
    int  Ndecays();
    int  GetDPOffset(int);
    Process_Info* GetSubProcess(int n);
    Process_Info* GetDecay(int &n);
    void Reshuffle();
    std::string GenerateName();
    void MergePointList(Point**,Point*,int);

    void Print();
    void FullPrint();
  };

  class Check_External_Flavours {
  public:
    bool ValidProcess(int _nin,ATOOLS::Flavour * _in,int _nout,ATOOLS::Flavour * _out);
  };
  extern Check_External_Flavours CF;                                                  
}



#endif
