#ifndef Amplitude_Group_h
#define Amplitude_Group_h

#include "Point.H"
#include "String_Handler.H"
#include "Amplitude_Base.H"
#include "Message.H"

#include <map>

namespace AMEGIC {

  typedef std::vector<Amplitude_Base *>          Amplitude_List;
  typedef std::map<std::string, Amplitude_List>  Graph_Table;

  struct Graph_Family {
    int             was_clustered;
    int             znumber;
    int             topnumber;
    int             permnumber;
    Amplitude_List  graphs;
    std::string     banner;
  };

  typedef std::vector<Graph_Family *>    Graph_Families;
  
  class Amplitude_Group: public virtual Amplitude_Base {
    Amplitude_Base* GetSingleGraph(std::list<sknot*>&);
    Graph_Family * FindFamily(int zn, int tn, int pn);
  protected:
    std::vector<Amplitude_Base*> graphs;
    Graph_Table             graph_table;
    Graph_Families          family_table;
    std::string             groupname;          
  public:
    Complex Zvalue(String_Handler * sh,int ihel) {
      Complex M(0.,0.);
      for (int i=0;i<graphs.size();i++) M += graphs[i]->Zvalue(sh, ihel);
      return M;
    }
    Complex Zvalue(int ihel,int* signlist) {
      Complex M(0.,0.);
      for (int i=0;i<graphs.size();i++) M += graphs[i]->Zvalue(ihel,signlist);
      return M;
    }

    void    PrintGraph() {
      AORGTOOLS::msg.Tracking()<<"Group: "<<groupname<<std::endl;
      for (int i=0;i<graphs.size();i++) graphs[i]->PrintGraph();
    }

    void FillCoupling(String_Handler* shand) 
       {for (int i=0;i<graphs.size();i++) graphs[i]->FillCoupling(shand);}
    void ClearCalcList() {for (int i=0;i<graphs.size();i++) graphs[i]->ClearCalcList();}

    void SetStringOn()   {buildstring=1;for (int i=0;i<graphs.size();i++) graphs[i]->SetStringOn();}
    void SetStringOff()  {buildstring=0;for (int i=0;i<graphs.size();i++) graphs[i]->SetStringOff();}
    void SetNumber(int& n) {for (int i=0;i<graphs.size();i++) graphs[i]->SetNumber(n);}
    void Add(Amplitude_Base* ab, int sign =1) {
      if (sign==-1)
	ab->SetSign(sign*ab->GetSign());
      graphs.push_back(ab);
    }
    int  Size()                  {return graphs.size();}
    Amplitude_Base* operator[](const int i) {return (i<graphs.size()) ? graphs[i] : 0;}

    Amplitude_Base* GetAmplitude(const int n) {
      for (int i=0;i<graphs.size();i++) {
	Amplitude_Base* help = graphs[i]->GetAmplitude(n);
	if (help!=0) return help;
      }
      return 0;
    }

    void BuildGlobalString(int*,int,Basic_Sfuncs*,APHYTOOLS::Flavour*,String_Handler*);
  };

  /*! 
    \file
    \brief this file contains the class AMEGIC::Amplitude_Group
  
  */

  /*! \class Amplitude_Group
      \brief is the basis class for grouping Feynman amplitudes in
             Color_Group's or Super_Amplitudes
  */

  /*!
    \fn Amplitude_Base* GetSingleGraph(std::list<sknot*>&)
    \brief
  */

  /*!
    \fn Complex Zvalue(String_Handler * sh,int ihel) 
    \brief
  */

  /*!
    \fn Complex Zvalue(int ihel,int* signlist) 

    void    PrintGraph()
    void FillCoupling(String_Handler* shand) 
    void ClearCalcList()

    void SetStringOn()  
    void SetStringOff() 
    void SetNumber(int& n) 
    void Add(Amplitude_Base* ab, int sign =1)
    int  Size() 
    Amplitude_Base* operator[](const int i) 

    Amplitude_Base* GetAmplitude(const int n)
    void BuildGlobalString(int*,int,Basic_Sfuncs*,APHYTOOLS::Flavour*,String_Handler*)

  */

}
#endif












