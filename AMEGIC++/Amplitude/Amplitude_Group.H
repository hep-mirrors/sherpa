#ifndef AMEGIC_Amplitude_Amplitude_Group_H
#define AMEGIC_Amplitude_Amplitude_Group_H

#include "AMEGIC++/String/String_Handler.H"
#include "AMEGIC++/Amplitude/Amplitude_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

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
    ~Amplitude_Group();
    bool IsGroup() { return true; }
    Complex Zvalue(String_Handler*,int);
    Complex Zvalue(int);
    Complex Zvalue(int,int*);
    void    PrintGraph();
    void FillCoupling(String_Handler*); 
    void ClearCalcList(); 
    void KillZList(); 
    void SetStringOn();   
    void SetStringOff();  
    void SetNumber(int& n); 
    void Add(Amplitude_Base*, int sign =1);
    int  Size();
    Amplitude_Base* operator[](const int i) { return (i<(int)graphs.size()) ? graphs[i] : 0; }

    Amplitude_Base* GetAmplitude(const int n);
    void BuildGlobalString(int*,int,Basic_Sfuncs*,ATOOLS::Flavour*,String_Handler*);
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
    void BuildGlobalString(int*,int,Basic_Sfuncs*,ATOOLS::Flavour*,String_Handler*)

  */

}
#endif












