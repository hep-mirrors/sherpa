#ifndef AMEGIC_Phasespace_Channel_Generator_Decays_H
#define AMEGIC_Phasespace_Channel_Generator_Decays_H

#include "AMEGIC++/Phasespace/Channel_Generator_Base.H"


namespace AMEGIC {

typedef std::map<std::string,std::string> Decls;
  
  typedef std::vector<std::string> String_List;

  class Channel_Generator_Decays : public Channel_Generator_Base {
    Decls declarations;
    
    int         extrachannelflag,newchannel;
    int         acount,m_aid;
    String_List m_idc;
    std::string m_idstr;

    void        Step0(int,Point *,int&,std::ofstream&,ATOOLS::Flavour*,int&);
    bool        StepS(int,Point *,int&,std::ofstream&,ATOOLS::Flavour*,int&);

    void        GenerateMasses(int,Point **,int,int&,std::ofstream&);
    std::string LinkedMasses(Point *);
    void        CalcSmin(int,const char*,std::string,std::ofstream&,Point*);
    bool        Massive(ATOOLS::Flavour fl) { return (fl.Mass()!=0.) ? 1:0; }
    std::string Order(std::string);
    std::string IString(int);
    void        AddToVariables(int,const std::string&,const std::string&,const int&,std::ofstream&); 
    void        ClearDeclarations() { declarations.clear(); } 
    double      PMassSum(Point*,std::vector<int>*);
  public:
    Channel_Generator_Decays(int,int,Point *,int);
    ~Channel_Generator_Decays();
    int         MakeChannel(int&,int,std::string&,std::string&);
    Point     * GetPointlist()             { return plist; }
    std::string CreateChannelID(int);
  };
}  
#endif

