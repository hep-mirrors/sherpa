#ifndef AMEGIC_Phasespace_Phase_Space_Generator_H
#define AMEGIC_Phasespace_Phase_Space_Generator_H

//#include "PHASIC++/Channels/Multi_Channel.H"
#include "AMEGIC++/String/MyString.H"
#include "AMEGIC++/Main/Point.H"
#include <list>

namespace AMEGIC {
  class Process_Base;

  class Phase_Space_Generator {
    static const int buffersize = 512*32;
    static const int maxchannels = 2048;
    int         nin,nout;
    char        filename[12], procname[12];
    std::string path,pathID,pID;
    int m_mode;

    //bool  LoadChannels(ATOOLS::Flavour *,PHASIC::Multi_Channel *,Process_Base *);
    bool  GetLibList(std::list<std::string>* liblist);
    bool  IsFile(std::string &);
    bool  Search(std::ifstream &,std::string);
    bool  RSearch(std::ifstream &,std::string &);
    int   Search(std::string file,std::string search);
    int   RSearch(std::string file,std::string &search);
    void  Copy(std::string,std::string);
    void  AddToMakefile(std::string Makefile,std::string pathID,std::string fileID);
    void  AddToMakefileAM(std::string, std::string, std::string);
  public:
    Phase_Space_Generator(int,int);
    bool Construct(std::list<std::string>*,std::string,std::string,ATOOLS::Flavour *, 
		   Process_Base *);
    //PHASIC::Single_Channel * SetChannel(int,int,ATOOLS::Flavour* fl,std::string&,ATOOLS::Integration_Info * const);
  };

}
#endif


