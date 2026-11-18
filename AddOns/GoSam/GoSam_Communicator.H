#ifndef GoSam_Communicator_H
#define GoSam_Communicator_H

#include <vector>
#include <string>
#include "ATOOLS/Phys/Flavour.H"

namespace GoSam {

  class GoSam_Communicator {
    int m_filestatus;
    std::string m_name;

    bool GetPLine(std::ifstream&,std::string&,std::string&);
    static std::string StripProcessName(const std::string &s);
  public:
    GoSam_Communicator(std::string name);
    ~GoSam_Communicator();

    int  FileStatus() { return m_filestatus; }
    void AddParameter(std::string);
    int  CheckParameterStatus();
    int  CheckProcess(int,int,const ATOOLS::Flavour_Vector&);
    void AddProcess(int,int,const ATOOLS::Flavour_Vector&);
    int  GetID(int,int,const ATOOLS::Flavour_Vector&,int);
    std::string& GetName() { return m_name; }

  }; // end of class GoSam_Communicator

}// end of namespace GoSam

#endif
