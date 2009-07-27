#ifndef LH_OLE_Communicator_H
#define LH_OLE_Communicator_H

#include <vector>
#include <string>
#include "ATOOLS/Phys/Flavour.H"

namespace OLE {

  class LH_OLE_Communicator {
    int m_filestatus;
    std::string m_name;

  public:
    LH_OLE_Communicator(std::string name);
    ~LH_OLE_Communicator();

    int  FileStatus() { return m_filestatus; }
    void AddParameter(std::string);
    int  CheckParameterStatus();
    int  CheckProcess(int,int,const ATOOLS::Flavour_Vector&);
    void AddProcess(int,int,const ATOOLS::Flavour_Vector&);
    int  GetID(int,int,const ATOOLS::Flavour_Vector&,int);
    std::string& GetName() { return m_name; }

  }; // end of class LH_OLE_Communicator

}// end of namespace OLE

#endif
