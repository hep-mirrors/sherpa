#ifndef Data_Collector_H
#define Data_Collector_H

#include <string>
#include "Blob.H"

namespace ATOOLS {  

  struct Process_Info {
    std::string texname;
    double  xsec, xsec_err;
    Process_Info(const std::string & tname, double xs, double xserr)
      : texname(tname), xsec(xs), xsec_err(xserr)
    {}
    friend std::ostream & operator<<(std::ostream & s, const Process_Info & wi);
  };

  class Data_Collector {
    static int s_counter;
    static ATOOLS::String_BlobDataBase_Map s_datacontainer;
  public:
    Data_Collector();
    static void AddData(const std::string name, Blob_Data_Base * data);
    ATOOLS::Blob_Data_Base * operator[](const std::string name);
    static void Print();
    ~Data_Collector();
  };


}

#endif
