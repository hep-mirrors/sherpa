#ifndef MODEL__UFO__UFO_Param_Reader_H
#define MODEL__UFO__UFO_Param_Reader_H

#include <string>

namespace ATOOLS
{
  class Data_Reader;
}

namespace UFO{

  class UFO_Param_Reader{

  public:
    UFO_Param_Reader(const std::string& filename, const std::string& path=std::string("./"));
    ~UFO_Param_Reader();

    template<class Read_Type> Read_Type GetEntry(const std::string& block, const unsigned int& n, const unsigned int& m);
    template<class Read_Type> Read_Type GetEntry(const std::string& block, const unsigned int& n);

  private:
    ATOOLS::Data_Reader* p_dataread;
    template<class Read_Type> Read_Type GetWidth(const unsigned int& n);
    void NotFound(const std::string& block, const unsigned int& n, const unsigned int& m);
    void NotFound(const std::string& block, const unsigned int& n);
    
  };
}

#endif
