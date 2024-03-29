#ifndef MODEL__UFO__UFO_Param_Reader_H
#define MODEL__UFO__UFO_Param_Reader_H

#include <string>
#include <vector>

namespace UFO{

  class UFO_Param_Reader{

  public:
    UFO_Param_Reader(const std::string& filepath);
    ~UFO_Param_Reader() {};

    template<class Read_Type> Read_Type GetEntry(const std::string& block,
                                                 const unsigned int& n,
                                                 const unsigned int& m,
                                                 const Read_Type& def,
                                                 const bool& err=true);
    template<class Read_Type> Read_Type GetEntry(const std::string& block,
                                                 const unsigned int& n,
                                                 const Read_Type& def,
                                                 const bool& err=true);

  private:
    std::vector< std::vector<std::string> > m_lines;
    int m_use_runcard; // use run card for ufo input instead of separate param_card.dat
    template<class Read_Type> Read_Type GetWidth(const unsigned int& n,
                                                 const Read_Type& def,
                                                 const bool& err=true);
    template<class Read_Type> static Read_Type NotFound(const std::string& block,
                                                        const unsigned int& n,
                                                        const unsigned int& m,
                                                        const Read_Type& def,
                                                        const bool& err);
    template<class Read_Type> static Read_Type NotFound(const std::string& block,
                                                        const unsigned int& n,
                                                        const Read_Type& def,
                                                        const bool& err);
    static bool IgnoreCaseCompare(const std::string& a, const std::string& b);
    std::vector< std::vector<std::string> >::const_iterator FindBlock(const std::string& block);
  };
}

#endif
