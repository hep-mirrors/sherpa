#ifndef ATOOLS_Org_Yaml_Reader_H
#define ATOOLS_Org_Yaml_Reader_H

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/YAML/yaml-cpp/yaml.h"
#include <iostream>

//#ifdef YAML_NAMESPACE
//#define YAML YAML_NAMESPACE
//#endif


namespace ATOOLS {

  class Yaml_Reader {
  protected:
   
    SHERPA_YAML::Node m_node; 
    std::string m_string;

  public:

    // constructors
    Yaml_Reader();

    // destructor
    ~Yaml_Reader();

    void Read_File(char*);

    void Print() {std::cout << m_node <<std::endl;};


    template <typename T >
    T Get(std::string key) {
      if (m_node[key]) return m_node[key].as<T>();
      else {
        msg_Error() << "Requested key '" << key << "' not found\n";
        abort();
      }
    }

    template <typename T >
    void Set(std::string key, T value, bool allowOverRide=false) {
      SHERPA_YAML::Node n = SHERPA_YAML::Load(key + ": " + ToString(value));
      for(SHERPA_YAML::const_iterator it=n.begin();it != n.end();++it) {
        std::string key = it->first.as<std::string>();
        if (m_node[key]) {
          if (allowOverRide) {
            msg_Info() << "Overwriting " << key << " from " << m_node[key].as<std::string>() << " to " << it->second.as<std::string>() << std::endl;
            m_node[key] = it->second;
          }
          else msg_Info() << "Not Overwriting " << key << " from " << m_node[key].as<std::string>() << " to " << it->second.as<std::string>();
        }
        else {
          msg_Info() << "setting " << key << " to " << it->second.as<std::string>();
          m_node[key] = it->second;
        }
      }
    }


  }; // end of class Yaml_Reader

} // end of namespace ATOOLS

#endif
