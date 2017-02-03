// a sketch of what the new API might look like

#include "yaml-cpp/yaml.h"
#include <iostream>

using namespace SHERPA_YAML;
using namespace std;

#ifdef YAML_NAMESPACE
#define YAML YAML_NAMESPACE
#endif


bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

int main(int argc, char *argv[]) {
  
    // test.yaml
    // - foo
    // - primes: [2, 3, 5, 7, 11]
    //   odds: [1, 3, 5, 7, 9, 11]
    // - [x, y]

    // move-like semantics

    Node config = SHERPA_YAML::LoadFile(argv[1]);

    if (config["EVENTS"]) {
      std::cout << "about to generate " << config["EVENTS"].as<std::string>() << " events\n";
    }

    for (unsigned int i=2;i<argc;i++) {
       std::string tt = argv[i];
       replace(tt, ":", ": ");

       SHERPA_YAML::Node cl = SHERPA_YAML::Load(tt);
       cout << cl << " " << tt << std::endl;
       //SHERPA_YAML::Node::const_iterator sub_it = cl.begin();
       //std::cout << sub_it->second << std::endl;
       //std::cout << name << std::endl;

       //std::cout  << "Found a node of type " >> cl.Type().as<std::str> << std::endl;
       for(SHERPA_YAML::const_iterator it=cl.begin();it != cl.end();++it) {
         std::string key = it->first.as<std::string>();
         if (config[key]) {
            std::cout << "Overwriting " << key << " from " << config[key].as<std::string>() << " to " << it->second.as<std::string>() << std::endl;
            config[key] = it->second;
         }
       }
    }

    std::cout << config << std::endl;

    //std::cout << root[0].as<std::string>();       // "foo"
    //std::cout << str(root[0]);                    // "foo", shorthand?
    //std::cout << root[1]["primes"][3].as<int>();  // "7"
    //std::cout << root[1]["odds"][6].as<int>();    // throws?

    //root[2].push_back(5);
    //root[3] = "Hello, World";
    //root[0].reset();
    //root[0]["key"] = "value";

    //std::cout << root;
    // # not sure about formatting
    // - {key: value}
    // - primes: [2, 3, 5, 7, 11]
    //   odds: [1, 3, 5, 7, 9, 11]
    // - [x, y, 5]
    // - Hello, World

  return 0;
}
