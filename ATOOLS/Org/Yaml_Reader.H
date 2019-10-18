#ifndef ATOOLS_Org_Yaml_Reader_H
#define ATOOLS_Org_Yaml_Reader_H

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/YAML/yaml-cpp/yaml.h"
#include <iostream>

namespace ATOOLS {

  class Settings_Keys;

  class Yaml_Reader {

  public:

    // constructors
    Yaml_Reader();
    Yaml_Reader(std::istream&);
    Yaml_Reader(const std::string& filename);

    std::vector<std::string> GetKeys(const Settings_Keys& scopekeys);
    bool IsList(const Settings_Keys& scopekeys);
    size_t GetItemsCount(const Settings_Keys& scopekeys);

    bool IsParameterCustomised(const Settings_Keys& scopekeys);

    template <typename T>
    T GetScalar(const Settings_Keys& keys)
    {
      T value;
      const auto node = NodeForKeys(keys);
      if (!node.IsNull())
        value = node.as<T>();
      return value;
    }

    template <typename T>
    std::vector<T> GetVector(const Settings_Keys& keys)
    {
      std::vector<T> values;
      const auto node = NodeForKeys(keys);
      if (node.IsNull())
        return values;

      // auto-wrap scalars in a vector
      if (node.Type() == SHERPA_YAML::NodeType::Scalar) {
        values.push_back(node.as<T>());
      } else {
        values = node.as<std::vector<T>>();
      }

      return values;
    }

    template <typename T>
    std::vector<std::vector<T> > GetMatrix(const Settings_Keys& keys)
    {
      std::vector<std::vector<T> > values;
      const auto node = NodeForKeys(keys);
      if (node.IsNull())
        return values;

      if (node.Type() == SHERPA_YAML::NodeType::Scalar) {
        // auto-wrap scalar
        const auto value = node.as<T>();
        values.push_back(std::vector<T>{value});
      } else if (node.Type() == SHERPA_YAML::NodeType::Sequence) {
        auto ismatrix = true;
        auto isvector = true;
        for (const auto& subnode : node) {
          if (!subnode.IsSequence())
            ismatrix = false;
          if (subnode.Type() != SHERPA_YAML::NodeType::Scalar)
            isvector = false;
          if (!isvector && !ismatrix)
            break;
        }
        if (ismatrix) {
          for (const auto& subnode : node)
            values.push_back(subnode.as<std::vector<T>>());
        } else if (isvector) {
          // auto-wrap vector
          auto valuesvec = node.as<std::vector<T>>();
          values.push_back(valuesvec);
        }
      }
      return values;
    }

  protected:

    void Parse(std::istream&);

  private:

    SHERPA_YAML::Node m_node;

    SHERPA_YAML::Node NodeForKeys(const Settings_Keys& scopekeys);

  };

}

#endif
