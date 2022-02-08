#ifndef ATOOLS_Org_Yaml_Reader_H
#define ATOOLS_Org_Yaml_Reader_H

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Settings_Keys.H"
#include "ATOOLS/YAML/yaml-cpp/yaml.h"
#include <iostream>

namespace ATOOLS {

  class Settings_Keys;

  class Yaml_Reader {

  public:

    // constructors
    Yaml_Reader(const std::string& name="");
    Yaml_Reader(std::istream&);
    Yaml_Reader(const std::string& path, const std::string& filename);

    std::string Name() const;

    std::vector<Settings_Keys> AllSettingsKeys();

    std::vector<std::string> GetKeys(const Settings_Keys& scopekeys);
    bool IsScalar(const Settings_Keys& scopekeys);
    bool IsList(const Settings_Keys& scopekeys);
    bool IsMap(const Settings_Keys& scopekeys);
    size_t GetItemsCount(const Settings_Keys& scopekeys);

    bool IsParameterCustomised(const Settings_Keys& scopekeys);

    std::vector<std::string> GetFlattenedStringVectorWithDelimiters(
        const Settings_Keys&,
        const std::string& open_delimiter="{{",
        const std::string& close_delimiter="}}");

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
          if (subnode.IsSequence()) {
            for (const auto& subsubnode : subnode) {
              if (!subsubnode.IsScalar())
                THROW(invalid_input,
                    "Attempting to read a more than 2-dimensional setting: "
                    + keys.Name());
            }
          } else {
            ismatrix = false;
          }
          if (!subnode.IsScalar())
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

    std::string m_name;

  private:

    std::vector<SHERPA_YAML::Node> m_nodes;

    SHERPA_YAML::Node NodeForKeys(const Settings_Keys& scopekeys);
    SHERPA_YAML::Node NodeForKeysInNode(const Settings_Keys& scopekeys,
                                        const SHERPA_YAML::Node&);

    void AddSettingsKeys(
        std::vector<Settings_Keys>&,
        Settings_Keys&,
        const SHERPA_YAML::Node&);

  };

}

#endif
