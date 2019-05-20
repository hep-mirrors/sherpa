#include "ATOOLS/Org/Yaml_Reader.H"

#include "ATOOLS/Org/Settings_Keys.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <cassert>

using namespace ATOOLS;

Yaml_Reader::Yaml_Reader()
{ }

Yaml_Reader::Yaml_Reader(std::istream& s)
{
  Parse(s);
}

Yaml_Reader::Yaml_Reader(const std::string& filename)
{
  assert(filename != "");
  std::ifstream file(filename);
  try {
    Parse(file);
  } catch (const std::exception& e) {
    MyStrStream str;
    str << filename << " appears to contain a syntax ";
    // append yaml-cpp error wihtout the "yaml-cpp: " prefix
    str << std::string{e.what()}.substr(10);
    THROW(fatal_error, str.str());
  }
}

void Yaml_Reader::Parse(std::istream& s)
{
  m_node = SHERPA_YAML::Load(s);
}

bool Yaml_Reader::IsParameterCustomised(const Settings_Keys& keys)
{
  const auto node = NodeForKeys(keys);
  return !node.IsNull();
}

std::vector<std::string> Yaml_Reader::GetKeys(const Settings_Keys& scopekeys)
{
  std::vector<std::string> keys;
  const auto node = NodeForKeys(scopekeys);
  if (node.IsNull())
    return keys;
  assert(node.IsMap());
  for (const auto& subnode : node) {
    keys.push_back(subnode.first.as<std::string>());
  }
  return keys;
}

size_t Yaml_Reader::GetItemsCount(const Settings_Keys& scopekeys)
{
  const auto node = NodeForKeys(scopekeys);
  if (node.IsNull())
    return 0;
  else if (node.IsSequence())
    return node.size();
  else if (node.IsMap())
    return 0;
  else
    return 1;
}

SHERPA_YAML::Node Yaml_Reader::NodeForKeys(const Settings_Keys& keys)
{
  static const SHERPA_YAML::Node NullNode{ SHERPA_YAML::NodeType::Null };
  if (!m_node)
    return NullNode;
  // we can not use assigment, instead we use reset(),
  // cf. https://github.com/jbeder/yaml-cpp/issues/208
  SHERPA_YAML::Node currentnode;
  currentnode.reset(m_node);
  for (const auto& key : keys) {
    if (key.IsIndex()) {
      if (currentnode.IsSequence()) {
        if (key.GetIndex() < currentnode.size()) {
          currentnode.reset(currentnode[key.GetIndex()]);
        } else {
          THROW(fatal_error, "There is no entry at position "
                             + ToString(key.GetIndex()) + ".");
        }
      } else if (key.GetIndex() != 0) {
        THROW(fatal_error, "The current node has no entries.");
      }
      // Note that we ignore indizes that are zero in the case of non-sequences,
      // i.e. a zero index is an identity operator for these cases, leaving
      // currentnode untouched
    } else {
      if (!currentnode.IsMap())
        return NullNode;
      const auto child = currentnode[key.GetName()];
      if (child)
        currentnode.reset(child);
      else
        return NullNode;
    }
  }
  return currentnode;
}
