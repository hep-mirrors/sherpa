#include "ATOOLS/Org/Yaml_Reader.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_File.H"

#include <cassert>

using namespace ATOOLS;

Yaml_Reader::Yaml_Reader(const std::string& name)
  : m_name{name}
{ }

Yaml_Reader::Yaml_Reader(std::istream& s)
{
  Parse(s);
}

Yaml_Reader::Yaml_Reader(const std::string& path, const std::string& filename)
  : m_name {(path.empty() ? "" : path + "/") + filename}
{
  assert(filename != "");
  My_File<std::ifstream> file {path, filename};
  if (!file.Open()) {
    THROW(invalid_input, filename + " could not be opened.");
  }
  try {
    Parse(*file);
  } catch (const std::exception& e) {
    MyStrStream str;
    str << path << '/' << filename << " appears to contain a syntax ";
    // append yaml-cpp error without the "yaml-cpp: " prefix
    str << std::string{e.what()}.substr(10);
    THROW(fatal_error, str.str());
  }
}

std::string Yaml_Reader::Name() const
{
  return m_name;
}

void Yaml_Reader::Parse(std::istream& s)
{
  m_nodes = SHERPA_YAML::LoadAll(s);
}

bool Yaml_Reader::IsParameterCustomised(const Settings_Keys& keys)
{
  const auto node = NodeForKeys(keys);
  return !node.IsNull();
}

std::vector<Settings_Keys> Yaml_Reader::AllSettingsKeys()
{
  std::vector<Settings_Keys> keys_vec;
  Settings_Keys base_keys;
  for (const auto& node : m_nodes) {
    AddSettingsKeys(keys_vec, base_keys, node);
  }
  return keys_vec;
}

void Yaml_Reader::AddSettingsKeys(
  std::vector<Settings_Keys>& keys_vec,
  Settings_Keys& current_keys,
  const SHERPA_YAML::Node& node)
{
  switch (node.Type()) {
    case SHERPA_YAML::NodeType::Null:
      break;
    case SHERPA_YAML::NodeType::Scalar:
      if (keys_vec.size() == 0 || keys_vec.back() != current_keys)
        keys_vec.push_back(current_keys);
      break;
    case SHERPA_YAML::NodeType::Sequence:
      for (auto it = node.begin(); it != node.end(); ++it) {
        auto element = *it;
        AddSettingsKeys(keys_vec, current_keys, element);
      }
      break;
    case SHERPA_YAML::NodeType::Map:
      for (auto it = node.begin(); it != node.end(); ++it) {
        auto key = it->first;
        auto value = it->second;
        current_keys.push_back(key.as<std::string>());
        if (keys_vec.size() == 0 || keys_vec.back() != current_keys)
          AddSettingsKeys(keys_vec, current_keys, value);
        current_keys.pop_back();
      }
      break;
    case SHERPA_YAML::NodeType::Undefined:
      break;
  }
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

bool Yaml_Reader::IsScalar(const Settings_Keys& scopekeys)
{
  const auto node = NodeForKeys(scopekeys);
  return node.IsScalar();
}

bool Yaml_Reader::IsList(const Settings_Keys& scopekeys)
{
  const auto node = NodeForKeys(scopekeys);
  return node.IsSequence();
}

bool Yaml_Reader::IsMap(const Settings_Keys& scopekeys)
{
  const auto node = NodeForKeys(scopekeys);
  return node.IsMap();
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

std::vector<std::string> Yaml_Reader::GetFlattenedStringVectorWithDelimiters(
    const Settings_Keys& keys,
    const std::string& open_delimiter,
    const std::string& close_delimiter)
{
  std::vector<std::string> values;
  const auto node = NodeForKeys(keys);
  if (node.IsNull())
    return values;

  // auto-wrap scalars in a vector
  if (node.IsScalar()) {
    values.push_back(node.as<std::string>());
  } else if (node.IsSequence()) {
    const auto items_count = GetItemsCount(keys);
    for (int i {0}; i < items_count; ++i) {
      auto subkeys = keys;
      subkeys.emplace_back(i);
      const auto newvalues =
        GetFlattenedStringVectorWithDelimiters(subkeys,
                                               open_delimiter,
                                               close_delimiter);
      values.push_back(open_delimiter);
      values.insert(values.end(), newvalues.begin(), newvalues.end());
      values.push_back(close_delimiter);
    }
  }

  return values;
}

SHERPA_YAML::Node Yaml_Reader::NodeForKeys(const Settings_Keys& keys)
{
  SHERPA_YAML::Node node;
  for (const auto& subnode : m_nodes) {
    node.reset(NodeForKeysInNode(keys, subnode));
    if (!node.IsNull()) {
      return node;
    }
  }
  return SHERPA_YAML::Node{SHERPA_YAML::NodeType::Null};
}

SHERPA_YAML::Node Yaml_Reader::NodeForKeysInNode(const Settings_Keys& keys,
                                                 const SHERPA_YAML::Node& node)
{
  static const SHERPA_YAML::Node NullNode{ SHERPA_YAML::NodeType::Null };
  if (!node)
    return NullNode;
  // we can not use assigment, instead we use reset(),
  // cf. https://github.com/jbeder/yaml-cpp/issues/208
  SHERPA_YAML::Node currentnode;
  currentnode.reset(node);
  for (const auto& key : keys) {
    if (key.IsIndex()) {
      if (currentnode.IsSequence()) {
        if (key.GetIndex() < currentnode.size()) {
          currentnode.reset(currentnode[key.GetIndex()]);
        } else {
          return NullNode;
        }
      } else if (key.GetIndex() != 0) {
        return NullNode;
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
