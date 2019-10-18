#include "ATOOLS/Org/Settings.H"

#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings_Writer.H"
#include "ATOOLS/Org/Data_Writer.H"

using namespace ATOOLS;

std::unique_ptr<Settings> Settings::mainsettings = nullptr;

void Settings::InitializeMainSettings(int argc, char* argv[])
{
  mainsettings = std::unique_ptr<Settings>(new Settings(argc, argv));
}

void Settings::FinalizeMainSettings()
{
  Settings_Writer().WriteSettings(*mainsettings);
}

Settings& Settings::GetMainSettings()
{
  return *mainsettings;
}

Settings::Settings()
{ }

Settings::Settings(const std::string& yaml)
{
  MyStrStream s{ yaml };
  m_yamlreaders.emplace_back(new Yaml_Reader{s});
}

Settings::Settings(int argc, char* argv[])
{
  m_yamlreaders.emplace_back(new Command_Line_Interface{argc, argv});
  const auto files = GetConfigFiles();
  for (auto it = files.rbegin(); it != files.rend(); ++it)
    m_yamlreaders.emplace_back(new Yaml_Reader{GetPath() + *it});

  Settings_Keys tagkeys{ Setting_Key{"TAGS"} };
  for (auto it = m_yamlreaders.rbegin(); it != m_yamlreaders.rend(); ++it) {
    const auto tags = (*it)->GetKeys(tagkeys);
    for (const auto& tag : tags) {
      tagkeys.push_back(Setting_Key{tag});
      AddTag(tag, (*it)->GetScalar<std::string>(tagkeys));
      tagkeys.pop_back();
    }
  }
}

Scoped_Settings Settings::operator[](const std::string& scope)
{
  return Scoped_Settings{ *this, scope };
}

std::vector<std::string> Settings::GetKeys()
{
  return GetKeys(Settings_Keys{});
}

void Settings::DeclareVectorSettingsWithEmptyDefault(
    const std::vector<std::string>& keys)
{
  DeclareVectorSettingsWithEmptyDefault(keys, Settings_Keys{});
}

void Settings::DeclareVectorSettingsWithEmptyDefault(
    const std::vector<std::string>& keys,
    const Settings_Keys& scopekeys)
{
  for (const auto& key : keys) {
    std::vector<std::string> fullkey{ scopekeys.IndizesRemoved() };
    fullkey.push_back(key);
    SetDefaultMatrix(fullkey, String_Matrix{{}});
  }
}

void Settings::DeclareMatrixSettingsWithEmptyDefault(
    const std::vector<std::string>& keys)
{
  DeclareMatrixSettingsWithEmptyDefault(keys, Settings_Keys{});
}

void Settings::DeclareMatrixSettingsWithEmptyDefault(
    const std::vector<std::string>& keys,
    const Settings_Keys& scopekeys)
{
  for (const auto& key : keys) {
    std::vector<std::string> fullkey{ scopekeys.IndizesRemoved() };
    fullkey.push_back(key);
    SetDefaultMatrix(fullkey, String_Matrix{});
  }
}

bool Settings::IsCustomised(const Settings_Keys& keys)
{
  for (auto& reader : m_yamlreaders)
    if (reader->IsParameterCustomised(keys))
      return true;
  return false;
}

std::vector<std::string> Settings::GetKeys(const Settings_Keys& scopekeys)
{
  std::vector<std::string> keys;
  for (auto& reader : m_yamlreaders) {
    std::vector<std::string> yamlkeys{ reader->GetKeys(scopekeys) };
    keys.insert(keys.end(), yamlkeys.begin(), yamlkeys.end());
  }
  keys.erase(AllUnique(keys.begin(), keys.end()), keys.end());
  return keys;
}

std::string Settings::GetPath()
{
  std::vector<std::string> possible_keys{ "STATUS_PATH", "PATH" };
  for (auto& reader : m_yamlreaders) {
    for (const auto& key : possible_keys) {
      std::string path = reader->GetScalar<std::string>(Settings_Keys{ key });
      if (path != "") {
        if (path[path.size() - 1] != '/') path += "/";
        return path;
      }
    }
  }
  return "";
}

String_Vector Settings::GetConfigFiles()
{
  for (auto& reader : m_yamlreaders) {
    const auto filenames
      = reader->GetVector<std::string>(Settings_Keys{"RUNDATA"});
    if (!filenames.empty())
      return filenames;
  }
  return {"Sherpa.yaml"};
}

bool Settings::IsList(const Settings_Keys& keys)
{
  for (auto& reader : m_yamlreaders) {
    if (reader->IsList(keys))
      return true;
  }
  return false;
}

size_t Settings::GetItemsCount(const Settings_Keys& keys)
{
  for (auto& reader : m_yamlreaders) {
    size_t count{ reader->GetItemsCount(keys) };
    if (count != 0)
      return count;
  }
  return 0;
}

void Settings::AddGlobalTag(const std::string& key, const std::string& value)
{
  if (m_globaltags.find(key) != m_globaltags.end()) {
    THROW(fatal_error, "The global `" + key + "' tag is already set.");
  }
  m_globaltags[key] = value;
}

void Settings::AddTag(const std::string& key, const std::string& value)
{
  m_tags[key] = value;
}

std::string Settings::GetScalarDefault(const Defaults_Key& k,
                                       const Defaults& d)
{
  String_Vector defaultvector{ GetVectorDefault(k, d) };
  if (defaultvector.empty())
    return "";
  if (defaultvector.size() != 1)
    THROW(fatal_error, "The default for " + k.back()
        + " does not have right dimensions.");
  return defaultvector[0];
}

String_Vector Settings::GetVectorDefault(const Defaults_Key& k,
                                         const Defaults& d)
{
  String_Matrix defaultmatrix{ GetMatrixDefault(k, d) };
  if (defaultmatrix.empty())
    return String_Vector{};
  if (defaultmatrix.size() != 1)
    THROW(fatal_error, "The default for " + k.back()
        + " does not have right dimensions.");
  return defaultmatrix[0];
}

String_Matrix Settings::GetMatrixDefault(const Defaults_Key& k,
                                         const Defaults& d)
{
  const auto it = d.find(k);
  if (it == d.end())
    THROW(fatal_error, "The default for " + k.back() + " has not been set.");
  return it->second;
}

void Settings::ReplaceTags(std::string& value)
{
  if (value.empty()) {
    return;
  }
  std::string original_value(value);
  // iterate through occurrences of `$(<tag>)'
  size_t begin_pos(0);
  while ((begin_pos = value.find("$(", begin_pos)) != std::string::npos) {
    size_t end_pos(value.find(")", begin_pos));
    if (end_pos == std::string::npos) {
      THROW(fatal_error, std::string("An unbalanced `$(' has been found in `")
                         + original_value
                         + "'. This is an input error, please check your"
                         + " command line input and your run card.");
    }
    std::string tag(value.substr(begin_pos + 2, end_pos - begin_pos - 2));

    // iterate through tag definitions to replace `$(<tag>)' with the
    // corresponding value
    bool found(false);
    for (auto const& defaulttag : m_tags) {
      if (defaulttag.first == tag) {
        value.replace(begin_pos, end_pos - begin_pos + 1, defaulttag.second);
        found = true;
        break;
      }
    }
    if (found == false) {
      THROW(fatal_error, std::string("Could not find a definition for `")
                         + value.substr(begin_pos, end_pos - begin_pos + 1)
                         + "'. This is an input error, please check your"
                         + " command line input and your run card.");
    }
  }

  // replace global tags
  for (auto const& tag : m_globaltags) {
    value = StringReplace(value, tag.first, tag.second);
  }
}

void Settings::SetDefault(const Settings_Keys& keys, const char* value)
{
  SetDefault<std::string>(keys, std::string(value));
}

std::string Settings::ApplyReplacements(const Settings_Keys& settings_keys,
                                        const std::string& value)
{
  const std::vector<std::string> keys{ settings_keys.IndizesRemoved() };
  const auto it = m_replacements.find(keys);
  if (it == m_replacements.end())
    return value;
  for (const auto& replacement : it->second) {
    if (value == replacement.first) {
      return replacement.second;
    }
  }
  return value;
}

void Settings::SetDefaultSynonyms(const Settings_Keys& settings_keys,
                                  const std::vector<std::string>& synonyms)
{
  const Defaults_Key keys{ settings_keys.IndizesRemoved() };
  const auto it = m_defaultsynonyms.find(keys);
  if (m_defaultsynonyms.find(keys) != m_defaultsynonyms.end())
    if (synonyms != it->second)
      THROW(fatal_error, "A different default synonyms list for "
          + keys.back() + " has already been set.");
  m_defaultsynonyms[keys] = synonyms;
}

bool Settings::IsDefaultSynonym(const Settings_Keys& settings_keys,
                                const std::string& value)
{
  const std::vector<std::string> keys{ settings_keys.IndizesRemoved() };
  const auto it = m_defaultsynonyms.find(keys);
  if (it == m_defaultsynonyms.end())
    return false;
  const auto& v = it->second;
  return (std::find(v.begin(), v.end(), value) != v.end());
}
