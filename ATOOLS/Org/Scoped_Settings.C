#include "ATOOLS/Org/Scoped_Settings.H"

using namespace ATOOLS;

Scoped_Settings::Scoped_Settings():
  m_rootsettings{ &Settings::GetMainSettings() }
{
}

Scoped_Settings::Scoped_Settings(const std::string& yamlstring):
  m_ownedsettings{ new Settings(yamlstring) },
  m_rootsettings{ m_ownedsettings.get() }
{
}

Scoped_Settings::Scoped_Settings(const Scoped_Settings& s):
  m_ownedsettings{ s.m_ownedsettings },
  m_rootsettings{ s.m_rootsettings },
  m_scopes{ s.m_scopes }
{
}

Scoped_Settings::Scoped_Settings(Settings &rootsettings,
                                 const std::string& scope):
  m_rootsettings{ &rootsettings },
  m_scopes{ Setting_Key{scope} }
{
}

Scoped_Settings Scoped_Settings::operator[](const std::string& scope) const
{
  return Scoped(Setting_Key{scope});
}
Scoped_Settings Scoped_Settings::operator[](size_t scope) const
{
  return Scoped(Setting_Key{scope});
}

Scoped_Settings& Scoped_Settings::operator=(Scoped_Settings other)
{
  std::swap(m_ownedsettings, other.m_ownedsettings);
  if (m_ownedsettings)
    m_rootsettings = m_ownedsettings.get();
  else
    m_rootsettings = other.m_rootsettings;
  std::swap(m_scopes, other.m_scopes);
  return *this;
}

std::vector<std::string> Scoped_Settings::GetKeys()
{
  return m_rootsettings->GetKeys(m_scopes);
}

void Scoped_Settings::DeclareVectorSettingsWithEmptyDefault(
    const std::vector<std::string>& keys)
{
  m_rootsettings->DeclareVectorSettingsWithEmptyDefault(keys, m_scopes);
}

void Scoped_Settings::DeclareMatrixSettingsWithEmptyDefault(
    const std::vector<std::string>& keys)
{
  m_rootsettings->DeclareMatrixSettingsWithEmptyDefault(keys, m_scopes);
}

Scoped_Settings& Scoped_Settings::UseNoneReplacements()
{
  static std::map<std::string, std::string> nonelist{
    {"Off", "None"},
    {"0", "None"},
    {"false", "None"},
    {"no", "None"}
  };
  return SetReplacementList(nonelist);
}

bool Scoped_Settings::IsCustomised()
{
  return m_rootsettings->IsCustomised(m_scopes);
}

bool Scoped_Settings::IsList() const
{
  return m_rootsettings->IsList(m_scopes);
}

std::vector<Scoped_Settings> Scoped_Settings::GetItems() const
{
  std::vector<Scoped_Settings> scoped_settings;
  const size_t count{ m_rootsettings->GetItemsCount(m_scopes) };
  scoped_settings.reserve(count);
  for (size_t i{ 0 }; i < count; ++i) {
    scoped_settings.push_back(Scoped(Setting_Key{i}));
  }
  return scoped_settings;
}

size_t Scoped_Settings::GetItemsCount()
{
  return m_rootsettings->GetItemsCount(m_scopes);
}

Scoped_Settings Scoped_Settings::GetItemAtIndex(const size_t& i) const
{
  return Scoped(Setting_Key{i});
}

size_t Scoped_Settings::GetIndex() const
{
  if (m_scopes.empty() || !m_scopes.back().IsIndex())
    return 0;
  return m_scopes.back().GetIndex();
}

Scoped_Settings Scoped_Settings::Scoped(const Setting_Key& scope) const
{
  Scoped_Settings settings{ *this };
  settings.AppendScope(scope);
  return settings;
}

void Scoped_Settings::AppendScope(const Setting_Key& scope)
{
  m_scopes.push_back(scope);
}
