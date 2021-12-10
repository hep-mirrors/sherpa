#include "ATOOLS/Org/Settings_Keys.H"
#include <algorithm>

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;

bool Setting_Key::IsIndex() const
{
  return (index != std::numeric_limits<size_t>::max());
}

std::string Setting_Key::GetName() const
{
  if (IsIndex()) THROW(fatal_error, "Settings_Key name undefined.");
  return name;
}

size_t Setting_Key::GetIndex() const
{
  if (!IsIndex()) THROW(fatal_error, "Settings_Key index undefined.");
  return index;
}

bool Setting_Key::operator==(const Setting_Key& rhs) const
{
  if (IsIndex() != rhs.IsIndex())
    return false;
  if (IsIndex())
    return index == rhs.index;
  else
    return name == rhs.name;
}

Settings_Keys::Settings_Keys(const std::vector<std::string>& strings)
{
  reserve(strings.size());
  std::transform(strings.begin(), strings.end(),
      std::back_inserter(*this),
      [](std::string s) -> Setting_Key { return Setting_Key{s}; });
}

std::string Settings_Keys::Name() const
{
  MyStrStream s;
  auto keys = IndizesRemoved();
  for (const auto& key : keys)
    s << key << ":";
  auto name = s.str();
  if (!name.empty())
    return name.substr(0, name.size() - 1);
  else
    return name;
}

std::vector<std::string> Settings_Keys::IndizesRemoved() const
{
  std::vector<std::string> filtered_keys;
  filtered_keys.reserve(size());
  for (const Setting_Key& k : *this)
    if (!k.IsIndex())
      filtered_keys.push_back(k.GetName());
  return filtered_keys;
}

bool Settings_Keys::ContainsNoIndizes() const
{
  const_iterator it{ std::find_if(begin(), end(),
      [](const Setting_Key& k) { return k.IsIndex(); }) };
  return (it == end());
}
