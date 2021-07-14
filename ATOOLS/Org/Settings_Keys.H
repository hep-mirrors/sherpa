#ifndef ATOOLS_Org_SettingsKeys_H
#define ATOOLS_Org_SettingsKeys_H

#include <string>
#include <vector>
#include <limits>

namespace ATOOLS {

  class Setting_Key {

  public:

    Setting_Key(size_t _index) :
      name{ "" },
      index{ _index }
    {};

    Setting_Key(const char* _name) : Setting_Key{ std::string{_name} } {}

    Setting_Key(const std::string& _name) :
      name{ _name },
      index{ std::numeric_limits<size_t>::max() }
    {};

    bool IsIndex() const;
    std::string GetName() const;
    size_t GetIndex() const;

    bool operator==(const Setting_Key&) const;

  private:

    const std::string name;
    const size_t index;

  };

  class Settings_Keys : public std::vector<Setting_Key> {

  public:

    Settings_Keys() {}

    using std::vector<Setting_Key>::vector;

    Settings_Keys(const std::vector<std::string>&);

    std::string Name() const;

    std::vector<std::string> IndizesRemoved() const;
    bool ContainsNoIndizes() const;

  };

}

#endif
