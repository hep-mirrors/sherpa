#ifndef ATOOLS_Org_Settings_H
#define ATOOLS_Org_Settings_H

#include "ATOOLS/Org/Command_Line_Interface.H"
#include "ATOOLS/Org/Settings_Keys.H"
#include "ATOOLS/Org/Yaml_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Read_Write_Base.H"

#include <vector>
#include <unordered_set>

namespace ATOOLS {

  class Settings_Keys;
  class Scoped_Settings;
  class Settings_Writer;

  /**
   * Settings organises the three ways run time parameters can be set.
   *
   * Overrides (`Override...`) have the highest priority, the command line
   * has the second highest priority, then the YAML run card is searched. If
   * the above do not provide the setting, the default set via SetDefault... is
   * returned.
   *
   * In fact, SetDefault... *must* be called before accessing a setting to make
   * sure that there are not suprises. Override... can be used even before
   * setting a default, but they will throw an exception if they conflict with
   * a previous override.
   *
   * To access a setting, use the [] operator (repeatedly for nested settings).
   * To set a default, use SetDefault..., e.g. for setting the default
   * YFS->MODE in the main settings singleton, use
   *
   *   auto& s = Settings::GetMainSettings();
   *   s["YFS"]["MODE"].SetDefault(1);
   *
   * Most functions support method chaining (i.e. they return a reference of
   * the called object), so if you want to read the setting immediately after
   * setting its default:
   *
   *   auto Nc = s["NCOLOURS"].SetDefault(3).GetScalar<int>();
   *
   * A default can be set repeatedly, but only to the same value. If a default
   * for a setting is set to one value, and later to a different value, an
   * exception is raised.
   *
   * To get around this, Get...WithOtherDefault can be called, which does the
   * same as SetDefault...().Get<T>() the above, but another default is passed
   * as an argument which is used instead of a default set by SetDefault...
   * earlier.  This can be useful e.g. for subclasses, that need to use a
   * different default for some reason.
   *
   * As we do not provide a caching functionality here, the retrieval of
   * defaults is rather expensive. Therefore make sure to set defaults and load
   * settings only during the set-up of the framework and certainly not from
   * within event loops.
   *
   * To access a list of subsettings, use `GetItems()`, e.g. like this:
   *
   *   for (const Scoped_Setting& s : mainsettings["PROCESSES"].GetItems()) {
   *     ...
   *   }
   *
   * When settings defaults and/or overrides, lists are ignored, i.e. for the
   * previous example, with
   * mainsettings["PROCESSES"][proc_name]["NIn"].SetDefault...(2), each
   * s[proc_name]["NIn"] will have its default set to 2.
   */
  class Settings {

    friend Scoped_Settings;
    friend Settings_Writer;

  public:

    static void InitializeMainSettings(int argc, char* argv[]);
    static void FinalizeMainSettings();
    static Settings& GetMainSettings();

    Settings();
    Settings(const std::string& yaml);

    Scoped_Settings operator[](const std::string& scope);

    std::vector<std::string> GetKeys();

    /// convenience function to declare (many) settings with empty defaults
    void DeclareVectorSettingsWithEmptyDefault(const std::vector<std::string>& keys);
    void DeclareMatrixSettingsWithEmptyDefault(const std::vector<std::string>& keys);

    /**
     * global tags will be replaced as-is, i.e. each occurrence of key will be
     * replaced by value; only values read as numbers are affected by this
     * replacement; an example is E_CMS; global tags are not defined by the
     * user
     */
    void AddGlobalTag(const std::string& key, const std::string& value);

    /**
     * tags will be replaced if they are found within $(...), i.e. each
     * occurrence of $(key) will be replaced by value; examples are NJETS or
     * SCF; all tags are defined by the user (on the command line or in the
     * config file)
     */
    void AddTag(const std::string& key, const std::string& value);
    const String_Map& GetTags()       { return m_tags; }
    void ReplaceTags(std::string&);

    std::string GetPath();
    String_Vector GetConfigFiles();

  private:

    static std::unique_ptr<Settings> mainsettings;

    // all defaults (and overrides) are stored as a matrix of strings
    typedef String_Vector Defaults_Key;
    typedef String_Matrix Defaults_Value;
    typedef std::map<Defaults_Key, Defaults_Value> Defaults;
    Defaults m_defaults;
    Defaults m_overrides;

    // keep records for writing out settings
    std::map<Defaults_Key, std::set<Defaults_Value>> m_usedvalues;
    std::map<Defaults_Key, std::unordered_set<std::string>> m_otherscalardefaults;

    std::map<Defaults_Key, String_Map> m_replacements;
    std::map<Defaults_Key, String_Vector> m_defaultsynonyms;
    std::map<Defaults_Key, String_Vector> m_synonyms;

    String_Map m_tags;
    String_Map m_globaltags;

    // all setting readers in their order of precedence (i.e. first one with a
    // match wins)
    std::vector<std::unique_ptr<Yaml_Reader>> m_yamlreaders;

    Algebra_Interpreter m_interpeter;
    bool m_interpreterenabled{ true };

    Settings(int argc, char* argv[]);

    std::vector<std::string> GetKeys(const Settings_Keys&);

    bool SetInterpreterEnabled(bool b) {
      const auto wasenabled = m_interpreterenabled;
      m_interpreterenabled = b;
      return wasenabled;
    }

    /**
     * the Get... functions already interprete strings as numbers; however, if
     * a value is first accessed as a string (e.g. when a list of different
     * value types read as a string vector first), Interprete can be used to
     * do this a-posteriori
     */
    template <typename T>
    T Interprete(std::string);

    void DeclareVectorSettingsWithEmptyDefault(const std::vector<std::string>& keys,
                                               const Settings_Keys&);
    void DeclareMatrixSettingsWithEmptyDefault(const std::vector<std::string>& keys,
                                               const Settings_Keys&);

    template <typename T>
    void SetDefault(const Settings_Keys&, const T& value);
    void SetDefault(const Settings_Keys&, const char* value);
    template <typename T>
    void SetDefault(const Settings_Keys&,
                    const std::vector<T>& values);
    template <typename T>
    void SetDefaultMatrix(const Settings_Keys&,
                          const std::vector<std::vector<T>>& values);
    template <typename T>
    void SetDefaultMatrix(const std::vector<std::string>& keys,
                          const std::vector<std::vector<T>>& values);

    bool IsCustomised(const Settings_Keys&);

    template <typename T>
    void OverrideScalar(const Settings_Keys&, const T& value);
    template <typename T>
    void OverrideVector(const Settings_Keys&, const std::vector<T>& values);
    template <typename T>
    void OverrideMatrix(const Settings_Keys&,
                        const std::vector<std::vector<T>>& values);

    std::string GetScalarDefault(const Defaults_Key&, const Defaults&);
    String_Vector GetVectorDefault(const Defaults_Key&, const Defaults&);
    String_Matrix GetMatrixDefault(const Defaults_Key&, const Defaults&);

    template <typename T>
    T GetScalar(const Settings_Keys&);
    template <typename T>
    std::vector<T> GetVector(const Settings_Keys&);
    template <typename T>
    std::vector<std::vector<T>> GetMatrix(const Settings_Keys&);

    template <typename T>
    T GetScalarWithOtherDefault(const Settings_Keys&, const T& otherdefault);

    template <typename T>
    void SetReplacementList(const Settings_Keys&,
                            const std::map<std::string, T>& list);

    void SetDefaultSynonyms(const Settings_Keys&,
                            const String_Vector& synonyms);

    void SetSynonyms(const Settings_Keys&,
                     const String_Vector& synonyms);

    bool IsScalar(const Settings_Keys&);
    bool IsList(const Settings_Keys&);
    bool IsMap(const Settings_Keys&);
    size_t GetItemsCount(const Settings_Keys&);

    template <typename T>
    T Convert(const Settings_Keys&, const std::string&);

    bool IsDefaultSynonym(const Settings_Keys&,
                          const std::string& value);

    std::string ApplyReplacements(const Settings_Keys&,
                                  const std::string& value);

    // check if path is absolute, only works on Unix
    inline bool is_absolute(const std::string &s) {return (s.find("/") ==0);}
  };

  template <typename T>
  void Settings::SetDefault(const Settings_Keys& keys,
                                  const T& value)
  {
    SetDefaultMatrix<T>(keys, {{value}});
  }

  template <typename T>
  void Settings::SetDefault(const Settings_Keys& keys,
                            const std::vector<T>& values)
  {
    SetDefaultMatrix<T>(keys, {values});
  }

  template <typename T>
  void Settings::SetDefaultMatrix(const Settings_Keys& keys,
                                  const std::vector<std::vector<T> >& values)
  {
    SetDefaultMatrix(keys.IndizesRemoved(), values);
  }

  template <typename T>
  void Settings::SetDefaultMatrix(const std::vector<std::string>& keys,
                                  const std::vector<std::vector<T> >& values)
  {
    String_Matrix stringvalues;
    for (const auto& valuerow : values) {
      std::vector<std::string> stringvaluerow;
      for (const auto& value : valuerow) {
        stringvaluerow.push_back(ToString(value));
      }
      stringvalues.push_back(stringvaluerow);
    }
    const auto it = m_defaults.find(keys);
    if (it != m_defaults.end()) {
      if (it->second != stringvalues) {
        THROW(fatal_error,
               "The default value for "
               + keys.back() + " is already set to a different value.");
      } else {
        return;
      }
    }
    m_defaults[keys] = stringvalues;
  }

  template <typename T>
  void Settings::OverrideScalar(const Settings_Keys& settings_keys,
                                const T& value)
  {
    OverrideVector<T>(settings_keys, {value});
  }

  template <typename T>
  void Settings::OverrideVector(const Settings_Keys& settings_keys,
                                const std::vector<T>& values)
  {
    OverrideMatrix<T>(settings_keys, {values});
  }

  template <typename T>
  void Settings::OverrideMatrix(const Settings_Keys& settings_keys,
                                const std::vector<std::vector<T>>& values)
  {
    String_Matrix stringvalues;
    for (const auto& valuerow : values) {
      std::vector<std::string> stringvaluerow;
      for (const auto& value : valuerow) {
        stringvaluerow.push_back(ToString(value));
      }
      stringvalues.push_back(stringvaluerow);
    }
    const std::vector<std::string> keys{ settings_keys.IndizesRemoved() };
    const auto it = m_overrides.find(keys);
    if (it != m_overrides.end()) {
      if (it->second != stringvalues) {
        THROW(fatal_error,
               "The override for "
               + keys.back() + " is already set to a different value.");
      } else {
        return;
      }
    }
    m_overrides[keys] = stringvalues;
  }

  template <typename T>
  T Settings::GetScalar(const Settings_Keys& settings_keys)
  {
    Defaults_Key keys{ settings_keys.IndizesRemoved() };
    std::string defaultvalue;
    defaultvalue = GetScalarDefault(keys, m_defaults);
    std::string value;
    if (m_overrides.find(keys) != m_overrides.end()) {
      value = GetScalarDefault(keys, m_overrides);
    } else {
      const auto it = m_synonyms.find(keys);
      for (auto& reader : m_yamlreaders) {
        value = reader->GetScalar<std::string>(settings_keys);
        if (value.empty() && it != m_synonyms.end()) {
          auto synonym_settings_keys = settings_keys;
          for (const auto& synonym : it->second) {
            synonym_settings_keys.back() = synonym;
            value = reader->GetScalar<std::string>(synonym_settings_keys);
            if (!value.empty()) {
              keys = synonym_settings_keys.IndizesRemoved();
              break;
            }
          }
        }
        if (!value.empty())
          break;
      }
    }
    if (value.empty() || IsDefaultSynonym(settings_keys, value)) {
      value = defaultvalue;
    }
    const auto convertedvalue = Convert<T>(settings_keys, value);
    if (value.empty())
      m_usedvalues[keys].insert(String_Matrix{{""}});
    else
      m_usedvalues[keys].insert(String_Matrix{{ToString<T>(convertedvalue)}});
    return Convert<T>(settings_keys, value);
  }

  template <typename T>
  std::vector<T> Settings::GetVector(const Settings_Keys& settings_keys)
  {
    Defaults_Key keys{ settings_keys.IndizesRemoved() };
    String_Vector defaultvalues;
    defaultvalues = GetVectorDefault(keys, m_defaults);
    std::vector<std::string> values;
    if (m_overrides.find(keys) != m_overrides.end()) {
      values = GetVectorDefault(keys, m_overrides);
    } else {
      const auto it = m_synonyms.find(keys);
      for (auto& reader : m_yamlreaders) {
        values = reader->GetVector<std::string>(settings_keys);
        if (values.empty() && it != m_synonyms.end()) {
          auto synonym_settings_keys = settings_keys;
          for (const auto& synonym : it->second) {
            synonym_settings_keys.back() = synonym;
            values = reader->GetVector<std::string>(synonym_settings_keys);
            if (!values.empty()) {
              keys = synonym_settings_keys.IndizesRemoved();
              break;
            }
          }
        }
        if (!values.empty())
          break;
      }
    }
    if (values.empty())
      values = defaultvalues;
    std::vector<T> convertedvalues;
    std::vector<std::string> convertedstringvalues;
    for (std::vector<std::string>::const_iterator it(values.begin());
         it != values.end();
         ++it) {
      convertedvalues.push_back(Convert<T>(settings_keys, *it));
      convertedstringvalues.push_back(ToString<T>(convertedvalues.back()));
    }
    if (values.empty())
      m_usedvalues[keys].insert(String_Matrix{{}});
    else
      m_usedvalues[keys].insert(String_Matrix{convertedstringvalues});
    return convertedvalues;
  }

  template <typename T>
  std::vector<std::vector<T> > Settings::GetMatrix(
      const Settings_Keys& settings_keys)
  {
    Defaults_Key keys{ settings_keys.IndizesRemoved() };
    String_Matrix defaultvalues;
    defaultvalues = GetMatrixDefault(keys, m_defaults);
    String_Matrix values;
    if (m_overrides.find(keys) != m_overrides.end()) {
      values = GetMatrixDefault(keys, m_overrides);
    } else {
      const auto it = m_synonyms.find(keys);
      for (auto& reader : m_yamlreaders) {
        values = reader->GetMatrix<std::string>(settings_keys);
        if (values.empty() && it != m_synonyms.end()) {
          auto synonym_settings_keys = settings_keys;
          for (const auto& synonym : it->second) {
            synonym_settings_keys.back() = synonym;
            values = reader->GetMatrix<std::string>(synonym_settings_keys);
            if (!values.empty()) {
              keys = synonym_settings_keys.IndizesRemoved();
              break;
            }
          }
        }
        if (!values.empty())
          break;
      }
    }
    if (values.empty())
      values = defaultvalues;
    std::vector<std::vector<T>> convertedvalues;
    String_Matrix convertedstringvalues;
    for (const auto& valuerow : values) {
      std::vector<T> convertedvaluerow;
      String_Vector convertedstringvaluerow;
      for (const auto& value : valuerow) {
        convertedvaluerow.push_back(Convert<T>(settings_keys, value));
        convertedstringvaluerow.push_back(ToString<T>(convertedvaluerow.back()));
      }
      convertedvalues.push_back(convertedvaluerow);
      convertedstringvalues.push_back(convertedstringvaluerow);
    }
    m_usedvalues[keys].insert(convertedstringvalues);
    return convertedvalues;
  }

  template <typename T>
  T Settings::GetScalarWithOtherDefault(const Settings_Keys& settings_keys,
                                        const T& otherdefault)
  {
    const Defaults_Key keys{ settings_keys.IndizesRemoved() };
    // temporarily replace m_defaults entry before calling GetScalar
    const auto it = m_defaults.find(keys);
    String_Matrix oldvalue;
    bool wasempty{ false };
    if (it == m_defaults.end()) {
      wasempty = true;
    } else {
      oldvalue = it->second;
      m_defaults.erase(it);
    }
    SetDefault(settings_keys, otherdefault);
    const T value{ GetScalar<T>(settings_keys) };
    if (wasempty)
      m_defaults.erase(m_defaults.find(keys));
    else
      m_defaults[keys] = oldvalue;
    m_otherscalardefaults[keys].insert(ToString<T>(otherdefault));
    return value;
  }

  template <typename T>
  void Settings::SetReplacementList(const Settings_Keys& settings_keys,
                                    const std::map<std::string, T>& list)
  {
    const Defaults_Key keys{ settings_keys.IndizesRemoved() };
    const auto it = m_replacements.find(keys);
    if (m_replacements.find(keys) != m_replacements.end())
      if (list != it->second)
        THROW(fatal_error, "A different scalar replacement list for "
                           + keys.back() + " has already been set.");
    m_replacements[keys] = list;
  }


  template <typename T>
  T Settings::Convert(const Settings_Keys& settings_keys,
                      const std::string& value)
  {
    std::string convertedvalue{ value };
    ReplaceTags(convertedvalue);
    convertedvalue = ApplyReplacements(settings_keys, convertedvalue);
    return Interprete<T>(convertedvalue);
  }

  template <typename T>
  T Settings::Interprete(std::string value)
  {
    T dummy;
    if (typeid(dummy) == typeid(int)
        || typeid(dummy)==typeid(unsigned int)
        || typeid(dummy)==typeid(long)
        || typeid(dummy)==typeid(float)
        || typeid(dummy)==typeid(double)
        || typeid(dummy)==typeid(unsigned long long)) {
      value = ReplaceUnits(value);
      if (m_interpreterenabled)
        value = m_interpeter.Interprete(value);
    }
    return ToType<T>(value);
  }

}

#endif
