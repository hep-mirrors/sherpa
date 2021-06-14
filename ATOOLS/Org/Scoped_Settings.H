#ifndef ATOOLS_Org_ScopedSettings_H
#define ATOOLS_Org_ScopedSettings_H

#include "ATOOLS/Org/Settings.H"

#include <cassert>

namespace ATOOLS {

  class Scoped_Settings {

    friend Settings;

  public:

    Scoped_Settings();
    Scoped_Settings(const std::string& yamlstring);
    Scoped_Settings(const Scoped_Settings&);

    Scoped_Settings& operator=(Scoped_Settings);

    Scoped_Settings operator[](const std::string& scope) const;
    Scoped_Settings operator[](size_t scope) const;

    std::vector<std::string> GetKeys();

    void DeclareVectorSettingsWithEmptyDefault(const std::vector<std::string>& keys);
    void DeclareMatrixSettingsWithEmptyDefault(const std::vector<std::string>& keys);

    template <typename T>
    Scoped_Settings& SetDefault(const T&);
    template <typename T>
    Scoped_Settings& SetDefault(const std::vector<T>&);
    template <typename T>
    Scoped_Settings& SetDefault(std::initializer_list<T>);
    template <typename T>
    Scoped_Settings& SetDefault(const std::vector<std::vector<T>>&);
    template <typename T>
    Scoped_Settings& SetDefault(std::initializer_list<std::initializer_list<T>>);

    template <typename T>
    Scoped_Settings& SetReplacementList(const std::map<std::string, T>& list);

    Scoped_Settings& SetDefaultSynonyms(const String_Vector&);
    Scoped_Settings& SetSynonyms(const String_Vector&);

    /// set a common replacement list: [Off, false, 0, no] -> None
    Scoped_Settings& UseNoneReplacements();

    /// set a common replacement list: [None] -> max double
    Scoped_Settings& UseMaxDoubleReplacements();

    /// set a common replacement list: [None] -> 0
    Scoped_Settings& UseZeroReplacements();

    bool IsCustomised();

    void SetInterpreterEnabled(bool b) { m_interpreterenabled = b; }
    template <typename T>
    T Interprete(const std::string&);

    template <typename T>
    void OverrideScalar(const T&);
    template <typename T>
    void OverrideVector(const std::vector<T>&);
    template <typename T>
    void OverrideMatrix(const std::vector<std::vector<T>>&);

    template <typename T>
    T Get(); // alias for GetScalar
    template <typename T>
    T GetScalar();
    template <typename T>
    std::vector<T> GetVector();
    template <typename T>
    std::vector<std::vector<T>> GetMatrix();

    template <typename T>
    T GetScalarWithOtherDefault(const T& otherdefault);

    /**
     * convenience method to read two-valued settings
     *
     * this is often useful in the context of beam-specific settings
     * @param expandsinglevalues If true and only a single value is
     * encountered, this value is automatically repeated to return a two-valued
     * vector.
     */
    template <typename T>
    std::vector<T> GetTwoVector(bool expandsinglevalues=true);

    void ReplaceTags(std::string& s) { return m_rootsettings->ReplaceTags(s); }

    std::string GetPath() { return m_rootsettings->GetPath(); }

    bool IsScalar() const;
    bool IsList() const;
    bool IsMap() const;
    std::vector<Scoped_Settings> GetItems() const;
    size_t GetItemsCount();
    Scoped_Settings GetItemAtIndex(const size_t&) const;

    /**
     * obtain the index of the current item
     *
     * If the scope does not point to a sequence of sub-settings, an exception
     * is thrown.
     */
    size_t GetIndex() const;

  private:

    // use this if it is not otherwise guaranteed that rootsettings is kept
    // alive
    std::shared_ptr<Settings> m_ownedsettings;
    Settings* m_rootsettings;
    Settings_Keys m_scopes;

    Scoped_Settings(Settings &rootsettings, const std::string& scope);
    Scoped_Settings Scoped(const Setting_Key& scope) const;
    void AppendScope(const Setting_Key& scope);

    bool m_interpreterenabled{ true };

  };

  template <typename T>
  Scoped_Settings& Scoped_Settings::SetDefault(const T& value)
  {
    m_rootsettings->SetDefault(m_scopes, value);
    return *this;
  }

  template <typename T>
  Scoped_Settings& Scoped_Settings::SetDefault(const std::vector<T>& values)
  {
    m_rootsettings->SetDefault<T>(m_scopes, values);
    return *this;
  }

  template <typename T>
  Scoped_Settings& Scoped_Settings::SetDefault(std::initializer_list<T> values)
  {
    m_rootsettings->SetDefault<T>(m_scopes, std::vector<T>(values));
    return *this;
  }

  template <typename T>
  Scoped_Settings& Scoped_Settings::SetDefault(
      const std::vector<std::vector<T>>& values)
  {
    m_rootsettings->SetDefaultMatrix<T>(m_scopes, values);
    return *this;
  }

  template <typename T>
  Scoped_Settings& Scoped_Settings::SetDefault(std::initializer_list<std::initializer_list<T>> l)
  {
    std::vector<std::vector<T>> m;
    for (const auto& row : l)
      m.push_back(row);
    m_rootsettings->SetDefaultMatrix<T>(m_scopes, m);
    return *this;
  }

  template <typename T>
  Scoped_Settings& Scoped_Settings::SetReplacementList(
      const std::map<std::string, T>& list)
  {
    m_rootsettings->SetReplacementList(m_scopes, list);
    return *this;
  }

  template <typename T>
  void Scoped_Settings::OverrideScalar(const T& value)
  {
    m_rootsettings->OverrideScalar<T>(m_scopes, value);
  }

  template <typename T>
  void Scoped_Settings::OverrideVector(const std::vector<T>& values)
  {
    m_rootsettings->OverrideVector<T>(m_scopes, values);
  }

  template <typename T>
  void Scoped_Settings::OverrideMatrix(
      const std::vector<std::vector<T>>& values)
  {
    m_rootsettings->OverrideMatrix<T>(m_scopes, values);
  }

  template <typename T>
  T Scoped_Settings::Get()
  {
    return GetScalar<T>();
  }

  template <typename T>
  T Scoped_Settings::GetScalar()
  {
    auto wasenabled = m_rootsettings->SetInterpreterEnabled(m_interpreterenabled);
    const auto r = m_rootsettings->GetScalar<T>(m_scopes);
    m_rootsettings->SetInterpreterEnabled(wasenabled);
    return r;
  }

  template <typename T>
  std::vector<T> Scoped_Settings::GetVector()
  {
    auto wasenabled = m_rootsettings->SetInterpreterEnabled(m_interpreterenabled);
    const auto r = m_rootsettings->GetVector<T>(m_scopes);
    m_rootsettings->SetInterpreterEnabled(wasenabled);
    return r;
  }

  template <typename T>
  std::vector<std::vector<T>> Scoped_Settings::GetMatrix()
  {
    auto wasenabled = m_rootsettings->SetInterpreterEnabled(m_interpreterenabled);
    const auto r = m_rootsettings->GetMatrix<T>(m_scopes);
    m_rootsettings->SetInterpreterEnabled(wasenabled);
    return r;
  }

  template <typename T>
  T Scoped_Settings::GetScalarWithOtherDefault(const T& otherdefault)
  {
    auto wasenabled = m_rootsettings->SetInterpreterEnabled(m_interpreterenabled);
    const auto r =
      m_rootsettings->GetScalarWithOtherDefault<T>(m_scopes, otherdefault);
    m_rootsettings->SetInterpreterEnabled(wasenabled);
    return r;
  }

  template <typename T>
  std::vector<T> Scoped_Settings::GetTwoVector(bool expandsinglevalues)
  {
    auto values = GetVector<T>();
    if (values.size() == 1 && expandsinglevalues) {
      values.push_back(values[0]);
    } else if (values.size() != 2) {
      THROW(fatal_error,
            "The setting " + m_scopes.Name() + " must hold "
            + (expandsinglevalues ? "one or two" : "two")
            + " values.");
    }
    assert(values.size() == 2);
    return values;
  }

  template <typename T>
  T Scoped_Settings::Interprete(const std::string& value)
  {
    return m_rootsettings->Interprete<T>(value);
  }

}

#endif
