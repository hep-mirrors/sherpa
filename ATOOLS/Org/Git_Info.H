#ifndef ATOOLS_Org_Git_Info_H
#define ATOOLS_Org_Git_Info_H

#include <string>
#include <map>

namespace ATOOLS {

  class Git_Info {
  private:

    std::string m_name, m_branch, m_revision, m_checksum;

    static std::map<const std::string,const Git_Info*> *s_objects;

    static bool s_check;

  public:

    Git_Info(const std::string &name,
	     const std::string &branch,
	     const std::string &revision,
	     const std::string &checksum);
    ~Git_Info();

    inline const std::string &Name() const     { return m_name;     }
    inline const std::string &Branch() const   { return m_branch;   }
    inline const std::string &Revision() const { return m_revision; }
    inline const std::string &Checksum() const { return m_checksum; }

    inline static const std::map<const std::string,const Git_Info*> *
    Infos() { return s_objects; }

    inline static void SetCheck(const bool &check) { s_check=check; }

  };// end of class Git_Info

}// end of namespace ATOOLS

#endif
