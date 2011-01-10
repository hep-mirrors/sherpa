#ifndef ATOOLS_Org_SVN_Info_H
#define ATOOLS_Org_SVN_Info_H

#include <string>
#include <vector>

namespace ATOOLS {

  class SVN_Info {
  private:

    std::string m_name, m_branch, m_revision;

    static std::vector<const SVN_Info*> *s_objects;

  public:

    SVN_Info(const std::string &name,
	     const std::string &branch,
	     const std::string &revision);
    ~SVN_Info();

    inline const std::string &Name() const     { return m_name;     }
    inline const std::string &Branch() const   { return m_branch;   }
    inline const std::string &Revision() const { return m_revision; }

    inline static const std::vector<const SVN_Info*> *
    Infos() { return s_objects; }

  };// end of class SVN_Info

}// end of namespace ATOOLS

#endif
