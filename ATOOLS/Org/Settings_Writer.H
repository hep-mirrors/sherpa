#ifndef ATOOLS_Org_SettingsWriter_H
#define ATOOLS_Org_SettingsWriter_H

#include <string>

namespace ATOOLS {

  class Settings;

  class Settings_Writer {

  public:

    void WriteSettings(Settings&);

  private:

    std::string EncodeForMarkdown(const std::string &) const;

  };

}

#endif
