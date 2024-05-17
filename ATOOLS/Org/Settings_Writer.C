#include "ATOOLS/Org/Settings_Writer.H"

#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"

#include <algorithm>
#include <cassert>
#include <set>
#include <utility>

using namespace ATOOLS;

void Settings_Writer::WriteSettings(Settings& s)
{
  // check for settings that have not been used
  MyStrStream unused;
  bool did_find_unused {false};
  std::map<std::string, std::set<std::string>> unused_toplevel_keys;
  for (const auto& reader : s.m_yamlreaders) {
    // Ignore Decaydata.yaml, which is not under the control of the user,
    // and would just clutter the not-used section of the Settings report,
    // and/or confuse the user by the unused-settings WARNING issued below.
    if (reader->Name().rfind("Decaydata.yaml") != std::string::npos)
      continue;
    bool did_print_file_header {false};
    auto keys_vec = reader->AllSettingsKeys();
    for (const auto& keys : keys_vec) {
      // if a used key is a prefix of `keys`, we consider `keys` to have been
      // used, too
      const auto it = std::find_if(
          s.m_usedvalues.begin(), s.m_usedvalues.end(),
          [&keys](const std::pair<Settings_Keys,
                                  std::set<Settings::Defaults_Value>> &pair) {
            return keys.IsBeginningOf(pair.first) ||
                   pair.first.EndsWithParentIndexFor(keys);
          });

      if (it == s.m_usedvalues.end()) {
        did_find_unused = true;
        if (!did_print_file_header) {
          unused << "### " << reader->Name() << '\n';
          did_print_file_header = true;
        }
        unused << "- ";
        unused << keys << " = ";
        const auto vals = reader->GetFlattenedStringVectorWithDelimiters(keys);
        if (vals.size() == 1)
          unused << vals[0];
        else
          unused << vals;
        unused << '\n';
        unused_toplevel_keys[reader->Name()].insert(keys[0].GetName());
      }
    }
    if (did_print_file_header)
      unused << '\n';
  }
  unused << '\n';
  if (did_find_unused) {
    msg_Info()
        << om::brown << om::bold << "WARNING" << om::reset
        << ": Some settings that have been defined in the input files and/or\n"
        << "  the command line have not been used. The following top-level settings\n"
        << "  have not been used or include subsettings that have not been used:\n";
    for (const auto &kv : unused_toplevel_keys) {
      msg_Out() << "  - " << kv.first << ": ";
      for (const auto &key : kv.second)
        msg_Out() << key << " ";
      msg_Out() << "\n";
    }
    msg_Out()
        << "  For more details, inspect the Settings Report: Run `make` within\n"
        << "  the Settings_Report directory (requires `pandoc` to be installed),\n"
        << "  then open 'Settings_Report.html' with your browser.\n";
  }


  // order output in rows of customised settings and uncustomised settings
  MyStrStream customised, uncustomised;
  for (const auto& keysetpair : s.m_usedvalues) {
    std::vector<String_Matrix> vals;
    const Settings_Keys keys {keysetpair.first};
    const auto finalvals = keysetpair.second;
    assert(!finalvals.empty());

    // put all values for the table rows in `vals`, if a value has multiple
    // entries, then these are separated by "-- AND --"; begin with the defaults
    vals.push_back(s.m_defaults[keysetpair.first.IndizesRemoved()]);

    // replace defaults with file provided defaults in the case of
    // Decaydata.yaml, which is outside of user control
    bool is_set_by_decaydata {false};
    for (const auto& reader : s.m_yamlreaders) {
      if (reader->Name().rfind("Decaydata.yaml") == std::string::npos)
        continue;
      if (!reader->IsParameterCustomised(keys))
        break;
      is_set_by_decaydata = true;
      vals.back() = reader->GetMatrix<std::string>(keys);
      break;
    }

    // take into account other alternative defaults
    const auto otherdefaultsit = s.m_otherscalardefaults.find(keysetpair.first.IndizesRemoved());
    for (const auto& v : s.m_otherscalardefaults[keysetpair.first.IndizesRemoved()]) {
      if (!vals.back().empty())
        vals.back().push_back({"-- AND --"});
      vals.back().push_back({v});
    }

    // figure out if the setting has been customized by the user
    bool iscustomised {false};
    if (is_set_by_decaydata) {
      for (const auto& reader : s.m_yamlreaders) {
        if (reader->Name().rfind("Decaydata.yaml") != std::string::npos)
          continue;
        Settings_Keys keys{ keysetpair.first };
        if (!reader->IsParameterCustomised(keys))
          continue;
        iscustomised = true;
      }
      if (iscustomised)
        iscustomised = (s.GetMatrix<std::string>(keysetpair.first) != vals[0]);
    }
    else {
      iscustomised =
          !(finalvals.size() == 1 && (*finalvals.begin() == vals[0]));
    }

    if (iscustomised) {
      vals.push_back(s.m_overrides[keysetpair.first.IndizesRemoved()]);
      Settings_Keys keys{ keysetpair.first };
      for (auto it = s.m_yamlreaders.rbegin();
           it != s.m_yamlreaders.rend();
           ++it) {
        auto new_vals =
          (*it)->GetFlattenedStringVectorWithDelimiters(keys);
        int maxdim = 0, curdim = 0;
        for (const auto& val : new_vals) {
          if (val == "{{") {
            maxdim = std::max(maxdim, ++curdim);
          } else if (val == "}}") {
            --curdim;
          }
        }
        String_Matrix pruned_new_vals;
        pruned_new_vals.emplace_back();
        for (const auto& val : new_vals) {
          if (val == "{{") {
            ++curdim;
          } else if (val == "}}") {
            if (curdim > 1 && curdim < maxdim) {
              pruned_new_vals.emplace_back();
            }
            --curdim;
            if (curdim == 0 && maxdim > 2 && &val != &new_vals.back()) {
              pruned_new_vals.back().push_back("-- AND --");
              pruned_new_vals.emplace_back();
            }
          } else {
            pruned_new_vals.back().push_back(val);
          }
        }
        vals.push_back(pruned_new_vals);
      }
      if (!finalvals.empty()) {
        vals.push_back(*finalvals.begin());
        for (auto it = ++finalvals.begin(); it != finalvals.end(); ++it) {
          vals.back().push_back({"-- AND --"});
          std::copy(it->begin(), it->end(), std::back_inserter(vals.back()));
        }
      }
    } else if (keysetpair.first[0] == "HADRON_DECAYS") {
      // if a HADRON_DECAYS setting is not customised, we omit its output
      // entirely, otherwise our Settings Report becomes very lengthy
      continue;
    }

    // write table body in the appropriate section
    MyStrStream& current = iscustomised ? customised : uncustomised;
    MyStrStream keystream;
    for (size_t i{ 0 }; i < keysetpair.first.size(); ++i) {
      keystream << keysetpair.first[i];
      if (i + 1 < keysetpair.first.size())
        keystream << ":";
    }
    current << EncodeForMarkdown(keystream.str());
    current << "| ";
    for (size_t i{ 0 }; i < vals.size(); ++i) {
      MyStrStream valstream;
      for (size_t j{ 0 }; j < vals[i].size(); ++j) {
        for (size_t k{ 0 }; k < vals[i][j].size(); ++k) {
          valstream << vals[i][j][k];
          if (k + 1 < vals[i][j].size())
            valstream << ", ";
        }
        if (j + 1 < vals[i].size())
          valstream << "\n";
      }
      current << EncodeForMarkdown(valstream.str()) << " | ";
      if (!iscustomised)
        break;
    }
    current << " |\n";
  }

  const auto path = rpa->gen.Variable("SHERPA_RUN_PATH") + "/Settings_Report";
  MakeDir(path, true);

  std::ofstream file(path + "/Settings_Report.md");
  file << "---\n";
  file << "title: Sherpa run-time settings\n";
  file << "date: " << rpa->gen.Timer().TimeString(0) << "\n";
  file << "...\n\n";

  if (did_find_unused) {
    file << "Unused settings\n";
    file << "-------------------\n";
    file << "Parameters that have never been read by Sherpa during its"
         << " run are listed here. If you did expect the setting to be used,"
         << " check its spelling, and note that Sherpa setting names are case"
         << "-sensitive.\n\n";
    file << unused.str();
  }

  file << "Customised settings\n";
  file << "-------------------\n";
  file << "The parameters listed here have been customised in some way and"
       << " they have been read by Sherpa during its run."
       << " The last column lists the actual value used"
       << " after taking into account all setting sources (default values,"
       << " overrides, input files and the command line).\n\n";
  file << "In some cases, an alternative default value is being used."
       << " These alternatives will be separated by \"`-- AND --`\" from the"
       << " standard default, which will always be listed on top.\n\n";
  file << "Note that parameters that can take on different values because they"
       << " are set within a list, for example `param: [{x: 1}, {x: 2}, ...]`,"
       << " will not appear in the config-file or command-line columns. They"
       << " will be listed in the final-value column, with each different value"
       << " separated by an \"`-- AND --`\" line.\n\n";

  file << "| parameter | default value | override by SHERPA";
  const auto files = s.GetConfigFiles();
  for (const auto& f : files)
     file << " | " << f;
  file << " | command line | final value |\n";
  file << "|-|-|-";
  for (int i {0}; i < files.size(); ++i)
    file << "|-";
  file << "|-|-|\n";
  file << customised.str();

  file << "Settings kept at their default value\n";
  file << "-------------------\n";
  file << "The parameter listed here have not been customised, but they have"
       << " been read by Sherpa during its run.\n\n";

  file << "| parameter | default value |\n";
  file << "|-|-|\n";
  file << uncustomised.str();
  file.close();

  std::ofstream cssfile(path + "/Style.css");
  cssfile << "tr:nth-of-type(odd) {"
          << "  background-color:#eef;"
          << "}";
  cssfile << "th {"
          << "  background-color:#fff;"
          << "}";
  cssfile.close();

  std::ofstream makefile(path + "/Makefile");
  makefile << "Settings_Report.html: Settings_Report.md\n"
           << "\tpandoc -s -o Settings_Report.html -c Style.css"
           << " Settings_Report.md\n\n"
	   << "Settings_Report.org: Settings_Report.md\n"
           << "\tpandoc -f markdown -t org -c Style.css -o Settings_Report.org"
           << " Settings_Report.md\n\n"
           << "clean:\n"
           << "\trm -f Settings_Report.html\n\n"
           << ".PHONY: clean";
  makefile.close();
}

std::string Settings_Writer::EncodeForMarkdown(const std::string &data) const
{
  std::string buffer;
  buffer.reserve(data.size());
  for(size_t pos = 0; pos != data.size(); ++pos) {
    switch(data[pos]) {
      case '\n':
        buffer.append("<br />");
        break;
      case '|':
      case '_':
      case '*':
      case '$':
      case '\\':
      case '`':
      case '{':
      case '}':
      case '[':
      case ']':
      case '(':
      case ')':
      case '#':
      case '+':
      case '-':
      case '.':
      case '!':
      case '<':
      case '>':
        buffer.append("\\");
      default:
        buffer.append(&data[pos], 1);
        break;
    }
  }
  return buffer;
}
