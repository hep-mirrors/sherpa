#include "ATOOLS/Org/MyStrStream.H"
#include <vector>

namespace ATOOLS {

  std::string StringReplace(const std::string &original,
                            const std::string &from, const std::string &to)
  { 
    if(from.length()==0) return original;
    std::string result=original;
    std::vector<int> matches;
    int pos=result.find(from);
    while(pos!=-1) {
      matches.push_back(pos);
      pos=result.find(from,pos+1);
    }

    int offset=0;
    size_t total=matches.size();
    int diff=to.length()-from.length();
    int fromlength=from.length();
    for(size_t i=0;i<total;++i) {
      result.erase(matches[i]+offset,fromlength);
      result.insert(matches[i]+offset,to);
      offset+=diff;
    }
    return result;
  }

}// end of namespace ATOOLS;
