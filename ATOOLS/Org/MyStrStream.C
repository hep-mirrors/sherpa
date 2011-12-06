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

  std::string ReplaceUnits(std::string v)
  {
    std::string f("");
    size_t l(v.length());
    for (size_t i(0);i<l;) {
      if (v[i]==' ' || v[i]=='\t');
      else if (v[i]=='k') f+=(i+1<l && v[i+1]=='B')?"*(1<<10)":"*1000";
      else if (v[i]=='M') f+=(i+1<l && v[i+1]=='B')?"*(1<<20)":"*1000000";
      else if (v[i]=='G') f+=(i+1<l && v[i+1]=='B')?"*(1<<30)":"*1000000000";
      else {
	++i;
	continue;
      }
      v.replace(i,(i+1<l && v[i+1]=='B')?2:1,f);
      f="";
    }
    return v;
  }

}// end of namespace ATOOLS;
