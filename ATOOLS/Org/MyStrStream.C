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

  template <class Type> Type Number(std::string v)
  {
    Type f(1);
    size_t l(v.length());
    for (size_t i(0);i<l;) {
      if (v[i]==' ' || v[i]=='\t');
      else if (v[i]=='k') f*=(i+1<l && v[i+1]=='B')?(1<<10):1000;
      else if (v[i]=='M') f*=(i+1<l && v[i+1]=='B')?(1<<20):1000000;
      else if (v[i]=='G') f*=(i+1<l && v[i+1]=='B')?(1<<30):1000000000;
      else {
	++i;
	continue;
      }
      v.erase(i,1);
      if (v[i]=='B') v.erase(i,1);
    }
    return ToType<Type>(v)*f;
  }

  template int Number<int>(std::string v);
  template long int Number<long int>(std::string v);
  template double Number<double>(std::string v);

}// end of namespace ATOOLS;
