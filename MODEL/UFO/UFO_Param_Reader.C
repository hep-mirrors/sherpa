#include "MODEL/UFO/UFO_Param_Reader.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <assert.h>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace UFO;

using std::string;
using std::vector;
using std::stringstream;
using ATOOLS::ToType;
using ATOOLS::Data_Reader;
using ATOOLS::rpa;

UFO_Param_Reader::UFO_Param_Reader(const string& filepath)
{
  // if filepath not explicitly given, use run card
  m_use_runcard = (filepath=="");
  string file_path = m_use_runcard ? rpa->gen.Variable("RUN_DATA_FILE") : filepath;
  string filename(""), path("");
  // split filepath into path and name
  size_t pos(file_path.find_last_of("/"));
  if (pos!=string::npos){
    path=file_path.substr(0,pos+1);
    filename=file_path.substr(pos+1);
  }
  else{
    path=string("");
    filename=file_path;
  }
  if (filename.find("|")!=string::npos) 
      filename=filename.substr(0,filename.find("|"));
  if(m_use_runcard) filename+="|(ufo){|}(ufo)";
  Data_Reader dataread(" ",";","#","=");
  dataread.AddWordSeparator("\t");
  dataread.SetInputPath(path);
  dataread.SetInputFile(filename);
  dataread.SetIgnoreCase(true);
  dataread.SetAddCommandLine(false);
  dataread.MatrixFromFile(m_lines);
}

template<class Read_Type> Read_Type 
UFO_Param_Reader::GetEntry(const string& block,
                           const unsigned int& n,
                           const unsigned int& m,
                           const Read_Type& def,
                           const bool& err)
{
  vector< vector<string> >::const_iterator line = FindBlock(block);
  for(++line; line!=m_lines.end(); ++line){
    if (line->empty()) continue;
    if (IgnoreCaseCompare((*line)[0],"block"))
      return NotFound<Read_Type>(block,n,m,def,err);
    if (line->size() < 3) continue;
    if (ToType<int>((*line)[0])==n && ToType<int>((*line)[1])==m)
      return ToType<Read_Type>((*line)[2]);
  }
  return NotFound<Read_Type>(block,n,m,def,err);
}

template<class Read_Type> Read_Type 
UFO_Param_Reader::GetEntry(const string& block,
                           const unsigned int& n,
                           const Read_Type& def,
                           const bool& err)
{
  // widths in UFO param cards are handled differently for some reason
  if (IgnoreCaseCompare(block, "decay")) return GetWidth<Read_Type>(n,def,err);
  vector< vector<string> >::const_iterator line = FindBlock(block);
  for(++line; line!=m_lines.end(); ++line){
    if (line->empty()) continue;
    if (IgnoreCaseCompare((*line)[0],"block"))
      return NotFound<Read_Type>(block,n,def,err);
    if (line->size() < 2) continue;
    if (ToType<int>((*line)[0]) == n)
      return ToType<Read_Type>((*line)[1]);
  }
  return NotFound<Read_Type>(block,n,def,err);
}

template<class Read_Type> Read_Type
UFO_Param_Reader::GetWidth(const unsigned int& n,
                           const Read_Type& def,
                           const bool& err)
{
  for(vector< vector<string> >::const_iterator line = m_lines.begin(); line != m_lines.end(); ++line){
    if (line->size() < 3) continue;
    if (IgnoreCaseCompare((*line)[0],"decay") && ToType<int>((*line)[1]) == n )
      return ToType<Read_Type>((*line)[2]);
  }
  return NotFound<Read_Type>(string("decay"),n,def,err);
}

vector< vector<string> >::const_iterator UFO_Param_Reader::FindBlock(const string& block){
  vector< vector<string> >::const_iterator ret=m_lines.begin();
  msg_Debugging() << m_lines.size() << "\n";
  for(; ret!=m_lines.end(); ++ret){
    msg_Debugging() << *ret << "\n";
    if(ret->size()<2) continue;
    if(IgnoreCaseCompare((*ret)[1],block)) return ret;
  }
  THROW(fatal_error, "Block "+block+" not found");
  // avoid compiler warnings concerning missing return statement
  return m_lines.end();
}

bool UFO_Param_Reader::IgnoreCaseCompare(const std::string& a, const std::string& b){
  msg_Debugging() << a << " " << b << "\n";
  if (a.size() != b.size())return false;
  for (string::const_iterator ia = a.begin(), ib = b.begin(); ia!=a.end(); ++ia, ++ib)
    if (tolower(*ia) != tolower(*ib)) 
      return false;
  return true;
}

template<class Read_Type> Read_Type
UFO_Param_Reader::NotFound(const string &block,
                           const unsigned int& n, const unsigned int& m,
                           const Read_Type& def, const bool& err)
{
  stringstream message;
  message << ("Entry [") << n << "," << m << "] " << "in block " << block << " not found.";
  if (err) { THROW(fatal_error, message.str().c_str() ); }
  else
    msg_Error()<<METHOD<<"() {\n"
               <<"  "<<message.str().c_str()<<"\n"
               <<"  reverting to default value of "<<def
               <<", if this parameter is relevant for your calculation, \n"
               <<"  make sure it is implemented in your UFO model.\n"
               <<"}"<<std::endl;
  // avoid compiler warnings concerning missing return statement
  return def;
}

template<class Read_Type> Read_Type
UFO_Param_Reader::NotFound(const string &block, const unsigned int& n,
                           const Read_Type& def, const bool& err)
{
  stringstream message;
  message << ("Entry [") << n << "] " << "in block " << block << " not found.";
  if (err) { THROW(fatal_error, message.str().c_str() ); }
  else
    msg_Error()<<METHOD<<"() {\n"
               <<"  "<<message.str().c_str()<<"\n"
               <<"  reverting to default value of "<<def
               <<", if this parameter is relevant for your calculation, \n"
               <<"  make sure it is implemented in your UFO model.\n"
               <<"}"<<std::endl;
  // avoid compiler warnings concerning missing return statement
  return def;
}

namespace UFO
{
template Complex UFO_Param_Reader::GetEntry(const string& block,
                                           const unsigned int& n,
                                           const unsigned int& m,
                                           const Complex& def,
                                           const bool& err);
template Complex UFO_Param_Reader::GetEntry(const string& block,
                                           const unsigned int& n,
                                           const Complex& def,
                                           const bool& err);
template double UFO_Param_Reader::GetEntry(const string& block,
                                           const unsigned int& n,
                                           const unsigned int& m,
                                           const double& def,
                                           const bool& err);
template double UFO_Param_Reader::GetEntry(const string& block,
                                           const unsigned int& n,
                                           const double& def,
                                           const bool& err);
  template int UFO_Param_Reader::GetEntry(const string& block,
                                          const unsigned int& n,
                                          const unsigned int& m,
                                          const int& def,
                                          const bool& err);
  template int UFO_Param_Reader::GetEntry(const string& block,
                                          const unsigned int& n,
                                          const int& def,
                                          const bool& err);
}
