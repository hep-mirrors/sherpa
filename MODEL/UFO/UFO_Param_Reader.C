#include "MODEL/UFO/UFO_Param_Reader.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Read_Write_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <assert.h>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace UFO;

using std::string;
using std::vector;
using std::stringstream;

UFO_Param_Reader::UFO_Param_Reader(const string& filepath)
{
  // if filepath not explicitly given, use run card
  bool given(filepath!="");
  string file_path = given ? filepath :ATOOLS::rpa->gen.Variable("RUN_DATA_FILE");
  string filename(""), path("");
  // split filepath into path and name
  size_t pos(file_path.find_last_of("/"));
  if (pos!=std::string::npos){
    path=file_path.substr(0,pos+1);
    filename=file_path.substr(pos+1);
  }
  else{
    path=string("./");
    filename=file_path;
  }
  if (filename.find("|")!=std::string::npos) 
      filename=filename.substr(0,filename.find("|"));
  p_dataread = new ATOOLS::Data_Reader(" ",";","#","=");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(path);
  p_dataread->SetInputFile(filename);
  p_dataread->SetIgnoreCase(true);
  p_dataread->SetAddCommandLine(false);
}

UFO_Param_Reader::~UFO_Param_Reader()
{
  delete p_dataread;
}


template<class Read_Type> Read_Type 
UFO_Param_Reader::GetEntry(const string& block, const unsigned int& n, const unsigned int& m)
{
  p_dataread->SetFileBegin(string("block ").append(block));
  p_dataread->RereadInFile();
  vector< vector<string> > vals;
  p_dataread->MatrixFromFile(vals);
  for(vector< vector<string> >::const_iterator it = vals.begin(); it != vals.end(); ++it)
    {
      for(vector<string>::const_iterator jt = it->begin(); jt != it->end(); ++jt)
	if (p_dataread->Find(*jt, "block") != string::npos)
	  NotFound(block, n, m);
      if (it->size() < 3)
	continue;
      if (ATOOLS::ToType<int>((*it)[0]) == n && ATOOLS::ToType<int>((*it)[1]) == m)
	return ATOOLS::ToType<Read_Type>((*it)[2]);
    }
  NotFound(block, n, m);
}

template<class Read_Type> Read_Type 
UFO_Param_Reader::GetEntry(const string& block, const unsigned int& n)
{
  string lc_block(block);
  for (size_t i=0;i<lc_block.length();++i) lc_block[i]=tolower(lc_block[i]);
  // width specifications don't follow UFO conventions for some reason
  if (lc_block == "decay") 
    return GetWidth<Read_Type>(n);
  p_dataread->SetFileBegin(string("block ").append(block));
  p_dataread->RereadInFile();
  vector< vector<string> > vals;
  p_dataread->MatrixFromFile(vals);
  for(vector< vector<string> >::const_iterator it = vals.begin(); it != vals.end(); ++it)
    {
      for(vector<string>::const_iterator jt = it->begin(); jt != it->end(); ++jt)
	if (p_dataread->Find(*jt, "block") != string::npos)
	    NotFound(block, n);
      if (it->size() < 2)
	continue;
      if (ATOOLS::ToType<int>((*it)[0]) == n)
	return ATOOLS::ToType<Read_Type>((*it)[1]);
    }
  NotFound(block, n);
}

template<class Read_Type> Read_Type
UFO_Param_Reader::GetWidth(const unsigned int& n)
{
  p_dataread->ClearFileBegin();
  p_dataread->RereadInFile();
  vector< vector<string> > vals;
  p_dataread->MatrixFromFile(vals);
  for(vector< vector<string> >::const_iterator it = vals.begin(); it != vals.end(); ++it)
    {
      assert(it->size());
      if (p_dataread->Find((*it)[0], "decay") != string::npos && it->size() > 2 && ATOOLS::ToType<int>((*it)[1]) == n )
	  return ATOOLS::ToType<Read_Type>((*it)[2]);
    }
  NotFound(string("decay"),n);
}

void UFO_Param_Reader::NotFound(const string &block, const unsigned int& n, const unsigned int& m)
{
  stringstream message;
  message << ("Entry [") << n << "," << m << "] " << "in block " << block << " not found.";
  THROW(fatal_error, message.str().c_str() );
}

void UFO_Param_Reader::NotFound(const string &block, const unsigned int& n)
{
  stringstream message;
  message << ("Entry [") << n << "] " << "in block " << block << " not found.";
  THROW(fatal_error, message.str().c_str() );
}


namespace UFO
{
  template double UFO_Param_Reader::GetEntry(const string& block, const unsigned int& n, const unsigned int& m);
  template double UFO_Param_Reader::GetEntry(const string& block, const unsigned int& n);
  template int UFO_Param_Reader::GetEntry(const string& block, const unsigned int& n, const unsigned int& m);
  template int UFO_Param_Reader::GetEntry(const string& block, const unsigned int& n);
}
