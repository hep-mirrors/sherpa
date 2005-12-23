#include "Read_Write_Base.H"

#include "MyStrStream.H"
#include "MathTools.H"

using namespace ATOOLS;

std::vector<std::string> Read_Write_Base::s_commandline;

Read_Write_Base::Read_Write_Base(const unsigned int infiles,
				 const unsigned int outfiles):
  File_IO_Base(infiles,outfiles),
  m_comment(std::vector<std::string>(1,defaultcom)), 
  m_ignore(std::vector<std::string>(1,defaultcut)),
  m_separator(std::vector<std::string>(1,defaultsep))
{
  Init();
}

Read_Write_Base::Read_Write_Base(const unsigned int infiles,
				 const unsigned int outfiles,
				 const std::string _m_cut,
				 const std::string _m_separator,
				 const std::string _m_comment):
  File_IO_Base(infiles,outfiles),
  m_comment(std::vector<std::string>(1,_m_comment)), 
  m_ignore(std::vector<std::string>(1,_m_cut)),
  m_separator(std::vector<std::string>(1,_m_separator))
{
  Init();
}

Read_Write_Base::~Read_Write_Base() 
{
  for (unsigned int i=0;i<m_infiles.size();++i) CloseInFile(i,true);
  m_infiles.clear();
  for (unsigned int i=0;i<m_outfiles.size();++i) CloseOutFile(i,true);
  m_outfiles.clear();
  delete p_interpreter;
}

void Read_Write_Base::Init()
{
  p_interpreter = new Algebra_Interpreter();
  m_blank.push_back(defaultblank);
  m_blank.push_back(defaulttab);
  m_vectortype=vtc::vertical;
  m_matrixtype=mtc::normal; 
  m_allownans=false;
  m_addcommandline=true;
  m_ignorecase=false;
  m_ignoreblanks=false;
  m_exactmatch=true;
  m_interprete=false;
  m_cmode=false;
  m_occurrence=std::string::npos;
  m_escape='\\';
}

size_t Read_Write_Base::Find(std::string input,std::string parameter,
			     size_t &length) const
{
  if (m_ignorecase) {
    for (size_t i=0;i<input.length();++i) input[i]=toupper(input[i]);
    for (size_t i=0;i<parameter.length();++i) 
      parameter[i]=toupper(parameter[i]);
  }
  size_t cutinputblanks=0;
  if (m_ignoreblanks) {
    for (size_t j=0;j<Blank().size();++j) {
      bool lastblank=true;
      for (size_t i=0;i<input.length();++i) {
	if (input[i]==Blank()[j]) {
	  input[i]=Blank()[0];
	  if (lastblank) {
	    input=input.substr(0,i)+input.substr(i+1,input.length());
	    ++cutinputblanks;
	  }
	  lastblank=true;
	}
	else {
	  lastblank=false;
	}
      }
    }
    for (size_t j=0;j<Blank().size();++j) {
      bool lastblank=true;
      for (size_t i=0;i<parameter.length();++i) {
	if (parameter[i]==Blank()[j]) {
	  parameter[i]=Blank()[0];
	  if (lastblank) {
	    parameter=parameter.substr(0,i)+
	      parameter.substr(i,parameter.length());
	  }
	  lastblank=true;
	}
	else {
	  lastblank=false;
	}
      }
    }
  }
  length=parameter.length()+cutinputblanks;
  size_t pos=input.find(parameter);
  if (m_exactmatch && pos!=0) {
    size_t i=0;
    for (;i<Blank().size();++i) if (input[pos-1]==Blank()[i]) break;
    if (i==Blank().size()) pos=std::string::npos;
  }
  if (pos==std::string::npos) length=0;
  return pos;
}

std::string &Read_Write_Base::Interprete(std::string &lastline) 
{
  if (!m_cmode) return lastline;
  if (!m_ifresults.empty()) {
    size_t pos=std::string::npos;
    for (size_t i=0;i<lastline.length();++i)
      if (lastline[i]=='{') ++m_ifresults.back().second;
      else if (lastline[i]=='}') {
	if (m_ifresults.back().second>1) --m_ifresults.back().second;
	else {
	  pos=i; 
	  break;
	}
      }
    if (pos!=std::string::npos) {
      if (m_ifresults.back().first) lastline.replace(pos,1,"");
      else lastline=lastline.substr(pos+1);
      m_ifresults.pop_back();
      return Interprete(lastline);
    }
    else {
      if (!m_ifresults.back().first) lastline="";
    }
  }
  size_t pos=Min(lastline.find("if "),lastline.find("if("));
  if (pos>0 && 
      lastline[pos-1]!=32 && lastline[pos-1]!=9) pos=std::string::npos;
  if (pos!=std::string::npos) {
    size_t opos=lastline.find("(",pos), cpos=lastline.find(")",pos);
    if (cpos!=std::string::npos && opos<cpos) {
      Algebra_Interpreter interpreter;
      interpreter.SetTagReplacer(this);
      bool result=ToType<int>
	(interpreter.Interprete(lastline.
				substr(opos+1,cpos-opos-1)));
      opos=lastline.find("{",cpos);
      for (opos=cpos+1;opos<lastline.length();++opos) 
	if (lastline[opos]=='{') {
	  m_ifresults.push_back(std::pair<int,size_t>(result,1));
	  ++opos;
	  break;
	}
	else {
	  bool blank=false;
	  for (size_t i=0;i<m_blank.size();++i) 
	    if (lastline[opos]==m_blank[i]) {
	      blank=true;
	      break;
	    }
	  if (!blank) break;
	}
      if (opos==lastline.length()) opos=cpos;
      if (result) lastline=lastline.substr(opos);
      else lastline="";
    }
  }
  return lastline;
}

bool Read_Write_Base::OpenInFile(const unsigned int i)
{  
  if (InputFile(i)==nullstring) {
    if (m_addcommandline && CommandLine().size()>0) {
      AddFileContent(CommandLine());
      return true;
    }
    return false;
  }
  if (InFileMode(i)==fom::unknown) SetInFileMode(fom::temporary);
  if (m_infiles[i]()==NULL) {
    m_infiles[i].Open();	
    m_filecontent.clear();
    std::string lastline;
    bool checkbegin=(bool)(m_filebegin.size()>0);
    bool checkend=(bool)(m_fileend.size()>0);
    int filebegin=0;
    unsigned int occurrence=0;
    if (*m_infiles[i]) {
      getline(*m_infiles[i],lastline);
      do {
	if (checkbegin) {
	  for (size_t length=0,j=0;j<m_filebegin.size();++j) {
	    size_t pos=Find(lastline,m_filebegin[j],length);
	    if (pos!=std::string::npos) {
	      if (filebegin==0) lastline=lastline.substr(pos+length);
	      if (occurrence==m_occurrence ||
		  m_occurrence==std::string::npos) ++filebegin;
	      if (filebegin==0) ++occurrence;
	      break;
	    }
	  }
	  if (filebegin==0) {
	    lastline=ATOOLS::nullstring;
	  }
	  else if (checkend) {
	    for (size_t length=0,j=0;j<m_fileend.size();++j) {
	      size_t pos=Find(lastline,m_fileend[j],length);
	      if (pos!=std::string::npos) {
		if (occurrence==m_occurrence ||
		    m_occurrence==std::string::npos) --filebegin;
		if (filebegin==0) {
		  lastline=lastline.substr(0,pos);
		  ++occurrence;
		}
		break;
	      }
	    }
	  }
	  Interprete(lastline);
	  if (lastline.length()>0) m_filecontent.push_back(lastline);
	}
	else if (lastline.length()>0) {
	  Interprete(lastline);
	  if (lastline.length()>0) m_filecontent.push_back(lastline);
	}
	getline(*m_infiles[i],lastline);
      } while (*m_infiles[i]);
    }
  }
  bool success=m_filecontent.size()>0;
  if (m_addcommandline) AddFileContent(CommandLine());
  if (!success) m_infiles[i].SetMode(fom::error);
  return success;
}

bool Read_Write_Base::OpenOutFile(const unsigned int i)
{  
  if (OutputFile(i)==nullstring) return false;
  if (OutFileMode(i)==fom::unknown) SetOutFileMode(fom::permanent);
  if (m_outfiles[i]()==NULL) {
    m_outfiles[i].Open();	
    if (m_filebegin.size()>0 && !m_outfiles[i]->bad()) {
      (*m_outfiles[i])<<m_filebegin[0]<<std::endl;
    }
  }
  return !m_outfiles[i]->bad();
}

void Read_Write_Base::CloseInFile(const unsigned int i,const bool force)
{ 
  if (m_infiles[i]()==NULL) return;
  if (m_infiles[i].Mode()==fom::permanent && !force) return;
  m_filecontent.clear();
  m_infiles[i].Close(); 
}

void Read_Write_Base::CloseOutFile(const unsigned int i,const bool force)
{ 
  if (m_outfiles[i]()==NULL) return;
  if (m_outfiles[i].Mode()==fom::permanent && !force) return;
  if (m_fileend.size()>0 && !m_outfiles[i]->bad()) {
    (*m_outfiles[i])<<m_fileend[0]<<std::endl;
  }
  m_outfiles[i].Close(); 
}

std::string Read_Write_Base::ReplaceTags(std::string &expr) const
{ 
  std::string tag=expr;
  bool success=false;
  for (std::map<std::string,std::string>::const_iterator 
	 tit=m_tags.begin();tit!=m_tags.end();++tit) {
    size_t pos=tag.find(tit->first);
    if (pos!=std::string::npos) {
      tag.replace(pos,tit->first.length(),tit->second);
      success=true;
    }
  }
  if (success && tag!=expr) return ReplaceTags(tag);
  return tag;
}

namespace ATOOLS {

  template <typename Type> Type Read_Write_Base::Default() 
  { 
    return std::numeric_limits<Type>::max(); 
  }
  
  template int Read_Write_Base::Default<int>();
  template unsigned int Read_Write_Base::Default<unsigned int>();
  template long int Read_Write_Base::Default<long int>();
  template float Read_Write_Base::Default<float>();
  template double Read_Write_Base::Default<double>();
  
  template <> std::string Read_Write_Base::Default<std::string>()
  {
    return "";
  }
  
}

