#include "Read_Write_Base.H"

//#define DEBUG__Read_Write_Base
#ifdef DEBUG__Read_Write_Base
#include <iostream>
#endif

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
#ifdef DEBUG__Read_Write_Base
  std::cout<<"~Read_Write_Base("<<m_infile.size()<<","<<m_outfile.size()<<")\n";
#endif
  for (unsigned int i=0;i<m_infile.size();++i) CloseInFile(i,true);
  m_infile.clear();
  for (unsigned int i=0;i<m_outfile.size();++i) CloseOutFile(i,true);
  m_outfile.clear();
}

void Read_Write_Base::Init()
{
#ifdef DEBUG__Read_Write_Base
  std::cout<<" Read_Write_Base:\n";
#endif
  m_blank.push_back(defaultblank);
  m_blank.push_back(defaulttab);
  m_vectortype=VVertical;
  m_matrixtype=MNormal; 
  m_allownans=false;
  m_addcommandline=true;
  m_ignorecase=false;
  m_ignoreblanks=false;
  m_exactmatch=true;
  m_occurrence=std::string::npos;
}

size_t Read_Write_Base::Find(std::string input,std::string parameter,size_t &length) const
{
#ifdef DEBUG__Read_Write_Base1
  std::cout<<"Read_Write_Base::Find("<<input<<","<<parameter<<"): "<<std::endl;
#endif
  if (m_ignorecase) {
    for (size_t i=0;i<input.length();++i) input[i]=toupper(input[i]);
    for (size_t i=0;i<parameter.length();++i) parameter[i]=toupper(parameter[i]);
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
	    parameter=parameter.substr(0,i)+parameter.substr(i,parameter.length());
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
#ifdef DEBUG__Read_Write_Base1
  std::cout<<"   input     = '"<<input<<"'("<<cutinputblanks<<")\n"
	   <<"   parameter = '"<<parameter<<"' at "<<(int)pos<<"\n"
	   <<"   exact     = "<<(pos!=std::string::npos)<<std::endl;
#endif
  if (pos==std::string::npos) length=0;
  return pos;
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
  if (InFileMode(i)==Unknown) SetInFileMode(Temporary);
  if (m_infile[i]==NULL) {
#ifdef DEBUG__Read_Write_Base
    std::cout<<"Read_Write_Base::OpenInFile("<<i<<"): "
	     <<"Opening file '"<<m_inputpath[i]+m_inputfile[i]<<"'."<<std::endl;
#endif
    m_infile[i]=new std::ifstream();	
    m_infile[i]->open((m_inputpath[i]+m_inputfile[i]).c_str()); 
    m_filecontent.clear();
    std::string lastline;
    bool checkbegin=(bool)(m_filebegin.size()>0);
    bool checkend=(bool)(m_fileend.size()>0);
    int filebegin=0;
    unsigned int occurrence=0;
    if (*m_infile[i]) {
      getline(*m_infile[i],lastline);
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
	  if (lastline.length()>0) m_filecontent.push_back(lastline);
	}
	else {
	  m_filecontent.push_back(lastline);
	}
	getline(*m_infile[i],lastline);
      } while (*m_infile[i]);
    }
  }
  bool sucess=m_filecontent.size()>0;
  if (m_addcommandline) AddFileContent(CommandLine());
  //  std::cout<<" filecontent.size()="<<m_filecontent.size()<<std::endl;
  return sucess;
  //  return !m_infile[i]->bad();
}

bool Read_Write_Base::OpenOutFile(const unsigned int i)
{  
  if (OutputFile(i)==nullstring) return false;
  if (OutFileMode(i)==Unknown) SetOutFileMode(Permanent);
  if (m_outfile[i]==NULL) {
#ifdef DEBUG__Read_Write_Base
    std::cout<<"Read_Write_Base::OpenOutFile("<<i<<"): "
	     <<"Opening file '"<<m_outputpath[i]+m_outputfile[i]<<"'."<<std::endl;
#endif
    m_outfile[i]=new std::ofstream();	
    m_outfile[i]->open((m_outputpath[i]+m_outputfile[i]).c_str());
    if (m_filebegin.size()>0 && !m_outfile[i]->bad()) {
      (*m_outfile[i])<<m_filebegin[0]<<std::endl;
    }
  }
  return !m_outfile[i]->bad();
}

void Read_Write_Base::CloseInFile(const unsigned int i,const bool force)
{ 
  if (m_infile[i]==NULL) return;
//   std::cout<<"calling Read_Write_Base::CloseInFile("<<i<<","<<force<<"): "
// 	   <<"Closing file '"<<m_inputpath[i]+m_inputfile[i]<<"'."<<m_infilemode[i]<<std::endl;
  if ((m_infilemode[i]==Permanent)&&(!force)) return;
#ifdef DEBUG__Read_Write_Base
  std::cout<<"Read_Write_Base::CloseInFile("<<i<<","<<force<<"): "
	   <<"Closing file '"<<m_inputpath[i]+m_inputfile[i]<<"'."<<m_infilemode[i]<<m_infile[i]<<std::endl;
#endif
  m_filecontent.clear();
  m_infile[i]->close(); 
  delete m_infile[i]; 
  m_infile[i]=NULL;
}

void Read_Write_Base::CloseOutFile(const unsigned int i,const bool force)
{ 
  if (m_outfile[i]==NULL) return;
  if ((m_outfilemode[i]==Permanent)&&(!force)) return;
#ifdef DEBUG__Read_Write_Base
  std::cout<<"Read_Write_Base::CloseOutFile("<<i<<","<<force<<"): "
	   <<"Closing file '"<<m_outputpath[i]+m_outputfile[i]<<"'."<<std::endl;
#endif
  if (m_fileend.size()>0 && !m_outfile[i]->bad()) {
    (*m_outfile[i])<<m_fileend[0]<<std::endl;
  }
  m_outfile[i]->close(); 
  delete m_outfile[i]; 
  m_outfile[i]=NULL;
}


