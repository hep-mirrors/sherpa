#include "Read_Write_Base.H"

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

Read_Write_Base::~Read_Write_Base() {}

void Read_Write_Base::Init()
{
  m_blank.push_back(defaultblank);
  m_blank.push_back(defaulttab);
  m_vectortype=VVertical;
  m_matrixtype=MNormal; 
  m_allownans=false;
  m_addcommandline=true;
}

bool Read_Write_Base::OpenInFile(const unsigned int i)
{  
  if (InputFile(i)==nullstring) return false;
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
    bool checkbegin=(bool)(m_filebegin!=ATOOLS::nullstring);
    bool checkend=(bool)(m_fileend!=ATOOLS::nullstring);
    int filebegin=0;
    if (*m_infile[i]) {
      getline(*m_infile[i],lastline);
      do {
	if (checkbegin) {
	  if (lastline.find(m_filebegin)!=std::string::npos) {
	    if (filebegin==0) {
	      lastline=lastline.substr(lastline.find(m_filebegin)+m_filebegin.length());
	    }
	    ++filebegin;
	  }
	  else {
	    if (filebegin==0) lastline=ATOOLS::nullstring;
	  }
	  if (checkend && filebegin>0) {
	    if (lastline.find(m_fileend)!=std::string::npos) {
	      --filebegin;
	      if (filebegin==0) {
		lastline=lastline.substr(0,lastline.find(m_fileend));
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
  if (m_addcommandline) AddFileContent(CommandLine());
  return !m_infile[i]->bad();
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
    if (m_filebegin!=std::string("") && !m_outfile[i]->bad()) {
      (*m_outfile[i])<<m_filebegin<<std::endl;
    }
  }
  return !m_outfile[i]->bad();
}

void Read_Write_Base::CloseInFile(const unsigned int i,const bool force)
{ 
  if (m_infile[i]==NULL) return;
  if ((m_infilemode[i]==Permanent)&&(!force)) return;
#ifdef DEBUG__Read_Write_Base
  std::cout<<"Read_Write_Base::CloseInFile("<<i<<","<<force<<"): "
	   <<"Closing file '"<<m_inputpath[i]+m_inputfile[i]<<"'."<<m_infilemode[i]<<std::endl;
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
  if (m_fileend!=std::string("") && !m_outfile[i]->bad()) {
    (*m_outfile[i])<<m_fileend<<std::endl;
  }
  m_outfile[i]->close(); 
  delete m_outfile[i]; 
  m_outfile[i]=NULL;
}


