#include "File_IO_Base.H"

#include "Message.H"

using namespace ATOOLS;

File_IO_Base::File_IO_Base(const unsigned int inputfiles,
			   const unsigned int outputfiles):
  m_inputpath(std::vector<std::string>(inputfiles)),
  m_inputfile(std::vector<std::string>(inputfiles)),
  m_outputpath(std::vector<std::string>(outputfiles)),
  m_outputfile(std::vector<std::string>(outputfiles)), 
  m_infile(std::vector<std::ifstream*>(inputfiles)),
  m_outfile(std::vector<std::ofstream*>(outputfiles)),
  m_infilemode(std::vector<OpenModeID>(inputfiles)),
  m_outfilemode(std::vector<OpenModeID>(outputfiles)) 
{
  for (unsigned int i=0;i<m_infile.size();m_infile[i++]=NULL);
  for (unsigned int i=0;i<m_outfile.size();m_outfile[i++]=NULL);
  for (unsigned int i=0;i<m_infilemode.size();m_infilemode[i++]=Unknown);
  for (unsigned int i=0;i<m_outfilemode.size();m_outfilemode[i++]=Unknown);
}

File_IO_Base::~File_IO_Base() 
{
  // note this always calls File_IO_Base::Close..File()
  //      ie. it is not possible to use virtual function can not 
  //      be called inside the constructor since the at this point
  //      the constructor of derived classes has already finished
  for (unsigned int i=0;i<m_infile.size();++i) CloseInFile(i,true);
  for (unsigned int i=0;i<m_outfile.size();++i) CloseOutFile(i,true);
}

bool File_IO_Base::OpenInFile(const unsigned int i)
{
  ATOOLS::msg.Error()<<"File_IO_Base::OpenInFile("<<i<<"): "
		     <<"Virtual function called!"<<std::endl;
  return false;
}

bool File_IO_Base::OpenOutFile(const unsigned int i)
{
  ATOOLS::msg.Error()<<"File_IO_Base::OpenOutFile("<<i<<"): "
		     <<"Virtual function called!"<<std::endl;
  return false;
}

void File_IO_Base::CloseInFile(const unsigned int i,const bool force)
{
  if (m_infile[i]==NULL) return;
  ATOOLS::msg.Error()<<"File_IO_Base::CloseInFile("<<i<<"): "
		     <<"Virtual function called!"<<std::endl;
  return;
}

void File_IO_Base::CloseOutFile(const unsigned int i,const bool force)
{
  if (m_outfile[i]==NULL) return;
  ATOOLS::msg.Error()<<"File_IO_Base::CloseOutFile("<<i<<"): "
		     <<"Virtual function called!"<<std::endl;
  return;
}

void File_IO_Base::SetInputPath(const std::string _m_inputpath,
				const unsigned int i)
{ 
  if (_m_inputpath!=m_inputpath[i]) {
    CloseInFile(i,true);
    m_inputpath[i]=_m_inputpath; 
  }
}

const std::string File_IO_Base::InputPath(const unsigned int i) const
{ 
  if (i<m_inputfile.size()) return m_inputpath[i]; 
  return nullstring; 
}

void File_IO_Base::SetInputFile(const std::string _m_inputfile,
				const unsigned int i)
{ 
  if (_m_inputfile!=m_inputfile[i]) {
    CloseInFile(i,true);
    m_inputfile[i]=_m_inputfile; 
  }
}

const std::string File_IO_Base::InputFile(const unsigned int i) const
{ 
  if (i<m_inputfile.size()) return m_inputfile[i]; 
  return nullstring; 
}

void File_IO_Base::SetOutputPath(const std::string _m_outputpath,
				 const unsigned int i)
{ 
  if (_m_outputpath!=m_outputpath[i]) {
    CloseOutFile(i,true);
    m_outputpath[i]=_m_outputpath; 
  }
}

const std::string File_IO_Base::OutputPath(const unsigned int i) const
{ 
  if (i<m_outputpath.size()) return m_outputpath[i]; 
  return nullstring; 
}

void File_IO_Base::SetOutputFile(const std::string _m_outputfile,
				 const unsigned int i)
{ 
  if (_m_outputfile!=m_outputfile[i]) {
    CloseOutFile(i,true);
    m_outputfile[i]=_m_outputfile; 
  }
}

const std::string File_IO_Base::OutputFile(const unsigned int i) const
{ 
  if (i<m_outputfile.size()) return m_outputfile[i]; 
  return nullstring; 
}

void File_IO_Base::SetInFileMode(const OpenModeID _m_infilemode,
				 const unsigned int i)
{
  m_infilemode[i]=_m_infilemode; 
  CloseInFile(i);
}

void File_IO_Base::SetOutFileMode(const OpenModeID _m_outfilemode,
				  const unsigned int i)
{ 
  m_outfilemode[i]=_m_outfilemode; 
  CloseOutFile(i);
}

bool File_IO_Base::CheckInputPath(const unsigned int i)
{	
  if (i>=m_inputpath.size()) return false;
  if (m_inputpath[i]==nullstring) {
    msg_Tracking()<<"File_IO_Base::CheckInputPath(): "
		  <<"No input path specified!"<<std::endl;
    return false;
  }
  return true;
}

bool File_IO_Base::CheckInputFile(const unsigned int i)
{	
  if (i>=m_inputfile.size()) return false;
  if (m_inputfile[i]==nullstring) {
    msg_Tracking()<<"File_IO_Base::CheckInputFile(): "
		  <<"No input file specified!"<<std::endl;
    return false;
  }
  return true;
}

bool File_IO_Base::CheckOutputPath(const unsigned int i)
{	
  if (i>=m_outputpath.size()) return false;
  if (m_outputpath[i]==nullstring) {
    msg_Tracking()<<"File_IO_Base::CheckOutputPath(): "
		  <<"No output path specified!"<<std::endl;
    return false;
  }
  return true;
}

bool File_IO_Base::CheckOutputFile(const unsigned int i)
{	
  if (i>=m_outputpath.size()) return false;
  if (m_outputfile[i]==nullstring) {
    msg_Tracking()<<"File_IO_Base::CheckOutputFile(): "
		  <<"No output file specified!"<<std::endl;
    return false;
  }
  return true;
}


