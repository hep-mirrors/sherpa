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
  return;
}

void File_IO_Base::CloseOutFile(const unsigned int i,const bool force)
{
  return;
}

