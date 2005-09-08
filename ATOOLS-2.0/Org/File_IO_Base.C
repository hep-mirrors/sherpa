#include "File_IO_Base.H"

#include "Message.H"

using namespace ATOOLS;

File_IO_Base::File_IO_Base(const unsigned int inputfiles,
			   const unsigned int outputfiles):
  m_infiles(std::vector<My_In_File>(inputfiles)),
  m_outfiles(std::vector<My_Out_File>(outputfiles)) {}

File_IO_Base::~File_IO_Base() 
{
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
  ATOOLS::msg.Error()<<"File_IO_Base::CloseInFile("<<i<<"): "
		     <<"Virtual function called!"<<std::endl;
  return;
}

void File_IO_Base::CloseOutFile(const unsigned int i,const bool force)
{
  ATOOLS::msg.Error()<<"File_IO_Base::CloseOutFile("<<i<<"): "
		     <<"Virtual function called!"<<std::endl;
  return;
}

