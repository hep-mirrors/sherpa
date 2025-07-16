#include "ATOOLS/Org/File_IO_Base.H"

#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

File_IO_Base::File_IO_Base(const unsigned int inputfiles,
			   const unsigned int outputfiles):
  m_infiles(std::vector<My_In_File>(inputfiles)),
  m_outfiles(std::vector<My_Out_File>(outputfiles)) {}

File_IO_Base::~File_IO_Base() 
{
}

