#include "File_IO_Base.H"

using namespace ATOOLS;

File_IO_Base::File_IO_Base():
  m_inputpath("./"),
  m_inputfile(""),
  m_outputpath("./"),
  m_outputfile("") {}

File_IO_Base::~File_IO_Base() {}
