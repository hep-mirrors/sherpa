#include "ATOOLS/Org/Yaml_Reader.H"

//#include "ATOOLS/Org/My_Limits.H"
//#include "ATOOLS/Org/Message.H"
//#include "ATOOLS/Org/MyStrStream.H"
//#include <typeinfo>
//#include <ctype.h>

using namespace ATOOLS;

Yaml_Reader::Yaml_Reader() {
  m_node = SHERPA_YAML::Node();
}

void Yaml_Reader::Read_File(char* fname) {
  m_node = SHERPA_YAML::LoadFile(fname);
}

Yaml_Reader::~Yaml_Reader() {
}


