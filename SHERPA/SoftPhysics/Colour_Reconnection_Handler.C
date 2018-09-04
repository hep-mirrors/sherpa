#include "SHERPA/SoftPhysics/Colour_Reconnection_Handler.H"

using namespace SHERPA;
using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Colour_Reconnection_Handler::
Colour_Reconnection_Handler(const string path,const string file) :
  p_reconnections(NULL)
{
  Default_Reader reader;
  reader.AddIgnore("[");
  reader.AddIgnore("]");
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  p_reconnections = new Reconnection_Handler();
  p_reconnections->Initialize(&reader);
}

Colour_Reconnection_Handler::~Colour_Reconnection_Handler() {
  msg_Out()<<METHOD<<"\n";
  if (p_reconnections) delete p_reconnections;
}

Return_Value::code Colour_Reconnection_Handler::operator()(Blob_List *const blobs) {
  if (!p_reconnections) return Return_Value::Nothing;
  return (*p_reconnections)(blobs);
}

void Colour_Reconnection_Handler::CleanUp(const size_t & mode) {
  p_reconnections->Reset();
}
