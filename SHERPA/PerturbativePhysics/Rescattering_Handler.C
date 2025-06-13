#include "SHERPA/PerturbativePhysics/Rescattering_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace SHERPA;
using namespace ALPACA;
using namespace ATOOLS;
using namespace std;

Rescattering_Handler::Rescattering_Handler() :
  m_name(std::string("None")) {
  if (true) {
    p_alpaca = new Alpaca();
    m_name = string("Alpaca");
  }
}

Rescattering_Handler::~Rescattering_Handler() {}

// member functions
Return_Value::code Rescattering_Handler::operator()(Blob_List * blobs) {
  Return_Value::code rval = Return_Value::Nothing;
  if (p_alpaca) {
    msg_Out()<<METHOD<<" to check "<<blobs->size()<<" blobs.\n";
    if ((*p_alpaca)(blobs)) {
      rval = Return_Value::Success;
      //msg_Out()<<"YES!\n"
	    //   <<(*blobs->back())<<"\n";
      msg_Out() << "YES!\n" <<"\n";
    } else{
      msg_Out() << METHOD << ": ERROR: !(*p_alpaca)(blobs), will exit" << endl;
      //exit(1.);
    }
  }
  return rval;
}
