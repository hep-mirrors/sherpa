#include "ALPACA/Main/Alpaca.H"
#include "ATOOLS/Org/Message.H"

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;


Alpaca::Alpaca() {}

bool Alpaca::operator()(ATOOLS::Blob_List * blobs) {
  if (!Harvest(blobs)) exit(1); 
}

bool Alpaca::Harvest(ATOOLS::Blob_List * blobs) {
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    if ((*bit)->Has(blob_status::needs_rescattering))
      msg_Out()<<(**bit)<<"\n";
  }
  return false;
}
