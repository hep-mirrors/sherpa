#ifndef Lund_Interface_H
#define Lund_Interface_H

#include "Event_Phase_Handler.H"
#include "Blob_List.H"
#include "My_Map.H"

namespace SHERPA {

  class Lund_Interface: public Event_Phase_Handler {
  private:

    std::string m_path, m_file;

    bool ConvertParticles(std::map<int,ATOOLS::Particle*> &converted);
    bool ConstructBlobs(ATOOLS::Blob_List *blobs,std::map<int,ATOOLS::Particle*> &converted);
    bool SetTypes(ATOOLS::Blob_List *blobs);
    bool DeleteObsolete(std::map<int,ATOOLS::Particle*> &converted);
    bool ConvertEvent(ATOOLS::Blob_List *blobs);

  public:

    // constructor
    Lund_Interface(std::string _m_path,std::string _m_file);

    // destructor
    ~Lund_Interface();

    // member functions
    void CleanUp();

    bool Treat(ATOOLS::Blob_List *, double &);

  };// end of class Lund_Interface

}// end of namespace SHERPA

#endif
