#ifndef Lund_Interface_H
#define Lund_Interface_H

#include "Particle_List.H"
#include "Blob_List.H"
#include <map>

namespace SHERPA {

  class Lund_Interface {
  private:

    std::string m_path, m_file;

    // tools for converting a whole event
    bool ConvertParticles(std::map<int,ATOOLS::Particle*> &converted);
    bool ConstructBlobs(ATOOLS::Blob_List *blobs,std::map<int,ATOOLS::Particle*> &converted);
    bool SetTypes(ATOOLS::Blob_List *blobs);
    bool DeleteObsolete(std::map<int,ATOOLS::Particle*> &converted);

    // tools for extracting hadrons
    void AddPartonToString(ATOOLS::Particle *particle,int &nhep);
    void FillPrimaryHadronsInBlob(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
				  ATOOLS::Particle_List *pl=NULL);
    void FillSecondaryHadronsInBlob(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
				    int daughter1,int daughter2,ATOOLS::Particle_List *pl=NULL);

  public:

    // constructor
    Lund_Interface(std::string _m_path,std::string _m_file);

    // destructor
    ~Lund_Interface();

    // member functions
    bool Hadronize(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,ATOOLS::Particle_List *pl=NULL);
    bool ConvertEvent(ATOOLS::Blob_List *bloblist);

  };// end of class Lund_Interface

}// end of namespace SHERPA

#endif
