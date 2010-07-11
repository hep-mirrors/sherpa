#ifndef ATOOLS_Phys_Blob_List_H
#define ATOOLS_Phys_Blob_List_H

#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle_List.H"
#include <deque>
#include <set>
#include <iostream>

namespace ATOOLS {

  class Blob_List: public std::deque<Blob*> {
  private:

    class Destructor {
    private:

      Blob_List *p_list;

    public:

      // destructor
      Destructor(): p_list(NULL) {}
      Destructor(Blob_List *const list): p_list(list) {}

      // destructor
      inline ~Destructor() { if (p_list!=NULL) p_list->Clear(); }

    };// end of class Destructor

    Destructor m_destructor;

    friend std::ostream &operator<<(std::ostream &str,const Blob_List &list);

    void FindConnected(Blob *blob,Blob_List &connected,
		       std::set<const Blob*> &selected);

    void DeleteConnected(Blob *blob,std::set<Blob*> &deleted);

    bool TotalFourMomentum(Blob *blob,std::set<Blob*> &summed,
			   Vec4D &inisum,Vec4D &finsum,const int mode) const;

  public:

    // constructor
    Blob_List();
    Blob_List(const bool destruct);

    // member functions
    Vec4D TotalFourMomentum() const;
    Vec4D IncomingFourMomentum() const;
    Vec4D OutgoingFourMomentum() const;

    bool FourMomentumConservation() const;
    bool ColorConservation() const;

    Blob *FindFirst(const btp::code code) const;
    Blob *FindLast(const btp::code code) const;

    Blob_List Find(const btp::code code) const;

    Blob_List FindConnected(const Blob *blob);
    Blob_List FindConnected(const Particle *blob);

    Blob *AddBlob(const btp::code &type);

    bool Delete(Blob *blob);

    size_t DeleteConnected(Blob *blob);
    size_t DeleteConnected(Particle *particle);

    Particle_List ExtractParticles(const int status,
				   const int mode=0) const;
    Particle_List ExtractLooseParticles(const int mode=0) const;

    bool MergeSubsequentType(btp::code,btp::code,long int &,long int &);
    void MergeSubsequentTypeRecursively(btp::code,btp::code,
					long int &,long int &);
    void ShortenRecursively(kf_code,
			    btp::code=btp::Unspecified,btp::code=btp::Unspecified);
    bool Shorten(kf_code,btp::code,btp::code);
    void Clear(Blob *blob=NULL);

    Blob_List Copy() const;

    double Weight() const;

  };// end of class Blob_List
  
  std::ostream &operator<<(std::ostream &str,const Blob_List &list);

}// end of namespace ATOOLS

#endif
