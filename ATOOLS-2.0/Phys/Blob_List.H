#ifndef Blob_List_H
#define Blob_List_H

#include "Blob.H"
#include <deque>
#include <set>
#include <iostream>

namespace ATOOLS {

  class Blob_List: public std::deque<Blob*> {
  private:

    friend std::ostream &operator<<(std::ostream &str,const Blob_List &list);

    void FindConnected(Blob *blob,Blob_List &connected,
		       std::set<const Blob*> &selected);

    void DeleteConnected(Blob *blob,std::set<Blob*> &deleted);

    void TotalFourMomentum(Blob *blob,std::set<Blob*> &summed,
			   Vec4D &inisum,Vec4D &finsum,const int mode) const;

  public:

    // member functions
    Vec4D TotalFourMomentum() const;
    Vec4D IncomingFourMomentum() const;
    Vec4D OutgoingFourMomentum() const;

    bool FourMomentumConservation() const;

    Blob *FindFirst(const btp::code code) const;
    Blob *FindLast(const btp::code code) const;

    Blob_List Find(const btp::code code) const;

    Blob_List FindConnected(const Blob *blob);
    Blob_List FindConnected(const Particle *blob);

    bool Delete(Blob *blob);

    size_t DeleteConnected(Blob *blob);
    size_t DeleteConnected(Particle *particle);
    
    void Clear();

  };// end of class Blob_List
  
  std::ostream &operator<<(std::ostream &str,const Blob_List &list);

}// end of namespace ATOOLS

#endif
