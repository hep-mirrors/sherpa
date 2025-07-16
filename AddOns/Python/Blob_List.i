//%module Blob_List
%include "std_deque.i"
%{
#include <ATOOLS/Phys/Blob_List.H>
#include <ATOOLS/Org/MyStrStream.H>
%}

%template(BlobDeque) std::deque<ATOOLS::Blob*>;

namespace ATOOLS {

  class Blob_List: public std::deque<Blob*> {
  public:

    // constructor
    Blob_List();
    Blob_List(const bool destruct);

    // member functions
    Vec4D TotalFourMomentum() const;
    Vec4D IncomingFourMomentum() const;
    Vec4D OutgoingFourMomentum() const;

    %extend {
      const Blob &GetFirst(const int code) const
      {
	return *self->FindFirst((btp::code)code);
      }
    }

    %extend{
      const Blob &__getitem__(unsigned int i){
	return *(*self)[i];
      };
    };

    %extend {
      std::string __str__() {
	MyStrStream conv;
	conv<<*self;
	return conv.str();
      };
    };

  };// end of class Blob_List
  
}// end of namespace ATOOLS
