#ifndef ATOOLS_Phys_Blob_Data_H
#define ATOOLS_Phys_Blob_Data_H

#include <string>
#include <vector>
#include <map>

namespace ATOOLS {

  class Blob_Data_Base {
  private:
    static long int s_number;
  public:
    template <class Type> Type &Get();
    template <class Type> void Set(const Type &data);
    Blob_Data_Base();
    Blob_Data_Base(const Blob_Data_Base &base);
    virtual ~Blob_Data_Base();
    virtual std::ostream & operator>>(std::ostream &) const =0;
    virtual Blob_Data_Base* ClonePtr() { return NULL; }
  };

  template <class Type> class Blob_Data : public Blob_Data_Base {
    Type m_data;
  public:
    Blob_Data(const Type & d);
    Type &Get() { return m_data; }
    void Set(const Type & d) { m_data=d; }
    std::ostream & operator>>(std::ostream &) const;
    Blob_Data_Base* ClonePtr();
    ~Blob_Data();
  };

  std::ostream& operator<<(std::ostream&,const Blob_Data_Base &);

  typedef std::map<std::string,Blob_Data_Base *> String_BlobDataBase_Map;

  template <class Type>
  Blob_Data<Type>::Blob_Data(const Type & d) : m_data(d) {}

  template <class Type>
  std::ostream & Blob_Data<Type>::operator>>(std::ostream & s) const 
  {
    s<<m_data;
    return s;
  }

  template <class Type>
  Type &Blob_Data_Base::Get() 
  {
    return ((Blob_Data<Type>*)this)->Get();
  }

  template <class Type>
  void Blob_Data_Base::Set(const Type &data) 
  {
    return ((Blob_Data<Type>*)this)->Set(data);
  }

}

#endif
  



