//%module Blob
%{
#include <ATOOLS/Phys/Blob.H>
#include <ATOOLS/Org/MyStrStream.H>
%}

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
  };

  template <class Type> class Blob_Data : public Blob_Data_Base {
    Type m_data;
  public:
    Blob_Data(const Type & d);
    Type &Get() { return m_data; }
    void Set(const Type & d) { m_data=d; }
    ~Blob_Data();
  };

  std::ostream& operator<<(std::ostream&,const Blob_Data_Base &);

  typedef std::map<std::string,Blob_Data_Base *> String_BlobDataBase_Map;

  class Blob {
  public:

    Blob(const Vec4D _pos = Vec4D(0.,0.,0.,0.), const int _id=-1);
    Blob(const Blob *,const bool);
    ~Blob();
    %extend{
        const Particle &InPart(const size_t i) const
	{
	  return *self->ConstInParticle(i);
        };
        const Particle &OutPart(const size_t i) const
	{
	  return *self->ConstOutParticle(i);
        };
    };

    inline bool Has(blob_status::code status)     { return (int(status&m_status)>0); }
    %extend{
        double __getitem__(const std::string name) 
	{
	  if ((*self)[name]==NULL) return 0.0;
	  return (*self)[name]->Get<double>();
        };
    };

    int      Id()                        const { return m_id; }
    int      Status()                    const { return m_status; }
    int      Beam()                      const { return m_beam; }
    int      NInP()                      const { return m_inparticles.size(); }
    int      NOutP()                     const { return m_outparticles.size(); }
 
    const Vec4D& Position()              const { return m_position; }
    const Vec4D& CMS()                   const { return m_cms_vec; }
    const btp::code& Type()              const { return m_type; }
    std::string const TypeSpec()         const { return m_typespec; }

    %extend {
      std::string __str__() {
	MyStrStream conv;
	conv<<*self;
	return conv.str();
      };
    };

  };// end of class Blob

  %template(Blob_Data_Double) Blob_Data<double>;

}
