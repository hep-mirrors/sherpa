#include "Object_Definition_Base.H"
#include "Exception.H"

using namespace ANALYSIS;
using namespace ATOOLS;


Object_Definition_Data::Object_Definition_Data(std::string mode)
{
  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) 
    THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
}

Object_Definition_Data::~Object_Definition_Data() { 
  m_particles.Clear();
}

void Object_Definition_Data::ResetPList() { 
  m_particles.Clear(); 
}

void Object_Definition_Data::AddPToPList(Particle * const part) { 
  m_particles.push_back(part); 
}

void Object_Definition_Data::SortPList() { 
  std::sort(m_particles.begin(),m_particles.end(),(*p_order));
}

Object_Definition_Base::Object_Definition_Base(const std::string name,
					       const kf::code code,
					       const std::string order="ET_UP") :
  m_name(name), m_code(code), p_data(new Object_Definition_Data(order))
{ }

Object_Definition_Base::~Object_Definition_Base() {
  if (p_data) { delete p_data; p_data = NULL; }
}



template <class Class>
Object_Definition_Base *const GetObjectDefinition(const std::string &parameter)
{									
  return new Class();
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Object_Definition_Base *						\
  NAME::operator()(const std::string &parameter) const			\
  { return GetObjectDefinition<CLASS>(parameter); }

#define DEFINE_PRINT_METHOD(NAME,PRINT)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<PRINT; }

#define DEFINE_ORDER_GETTER(CLASS,NAME,TAG,PRINT)			\
  DECLARE_GETTER(NAME,TAG,Object_Definition_Base,std::string);		\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME,PRINT)

