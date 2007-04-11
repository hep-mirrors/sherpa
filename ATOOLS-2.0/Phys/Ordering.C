#include "Ordering.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE ATOOLS::Order_Base
#define PARAMETER_TYPE std::string
#include "Getter_Function.C"

#include "Message.H"
#include "Exception.H"

using namespace ATOOLS;

Order_Base::~Order_Base() {}

bool Order_Base::operator()(const Vec4D &a,const Vec4D &b) const
{
  THROW(fatal_error,"Virtual function called");
}

bool Order_Base::operator()(const Particle &a,const Particle &b) const
{
  THROW(fatal_error,"Virtual function called");
}

bool Order_Base::operator()(Particle * const &a,Particle * const &b) const
{
  THROW(fatal_error,"Virtual function called");
}

void Order_Base::ShowOrders(const int mode)
{
  if (!msg.LevelIsInfo() || mode==0) return;
  msg.Out()<<"Order_Base::ShowOrders(): {\n\n";
  Order_Getter::PrintGetterInfo(msg.Out(),20);
  msg.Out()<<"\n}"<<std::endl;
}

template <class Class>
Order_Base *const GetOrder(const std::string &parameter)
{									
  return new Class();
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Order_Base *								\
  NAME::operator()(const std::string &parameter) const			\
  { return GetOrder<CLASS>(parameter); }

#define DEFINE_PRINT_METHOD(NAME,PRINT)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<PRINT; }

#define DEFINE_ORDER_GETTER(CLASS,NAME,TAG,PRINT)			\
  DECLARE_GETTER(NAME,TAG,Order_Base,std::string);			\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME,PRINT)

//-------------------------------------------------------------------------------

class Order_Up_E : public Order_Base {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return dabs(a[0])>dabs(b[0]); }
  bool operator()(const Particle &a,const Particle &b) const 
  { return dabs(a.Momentum()[0])>dabs(b.Momentum()[0]); }
  bool operator()(Particle * const &a,Particle * const &b) const 
  { return dabs(a->Momentum()[0])>dabs(b->Momentum()[0]); }
};

DEFINE_ORDER_GETTER(Order_Up_E,Order_Up_E_Getter,"E_UP","order E ascending")

//-------------------------------------------------------------------------------

class Order_Up_ET : public Order_Base {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return a.EPerp()>b.EPerp(); }
  bool operator()(const Particle &a,const Particle &b) const 
  { return a.Momentum().EPerp()>b.Momentum().EPerp(); }
  bool operator()(Particle * const &a,Particle * const &b) const 
  { return a->Momentum().EPerp()>b->Momentum().EPerp(); }
};

DEFINE_ORDER_GETTER(Order_Up_ET,Order_Up_ET_Getter,"ET_UP","order ET ascending")

//-------------------------------------------------------------------------------

class Order_Up_PT : public Order_Base {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return a.PPerp2()>b.PPerp2(); }
  bool operator()(const Particle &a,const Particle &b) const 
  { return a.Momentum().PPerp2()>b.Momentum().PPerp2(); }
  bool operator()(Particle * const &a,Particle * const &b) const 
  { return a->Momentum().PPerp2()>b->Momentum().PPerp2(); }
};

DEFINE_ORDER_GETTER(Order_Up_PT,Order_Up_PT_Getter,"PT_UP","order PT ascending")

//-------------------------------------------------------------------------------

class Order_Up_Eta : public Order_Base {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return dabs(a.Eta())>dabs(b.Eta()); }
  bool operator()(const Particle &a,const Particle &b) const 
  { return dabs(a.Momentum().Eta())>dabs(b.Momentum().Eta()); }
  bool operator()(Particle * const &a,Particle * const &b) const 
  { return dabs(a->Momentum().Eta())>dabs(b->Momentum().Eta()); }
};

DEFINE_ORDER_GETTER(Order_Up_Eta,Order_Up_Eta_Getter,"ETA_UP","order eta ascending")
