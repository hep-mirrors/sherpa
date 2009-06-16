#include "ATOOLS/Phys/Ordering.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE ATOOLS::Order_Base
#define PARAMETER_TYPE std::string
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

Order_Base::~Order_Base() {}

void Order_Base::ShowOrders(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<"Order_Base::ShowOrders(): {\n\n";
  Order_Getter::PrintGetterInfo(msg_Out(),20);
  msg_Out()<<"\n}"<<std::endl;
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
  static bool OrderV(const Vec4D &a,const Vec4D &b)
  { return dabs(a[0])>dabs(b[0]); }
  static bool OrderP(const Particle &a,const Particle &b)
  { return dabs(a.Momentum()[0])>dabs(b.Momentum()[0]); }
  static bool OrderPP(Particle * const &a,Particle * const &b)
  { return dabs(a->Momentum()[0])>dabs(b->Momentum()[0]); }
  Order_Up_E(): Order_Base(OrderV,OrderP,OrderPP) {}
};

DEFINE_ORDER_GETTER(Order_Up_E,Order_Up_E_Getter,"E_UP","order E ascending")

//-------------------------------------------------------------------------------

class Order_Up_ET : public Order_Base {
public:
  static bool OrderV(const Vec4D &a,const Vec4D &b) 
  { return a.EPerp()>b.EPerp(); }
  static bool OrderP(const Particle &a,const Particle &b) 
  { return a.Momentum().EPerp()>b.Momentum().EPerp(); }
  static bool OrderPP(Particle * const &a,Particle * const &b) 
  { return a->Momentum().EPerp()>b->Momentum().EPerp(); }
  Order_Up_ET(): Order_Base(OrderV,OrderP,OrderPP) {}
};

DEFINE_ORDER_GETTER(Order_Up_ET,Order_Up_ET_Getter,"ET_UP","order ET ascending")

//-------------------------------------------------------------------------------

class Order_Up_PT : public Order_Base {
public:
  static bool OrderV(const Vec4D &a,const Vec4D &b) 
  { return a.PPerp2()>b.PPerp2(); }
  static bool OrderP(const Particle &a,const Particle &b) 
  { return a.Momentum().PPerp2()>b.Momentum().PPerp2(); }
  static bool OrderPP(Particle * const &a,Particle * const &b) 
  { return a->Momentum().PPerp2()>b->Momentum().PPerp2(); }
  Order_Up_PT(): Order_Base(OrderV,OrderP,OrderPP) {}
};

DEFINE_ORDER_GETTER(Order_Up_PT,Order_Up_PT_Getter,"PT_UP","order PT ascending")

//-------------------------------------------------------------------------------

class Order_Up_Eta : public Order_Base {
public:
  static bool OrderV(const Vec4D &a,const Vec4D &b) 
  { return dabs(a.Eta())>dabs(b.Eta()); }
  static bool OrderP(const Particle &a,const Particle &b) 
  { return dabs(a.Momentum().Eta())>dabs(b.Momentum().Eta()); }
  static bool OrderPP(Particle * const &a,Particle * const &b) 
  { return dabs(a->Momentum().Eta())>dabs(b->Momentum().Eta()); }
  Order_Up_Eta(): Order_Base(OrderV,OrderP,OrderPP) {}
};

DEFINE_ORDER_GETTER(Order_Up_Eta,Order_Up_Eta_Getter,"ETA_UP","order eta ascending")
