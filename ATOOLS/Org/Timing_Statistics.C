#include "ATOOLS/Org/Timing_Statistics.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace ATOOLS;

Timing_Statistics::Node Timing_Statistics::s_root {"Event generation"};
Timing_Statistics::Node *Timing_Statistics::p_current(&Timing_Statistics::s_root);
Timing_Statistics::Clock::time_point Timing_Statistics::s_start;
bool Timing_Statistics::s_running(false);
int Timing_Statistics::s_mode(0);

namespace {
  const std::string s_other("(other)");
  const size_t s_indent(2);
}

Timing_Statistics::Node *Timing_Statistics::Node::Child(const std::string &name)
{
  auto it(m_children.find(name));
  if (it!=m_children.end()) return it->second.get();
  Node *child(new Node{name});
  child->p_parent=this;
  m_children.emplace(name,std::unique_ptr<Node>(child));
  return child;
}

void Timing_Statistics::Start()
{
  if (s_mode==0 || s_running) return;
  s_start=Clock::now();
  s_running=true;
}

void Timing_Statistics::Stop()
{
  if (!s_running) return;
  s_root.m_time+=std::chrono::duration<double>(Clock::now()-s_start).count();
  s_running=false;
}

Scoped_Timer::Scoped_Timer(const std::string &name):
  p_node(nullptr)
{
  if (!Timing_Statistics::IsRunning()) return;
  p_node=Timing_Statistics::Current()->Child(name);
  Timing_Statistics::SetCurrent(p_node);
  m_start=Timing_Statistics::Clock::now();
}

Scoped_Timer::~Scoped_Timer()
{
  if (p_node==nullptr) return;
  p_node->m_time+=std::chrono::duration<double>
    (Timing_Statistics::Clock::now()-m_start).count();
  ++p_node->m_calls;
  Timing_Statistics::SetCurrent(p_node->p_parent);
}

size_t Timing_Statistics::NameWidth(const Node &node,const size_t depth)
{
  size_t width(s_indent*depth+node.m_name.size());
  if (node.m_children.empty() || (s_mode<2 && depth>=3)) return width;
  width=std::max(width,s_indent*(depth+1)+s_other.size());
  for (const auto &child : node.m_children)
    width=std::max(width,NameWidth(*child.second,depth+1));
  return width;
}

void Timing_Statistics::Print(std::ostream &str,const Node &node,
                              const size_t depth,const size_t width)
{
  const double total(s_root.m_time);
  auto line = [&str,width,total](const std::string &name,const double &time,
                                 const size_t depth,const long int calls) {
    const std::ios_base::fmtflags flags(str.flags());
    const std::streamsize precision(str.precision());
    str<<"  "<<std::string(s_indent*depth,' ')<<std::left
       <<std::setw(int(width-s_indent*depth))<<name<<std::right<<std::fixed
       <<std::setprecision(2)<<std::setw(12)<<time<<" s"
       <<std::setprecision(1)<<std::setw(8)
       <<(total>0.0?100.0*time/total:0.0)<<" %";
    if (calls>=0) str<<std::setw(14)<<calls<<" calls";
    str<<"\n";
    str.flags(flags);
    str.precision(precision);
  };
  line(node.m_name,node.m_time,depth,depth==0?-1:long(node.m_calls));
  // the partonic processes are only shown in mode 2
  if (node.m_children.empty() || (s_mode<2 && depth>=3)) return;
  std::vector<const Node*> children;
  double sum(0.0);
  for (const auto &child : node.m_children) {
    children.push_back(child.second.get());
    sum+=child.second->m_time;
  }
  std::sort(children.begin(),children.end(),
            [](const Node *a,const Node *b) { return a->m_time>b->m_time; });
  for (const Node *child : children) Print(str,*child,depth+1,width);
  // time not covered by the timers below, unless they cover all of it, as
  // for the chains filled by Add
  if (std::abs(node.m_time-sum)>=0.005) line(s_other,node.m_time-sum,depth+1,-1);
}

void Timing_Statistics::PrintStatistics(std::ostream &str)
{
  if (s_mode==0 || s_root.m_time==0.0) return;
  str<<"Event generation timing (wall clock):\n";
  Print(str,s_root,0,NameWidth(s_root,0));
}
