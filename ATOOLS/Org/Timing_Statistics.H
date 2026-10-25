#ifndef ATOOLS_Org_Timing_Statistics_H
#define ATOOLS_Org_Timing_Statistics_H

#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace ATOOLS {

  /*!
    Hierarchical wall-clock timers for the event generation. The timers form
    a tree whose root is the event generation itself, measured over the same
    period as the "Event N ( X s total )" summary. A Scoped_Timer adds a node
    below the innermost running one, such that nested timers form the tree.
    Timers only record while the event generation is running.
  */
  class Timing_Statistics {
  public:

    typedef std::chrono::steady_clock Clock;

    struct Node {
      std::string m_name;
      double m_time {0.0};
      unsigned long int m_calls {0};
      Node *p_parent {nullptr};
      std::map<std::string,std::unique_ptr<Node> > m_children;
      Node *Child(const std::string &name);
    };

  private:

    static Node s_root;
    static Node *p_current;
    static Clock::time_point s_start;
    static bool s_running;
    static int s_mode;

    static void Print(std::ostream &str,const Node &node,const size_t depth,
                      const size_t width);
    static size_t NameWidth(const Node &node,const size_t depth);

  public:

    // 0: off, 1: up to the processes of the integration results,
    // 2: also each partonic process
    static void SetMode(const int mode) { s_mode=mode; }
    static void Start();
    static void Stop();

    static bool IsRunning() { return s_running; }
    static Node *Current() { return p_current; }
    static void SetCurrent(Node *const node) { p_current=node; }

    // add a duration measured elsewhere to a chain of nodes below the
    // innermost running timer, e.g. Add(dt,"LO",group,partonic)
    template <typename... Names>
    static void Add(const double &seconds,const Names &...names)
    {
      if (!s_running) return;
      Node *node(p_current);
      ((node=node->Child(names), node->m_time+=seconds, ++node->m_calls),
       ...);
    }

    static void PrintStatistics(std::ostream &str);

  };

  class Scoped_Timer {
  private:
    Timing_Statistics::Node *p_node;
    Timing_Statistics::Clock::time_point m_start;
  public:
    explicit Scoped_Timer(const std::string &name);
    ~Scoped_Timer();
    Scoped_Timer(const Scoped_Timer &)=delete;
    Scoped_Timer &operator=(const Scoped_Timer &)=delete;
  };

}

#endif
