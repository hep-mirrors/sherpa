// written by Thomas Kunert

#include "prof.hh"
#include <iomanip>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <vector>
#include <algorithm>
//#include <parameter.hh>
//#include <tensor.hh>
//#include <ext/hash_map>
#include <signal.h>
//#include <limits.h>
#include <string>
#include <iostream>

/* Mit -O2 : Zeitaufwand: ca. 450 Takte User + 330 Takte System pro Funktionsaufruf. 
   Davon entfallen lediglich 150 Takte auf diesen Code hier, der Rest ist zum Zeitmessen. */

using namespace std;

static bool do_prof = false;

static int proftime()
{
    //    static int dummy = start_timer();
    //    return tt;
    tms t;
    times( &t );
    return t.tms_utime;

}

/* const char* as keytype works because the origin of these strings is __PRETTY_FUNCTION__ */

typedef map<const char*, call_node> call_map_type;

/**
   information saved for every element of the function call tree
 */
struct call_info {
    unsigned int calls;
    unsigned int time;
    long memory;
    call_info( ) { memory = calls = time = 0; }
    void add_call( int t, long m ) {
        ++calls;
        time += t;
        memory += m;
    }
    friend ostream& operator<<( ostream& ost, const call_info& ci );
};

/**
   function call tree saved as a tree
 */

struct call_node{
    call_info info_;
    call_map_type  called_;
    call_node* find( const char* s ) { return &called_[s];}
    void add_call( int t, long m ) { info_.add_call( t, m ); }
};

static call_node root_call_node;
static call_node* current_call_node = &root_call_node;

class ProfCollect {
public:
    void add( const ProfStatic* ps );
    void print_profile( ostream& ost, int n );
    //    void print_call_tree( ostream& ost, const call_node* cn, int lv, 
    //			  const string& name, const string& treestr, bool last );
    ~ProfCollect();
private:
    vector<const ProfStatic*> vec;
};

static ProfCollect* pcol = 0;

/*    
Parameter<int> Pprofiling("profiling", 0 );
*/
const int Pprofiling() { return 100; };

template<typename T>
struct Validator {
  T * link;
};

struct ProfStatic : public Validator<ProfStatic> {
    ProfStatic( const char* s );
    ~ProfStatic();
    void add_call( int t, long m );
    int time;
    long memory, mem_min, mem_max;
    int n;
    const char* name;
};

Prof::Prof( ProfStatic* ps0 ) : ps( ps0 )
{
    if( do_prof ){
        start = proftime(); 
        start_mem = 0; //query_tensor_memory();
        last_call_node_ = current_call_node;
        current_call_node = current_call_node->find( ps->name );
    }
}

Prof::~Prof()
{ 
    if( do_prof ) { 
        int t = proftime() - start;
        long m = 0; //query_tensor_memory()-start_mem;
        ps->add_call( t, m ); 
        current_call_node->add_call( t, m );
        current_call_node = last_call_node_;
    } 
}

ProfStaticP::ProfStaticP( const char* s )
{
    p = new ProfStatic( s );
}

ProfStaticP::~ProfStaticP()
{
    delete p;
}

ProfStatic::ProfStatic( const char* s ) : name( s )
{
    if( pcol == 0 ) {
        pcol = new ProfCollect;
        if( Pprofiling() ) do_prof=true;
    }
    pcol->add( this );
    mem_max = -2147483646; //numeric_limits<long>::min();
    mem_min = 2147483646; //numeric_limits<long>::max();
    memory = n = time = 0;
}

ProfStatic::~ProfStatic()
{
    delete pcol; pcol=0;
}

void ProfStatic::add_call( int t, long m )
{ 
    ++n;
    time +=t;
    memory +=m; 
    if (m > mem_max) mem_max=m;
    if (m < mem_min) mem_min=m; 
}


void ProfCollect::add( const ProfStatic* ps )
{
    vec.push_back( ps );
}

bool prof_kleiner( const ProfStatic* a, const ProfStatic* b )
{
    return a->time > b->time;
}

class MemPrint
{
public:
    MemPrint( long a ) : a_(a){}
private:
    long a_;
    friend ostream& operator<<( ostream& ost, const MemPrint& m );
};

ostream& operator<<( ostream& ost, const MemPrint& m )
{
    static char units[] = {' ', 'k', 'M', 'G', 'T'};
    static long factor[] = {1000, 10000, 100000, 1000000, 10000000, 100000000,  1000000000};

    int uindex = 0;
    long d = abs(m.a_);
    int sd = ( m.a_ >= 0 ? 1 : -1 );
    int wi0 = ost.width() - 1; // for the unit
    int wi = wi0 - 3;  
    if ( sd < 0 ) wi--; // 1 zusaetzlich abziehen, fuers Vorzeichen
    if ( wi < 0 ) wi = 0; 
    const int fs = sizeof( factor );
    if ( wi >= fs ) wi = fs-1;
    while( d >= factor[wi] && uindex < static_cast<int> ( sizeof( units )) - 1 ){
        d = ( d + 511 ) / 1024;
        ++uindex;
    }
    ost << setw( wi0 ) << sd*d << units[ uindex ];
    return ost;
}

ostream& operator<<( ostream& ost, const call_info& ci )
{
    if( ci.calls > 0 ) {
        ost << setw(7) << ci.calls <<" ";
        ost << setw(5) << MemPrint( ci.memory ) << " ";
        ost << setw(8) << setprecision(2) << ci.time / 100.0;
    }  else  ost << setw( 22 ) << "";
    return ost;
}

typedef call_map_type::const_iterator node_iter;

bool time_comp( const node_iter& a, const node_iter& b )
{
    return a -> second.info_.time > b -> second.info_.time;
} 

void print_call_tree( ostream& ost, 
		      const call_node* cn, 
		      const string& treestr )
{
    vector< node_iter > v;
    for( node_iter p = cn -> called_.begin(); p != cn -> called_.end(); ++p )
        v.push_back( p );

    sort( v.begin(), v.end(), time_comp );

    const string branchstr=( cn->info_.calls ? "+-" : "" );
 
    for( unsigned int i=0; i<v.size(); ++i ){
        ost << v[ i ] -> second.info_ << "  " << treestr << branchstr << v[ i ] -> first << endl;
        string next_str = treestr + ( i + 1 < v.size() ? "| " : "  " );
        if( cn->info_.calls == 0 ) next_str = "";
        print_call_tree( ost, &(v[ i ]->second), next_str );
    }
} 
    
/*
void ProfCollect::print_call_tree( ostream& ost, const call_node* cn, int lv, 
				   const string& name, const string& treestr, bool last )
{
    if ( lv != 0 ) {
	//	ost.setf( ios::left, ios::adjustfield );
	ost << cn->info_ << "  " << treestr;
	if ( lv > 1 ) ost << "+--";
	ost << name << endl;
    }
    map< string, call_node >::const_iterator p;
    for ( p = cn->called_.begin() ; p != cn->called_.end() ; ){
	string ts = treestr;
	string nn = p->first;
        call_node cnn = p->second;
	++p;
	if (lv > 1) {
	    if ( last ) ts += "   "; 
	    else ts += "|  "; 
	}
	print_call_tree( ost, &cnn, lv+1, nn, ts, p == cn->called_.end() );
    }
}
*/

void ProfCollect::print_profile( ostream& ost, int n )
{
    sort( vec.begin(), vec.end(), prof_kleiner );
    ost << "Execution profile:" << endl;
    int m = vec.size();
    if( n < m && n>0 ) m = n;

    // show total memory consumption
    ost.setf(ios::right, ios::adjustfield);
    ost.setf(ios::fixed, ios::floatfield);
    //    ost << " user time consumed: " << infoutime() << endl;
    //    ost << " maximal memory consumption: " << setw(8) << MemPrint(query_max_tensor_memory()) << endl << endl;

    // first a few explaining words on top
    ost.setf(ios::right, ios::adjustfield);
    ost.setf(ios::fixed, ios::floatfield);
    ost << setw(7) << "calls";
    ost << " " << setw(5) << "mem";
    ost << " " << setw(5) << "min";
    ost << " " << setw(5) << "max";
    ost << " " << setw(8) << "time";
    ost.setf(ios::left, ios::adjustfield);
    ost << " " << "function name" << endl;

    // now the profiled data
    for( int i=0; i<m; ++i )  {
        ost.setf(ios::right, ios::adjustfield);
        ost.setf(ios::fixed, ios::floatfield);
        if( vec[i] -> n > 0 ) {
            ost << setw(7) << vec[i]->n;
            ost << " " << setw(5) << MemPrint(vec[i]->memory);
            ost << " " << setw(5) << MemPrint(vec[i]->mem_min);
            ost << " " << setw(5) << MemPrint(vec[i]->mem_max);
            ost << " " << setw(8) << setprecision(2) << vec[i]->time/100.0; // CLK_TCK
        } else ost << setw(34) << "";
        ost.setf(ios::left, ios::adjustfield);
        ost << "  " << vec[i]->name << endl;
    }
    
    // now the call-trees
    ost<<endl<<endl;
    ost.setf(ios::right, ios::adjustfield);
    //    print_call_tree( ost, current_call_node, 0, "", "", 1 );
    ::print_call_tree( ost, &root_call_node, "");
}
    

ProfCollect::~ProfCollect()
{
  //    if( Pprofiling() ) print_profile( cout, Pprofiling() );
}

void print_profile( ostream& ost )
{
    if( pcol ) pcol -> print_profile( ost, 0 );
}
 
void set_prof()
{
    do_prof = true;
}


#pragma instantiate std::vector<const ProfStatic*>
#pragma instantiate std::vector< node_iter >
