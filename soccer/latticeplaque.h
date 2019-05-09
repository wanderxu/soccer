//
// Writen by Xiao Yan Xu (wanderxu@gmail.com)
//
#ifndef __LATTICEPLAQUE_H__
#define __LATTICEPLAQUE_H__

#include <vector>
#include "itensor/global.h"

namespace itensor {

struct Lattice3Plaque;
struct Lattice4Plaque;
struct LatticeNeighbor;

using Lattice3PlaqueGraph = std::vector<Lattice3Plaque>;
using Lattice4PlaqueGraph = std::vector<Lattice4Plaque>;
using LatticeNeighborGraph = std::vector<LatticeNeighbor>;

struct Lattice3Plaque
    {
    int s1 = 0,
        s2 = 0,
        s3 = 0;
    std::string type;
    Real x1 = NAN,
         y1 = NAN,
         x2 = NAN,
         y2 = NAN,
         x3 = NAN,
         y3 = NAN;

    Lattice3Plaque() { }

    Lattice3Plaque(int s1_, int s2_, int s3_)
      : s1{s1_}, 
        s2{s2_}, 
        s3{s3_}
        { }

    Lattice3Plaque(int s1_, int s2_, int s3_,
                Real x1_,  Real y1_,
                Real x2_,  Real y2_,
                Real x3_,  Real y3_)
      : s1{s1_}, 
        s2{s2_},
        s3{s3_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_},
        x3{x3_},
        y3{y3_}
        { }

    Lattice3Plaque(int s1_, int s2_, int s3_, std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        s3{s3_}, 
        type{type_} 
        { }

    Lattice3Plaque(int s1_, int s2_, int s3_,
                Real x1_,  Real y1_,
                Real x2_,  Real y2_,
                Real x3_,  Real y3_,
                std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        s3{s3_}, 
        type{type_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_},
        x3{x3_},
        y3{y3_}
        { }
    };

inline std::ostream& 
operator<<(std::ostream & s, Lattice3Plaque const& b) 
    { 
    //s << format("(%*d,%*d",3,b.s1,3,b.s2);
    s << format("(%d,%d,%d",b.s1,b.s2,b.s3);
    if(b.type.size()!=0) s << "," << b.type;
    s << ")";
    if(!std::isnan(b.x1) && !std::isnan(b.y1))
        {
        s << format("[%s,%s",b.x1,b.y1);
        if(!std::isnan(b.x2) && !std::isnan(b.y2))
            {
            s << format(";%s,%s",b.x2,b.y2);
            if(!std::isnan(b.x3) && !std::isnan(b.y3))
                {
                s << format(";%s,%s]",b.x3,b.y3);
                }
            else
                {
                s << "]";
                }
            }
        else
            {
            s << "]";
            }
        }
    return s;
    }

inline std::ostream& 
operator<<(std::ostream& s, Lattice3PlaqueGraph const& G) 
    { 
    for(auto& b : G)
        {
        s << b << "\n";
        }
    return s;
    }

struct Lattice4Plaque
    {
    int s1 = 0,
        s2 = 0,
        s3 = 0,
        s4 = 0;
    std::string type;
    Real x1 = NAN,
         y1 = NAN,
         x2 = NAN,
         y2 = NAN,
         x3 = NAN,
         y3 = NAN,
         x4 = NAN,
         y4 = NAN;

    Lattice4Plaque() { }

    Lattice4Plaque(int s1_, int s2_, int s3_, int s4_)
      : s1{s1_}, 
        s2{s2_}, 
        s3{s3_}, 
        s4{s4_} 
        { }

    Lattice4Plaque(int s1_, int s2_, int s3_, int s4_,
                Real x1_,  Real y1_,
                Real x2_,  Real y2_,
                Real x3_,  Real y3_,
                Real x4_,  Real y4_)
      : s1{s1_}, 
        s2{s2_},
        s3{s3_},
        s4{s4_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_},
        x3{x3_},
        y3{y3_},
        x4{x4_},
        y4{y4_}
        { }

    Lattice4Plaque(int s1_, int s2_, int s3_, int s4_, std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        s3{s3_}, 
        s4{s4_}, 
        type{type_} 
        { }

    Lattice4Plaque(int s1_, int s2_, int s3_, int s4_,
                Real x1_,  Real y1_,
                Real x2_,  Real y2_,
                Real x3_,  Real y3_,
                Real x4_,  Real y4_,
                std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        s3{s3_}, 
        s4{s4_}, 
        type{type_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_},
        x3{x3_},
        y3{y3_},
        x4{x4_},
        y4{y4_}
        { }
    };

inline std::ostream& 
operator<<(std::ostream & s, Lattice4Plaque const& b) 
    { 
    //s << format("(%*d,%*d",3,b.s1,3,b.s2);
    s << format("(%d,%d,%d,%d",b.s1,b.s2,b.s3,b.s4);
    if(b.type.size()!=0) s << "," << b.type;
    s << ")";
    if(!std::isnan(b.x1) && !std::isnan(b.y1))
        {
        s << format("[%s,%s",b.x1,b.y1);
        if(!std::isnan(b.x2) && !std::isnan(b.y2))
            {
            s << format(";%s,%s",b.x2,b.y2);
            if(!std::isnan(b.x3) && !std::isnan(b.y3))
                {
                s << format(";%s,%s",b.x3,b.y3);
                if(!std::isnan(b.x4) && !std::isnan(b.y4))
                    {
                    s << format(";%s,%s]",b.x4,b.y4);
                    }
                else
                    {
                    s << "]";
                    }
                }
            else
                {
                s << "]";
                }
            }
        else
            {
            s << "]";
            }
        }
    return s;
    }

inline std::ostream& 
operator<<(std::ostream& s, Lattice4PlaqueGraph const& G) 
    { 
    for(auto& b : G)
        {
        s << b << "\n";
        }
    return s;
    }

struct LatticeNeighbor {
    int s0 = 0;
    std::vector<int> snn;
    LatticeNeighbor() { }
    LatticeNeighbor(int s0_, std::vector<int> snn_) : s0{s0_}, snn{snn_} { }
};

inline std::ostream& 
operator<<(std::ostream & s, LatticeNeighbor const& b) { 
    //s << format("(%*d,%*d",3,b.s1,3,b.s2);
    s << format("(%d : ", b.s0);
    for (std::vector<int>::const_iterator i = b.snn.begin(), j = --b.snn.end(); i != j; ++i)
    s << format("%d, ",*i);
    s << format("%d)", b.snn.back());
    return s;
}

inline std::ostream& 
operator<<(std::ostream& s, LatticeNeighborGraph const& G) { 
    for(auto& b : G) {
        s << b << "\n";
    }
    return s;
}

} //namespace itensor

#endif
