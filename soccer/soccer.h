//
// Writen by Xiao Yan Xu (wanderxu@gmail.com)
//
#ifndef __LATTICE_TRIANGULAR_MORE_H_
#define __LATTICE_TRIANGULAR_MORE_H_

#include "latticeplaque.h"

namespace itensor {
struct LatticeBondv2;

using LatticeGraphv2 = std::vector<LatticeBondv2>;

struct LatticeBondv2
    {
    int s1 = 0,
        s2 = 0;
    bool isbd = false;
    std::string type;
    Real x1 = NAN,
         y1 = NAN,
         x2 = NAN,
         y2 = NAN;

    LatticeBondv2() { }

    LatticeBondv2(int s1_, int s2_, bool isbd_)
      : s1{s1_}, 
        s2{s2_},
        isbd{isbd_}
        { }

    LatticeBondv2(int s1_, int s2_, bool isbd_,
                Real x1_,  Real y1_,
                Real x2_,  Real y2_)
      : s1{s1_}, 
        s2{s2_},
        isbd{isbd_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_}
        { }

    LatticeBondv2(int s1_, int s2_, bool isbd_, std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        isbd{isbd_},
        type{type_} 
        { }

    LatticeBondv2(int s1_, int s2_, bool isbd_,
                Real x1_,  Real y1_,
                Real x2_,  Real y2_,
                std::string type_)
      : s1{s1_}, 
        s2{s2_}, 
        isbd{isbd_},
        type{type_},
        x1{x1_},
        y1{y1_},
        x2{x2_},
        y2{y2_}
        { }
    };

inline std::ostream& 
operator<<(std::ostream & s, LatticeBondv2 const& b) 
    { 
    //s << format("(%*d,%*d",3,b.s1,3,b.s2);
    s << format("(%d,%d",b.s1,b.s2);
    if(b.type.size()!=0) s << "," << b.type;
    s << ")";
    if(!std::isnan(b.x1) && !std::isnan(b.y1))
        {
        s << format("[%s,%s",b.x1,b.y1);
        if(!std::isnan(b.x2) && !std::isnan(b.y2))
            {
            s << format(";%s,%s]",b.x2,b.y2);
            }
        else
            {
            s << "]";
            }
        }
    return s;
    }

inline std::ostream& 
operator<<(std::ostream& s, LatticeGraphv2 const& G) 
    { 
    for(auto& b : G)
        {
        s << b << "\n";
        }
    return s;
    }

LatticeGraphv2 inline
soccer( Args const& args = Args::global()) {
    auto N = 60;
    auto Nbond = 90;
    LatticeGraphv2 latt; 
    latt.reserve(Nbond);

    latt.emplace_back(1,2,false,0,0,0,0);
    latt.emplace_back(1,5,false,0,0,0,0);
    latt.emplace_back(1,9,false,0,0,0,0);
    latt.emplace_back(2,3,false,0,0,0,0);
    latt.emplace_back(2,12,false,0,0,0,0);
    latt.emplace_back(3,4,false,0,0,0,0);
    latt.emplace_back(3,15,false,0,0,0,0);
    latt.emplace_back(4,5,false,0,0,0,0);
    latt.emplace_back(4,18,false,0,0,0,0);
    latt.emplace_back(5,6,false,0,0,0,0);
    latt.emplace_back(6,7,false,0,0,0,0);
    latt.emplace_back(6,20,false,0,0,0,0);
    latt.emplace_back(7,8,false,0,0,0,0);
    latt.emplace_back(7,22,false,0,0,0,0);
    latt.emplace_back(8,9,false,0,0,0,0);
    latt.emplace_back(8,25,false,0,0,0,0);
    latt.emplace_back(9,10,false,0,0,0,0);
    latt.emplace_back(10,11,false,0,0,0,0);
    latt.emplace_back(10,26,false,0,0,0,0);
    latt.emplace_back(11,12,false,0,0,0,0);
    latt.emplace_back(11,29,false,0,0,0,0);
    latt.emplace_back(12,13,false,0,0,0,0);
    latt.emplace_back(13,14,false,0,0,0,0);
    latt.emplace_back(13,30,false,0,0,0,0);
    latt.emplace_back(14,15,false,0,0,0,0);
    latt.emplace_back(14,33,false,0,0,0,0);
    latt.emplace_back(15,16,false,0,0,0,0);
    latt.emplace_back(16,17,false,0,0,0,0);
    latt.emplace_back(16,34,false,0,0,0,0);
    latt.emplace_back(17,18,false,0,0,0,0);
    latt.emplace_back(17,37,false,0,0,0,0);
    latt.emplace_back(18,19,false,0,0,0,0);
    latt.emplace_back(19,20,false,0,0,0,0);
    latt.emplace_back(19,38,false,0,0,0,0);
    latt.emplace_back(20,21,false,0,0,0,0);
    latt.emplace_back(21,22,false,0,0,0,0);
    latt.emplace_back(21,40,false,0,0,0,0);
    latt.emplace_back(22,23,false,0,0,0,0);
    latt.emplace_back(23,24,false,0,0,0,0);
    latt.emplace_back(23,42,false,0,0,0,0);
    latt.emplace_back(24,25,false,0,0,0,0);
    latt.emplace_back(24,44,false,0,0,0,0);
    latt.emplace_back(25,26,false,0,0,0,0);
    latt.emplace_back(26,27,false,0,0,0,0);
    latt.emplace_back(27,28,false,0,0,0,0);
    latt.emplace_back(27,45,false,0,0,0,0);
    latt.emplace_back(28,29,false,0,0,0,0);
    latt.emplace_back(28,47,false,0,0,0,0);
    latt.emplace_back(29,30,false,0,0,0,0);
    latt.emplace_back(30,31,false,0,0,0,0);
    latt.emplace_back(31,32,false,0,0,0,0);
    latt.emplace_back(31,48,false,0,0,0,0);
    latt.emplace_back(32,33,false,0,0,0,0);
    latt.emplace_back(32,50,false,0,0,0,0);
    latt.emplace_back(33,34,false,0,0,0,0);
    latt.emplace_back(34,35,false,0,0,0,0);
    latt.emplace_back(35,36,false,0,0,0,0);
    latt.emplace_back(35,51,false,0,0,0,0);
    latt.emplace_back(36,37,false,0,0,0,0);
    latt.emplace_back(36,53,false,0,0,0,0);
    latt.emplace_back(37,38,false,0,0,0,0);
    latt.emplace_back(38,39,false,0,0,0,0);
    latt.emplace_back(39,40,false,0,0,0,0);
    latt.emplace_back(39,54,false,0,0,0,0);
    latt.emplace_back(40,41,false,0,0,0,0);
    latt.emplace_back(41,42,false,0,0,0,0);
    latt.emplace_back(41,55,false,0,0,0,0);
    latt.emplace_back(42,43,false,0,0,0,0);
    latt.emplace_back(43,44,false,0,0,0,0);
    latt.emplace_back(43,57,false,0,0,0,0);
    latt.emplace_back(44,45,false,0,0,0,0);
    latt.emplace_back(45,46,false,0,0,0,0);
    latt.emplace_back(46,47,false,0,0,0,0);
    latt.emplace_back(46,58,false,0,0,0,0);
    latt.emplace_back(47,48,false,0,0,0,0);
    latt.emplace_back(48,49,false,0,0,0,0);
    latt.emplace_back(49,50,false,0,0,0,0);
    latt.emplace_back(49,59,false,0,0,0,0);
    latt.emplace_back(50,51,false,0,0,0,0);
    latt.emplace_back(51,52,false,0,0,0,0);
    latt.emplace_back(52,53,false,0,0,0,0);
    latt.emplace_back(52,60,false,0,0,0,0);
    latt.emplace_back(53,54,false,0,0,0,0);
    latt.emplace_back(54,55,false,0,0,0,0);
    latt.emplace_back(55,56,false,0,0,0,0);
    latt.emplace_back(56,57,false,0,0,0,0);
    latt.emplace_back(56,60,false,0,0,0,0);
    latt.emplace_back(57,58,false,0,0,0,0);
    latt.emplace_back(58,59,false,0,0,0,0);
    latt.emplace_back(59,60,false,0,0,0,0);

    if(int(latt.size()) != Nbond) Error("Wrong number of bonds");

    return latt;
    }

} //namespace itensor

#endif
