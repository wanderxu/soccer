#include "itensor/all.h"
using namespace itensor;
int main(int argc, char* argv[]) {
    // a Ising model with four spin
    // s0-----s1
    //  |      |
    //  |      |
    // s3-----s2
    auto J1=1.0;

    // std::vector<Index> bonds = {};
    auto sqrt_coshj = std::sqrt( std::cosh( J1 )/2.0 );
    auto sqrt_sinhj = std::sqrt( std::sinh( J1 )/2.0 );
    println( "sqrt_coshj = ", sqrt_coshj );
    println( "sqrt_sinhj = ", sqrt_sinhj );
    Index a("a",2),
          b("b",2),
          c("c",2);
    // build tensors
    ITensor tQ1(a,b),
            tQ2(b,c),
            tQ3(c,a);

    println( "set tensors");

    // set tensors
    tQ1.set( a(1), b(1),  std::exp(-J1) );
    tQ1.set( a(2), b(1),  std::exp( J1) );
    tQ1.set( a(1), b(2),  std::exp( J1) );
    tQ1.set( a(2), b(2),  std::exp(-J1) );

    tQ2.set( b(1), c(1),  std::exp(-J1) );
    tQ2.set( b(2), c(1),  std::exp( J1) );
    tQ2.set( b(1), c(2),  std::exp( J1) );
    tQ2.set( b(2), c(2),  std::exp(-J1) );

    tQ3.set( c(1), a(1),  std::exp(-J1) );
    tQ3.set( c(2), a(1),  std::exp( J1) );
    tQ3.set( c(1), a(2),  std::exp( J1) );
    tQ3.set( c(2), a(2),  std::exp(-J1) );

    println( "contraction");

    // contraction
    //println( (tQ1*tQ2*tQ3) );
    println( (tQ1*tQ2*tQ3).real() );

}
