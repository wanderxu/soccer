#include "itensor/all.h"
using namespace itensor;
int main(int argc, char* argv[]) {

    double J=1.0;
    auto maxpoints=2000;
    double alpha=-66.0/90.0;
    std::vector<double> logz={};
    std::vector<double> jlist={};
    for(int i = 0; i <= maxpoints; ++i) {
        double J1 = J+i*0.01;
        //auto sqrt_coshj = std::sqrt( std::cosh( J1 )/2.0 );
        //auto sqrt_sinhj = std::sqrt( std::sinh( J1 )/2.0 );
        //println( "sqrt_coshj = ", sqrt_coshj );
        //println( "sqrt_sinhj = ", sqrt_sinhj );
        //auto q11 = sqrt_coshj + sqrt_sinhj*Cplx(0.0,1.0);
        //auto q12 = sqrt_coshj - sqrt_sinhj*Cplx(0.0,1.0);
        auto tmp1 = std::exp((-1.0+alpha)*J1);
        auto tmp2 = std::exp(( 1.0+alpha)*J1);
        auto tmp3 = std::sqrt( tmp2 - tmp1 );
        auto tmp4 = std::sqrt( tmp2 + tmp1 );
        auto q11 = (tmp4 + tmp3*Cplx(0.0,1.0))/2.0;
        auto q12 = (tmp4 - tmp3*Cplx(0.0,1.0))/2.0;
        println( "q11 = ", q11 );
        println( "q12 = ", q12 );
        Index a1("a1",2),
              a2("a2",2),
              a3("a3",2),
              a4("a4",2),
              a5("a5",2),
              b1("b1",2),
              b2("b2",2),
              b3("b3",2),
              b4("b4",2),
              b5("b5",2);
        // build tensors
        ITensor tQ1(a1,b1,b5),
                tQ2(a2,b2,b1),
                tQ3(a3,b3,b2),
                tQ4(a4,b4,b3),
                tQ5(a5,b5,b4);

        println( "set tensors");

        // set tensors
        tQ1.set( a1(1), b1(1), b5(1),  q11*q11*q11+q12*q12*q12 );
        tQ1.set( a1(2), b1(1), b5(1),  q11*q11*q12+q11*q12*q12 );
        tQ1.set( a1(1), b1(2), b5(1),  q11*q11*q12+q11*q12*q12 );
        tQ1.set( a1(2), b1(2), b5(1),  q11*q12*q12+q11*q11*q12 );
        tQ1.set( a1(1), b1(1), b5(2),  q11*q11*q12+q11*q12*q12 );
        tQ1.set( a1(2), b1(1), b5(2),  q11*q12*q12+q11*q11*q12 );
        tQ1.set( a1(1), b1(2), b5(2),  q11*q12*q12+q11*q11*q12 );
        tQ1.set( a1(2), b1(2), b5(2),  q11*q11*q11+q12*q12*q12 );

        tQ2 = reindex(tQ1, a1,a2, b1,b2, b5,b1);
        tQ3 = reindex(tQ1, a1,a3, b1,b3, b5,b2);
        tQ4 = reindex(tQ1, a1,a4, b1,b4, b5,b3);
        tQ5 = reindex(tQ1, a1,a5, b1,b5, b5,b4);

        ITensor fiveT1 = tQ1*tQ2*tQ3*tQ4*tQ5;
        println( fiveT1 );

        Index c1("c1",2),
              c2("c2",2),
              c3("c3",2),
              c4("c4",2),
              c5("c5",2),
              c6("c6",2),
              c7("c7",2),
              c8("c8",2),
              c9("c9",2),
              c10("c10",2);
        
        ITensor fiveT2(a1,b1,c1,c2,b2);
        ITensor fiveT3(a2,b2,c3,c4,b3);
        ITensor fiveT4(a3,b3,c5,c6,b4);
        ITensor fiveT5(a4,b4,c7,c8,b5);
        ITensor fiveT6(a5,b5,c9,c10,b1);

        fiveT2 = reindex(fiveT1,        a2,b1, a3,c1, a4,c2, a5,b2);
        fiveT3 = reindex(fiveT1, a1,a2, a2,b2, a3,c3, a4,c4, a5,b3);
        fiveT4 = reindex(fiveT1, a1,a3, a2,b3, a3,c5, a4,c6, a5,b4);
        fiveT5 = reindex(fiveT1, a1,a4, a2,b4, a3,c7, a4,c8, a5,b5);
        fiveT6 = reindex(fiveT1, a1,a5, a2,b5, a3,c9, a4,c10, a5,b1);

        ITensor halfT = fiveT1*fiveT2*fiveT3*fiveT4*fiveT5*fiveT6;
        println( halfT );

        println( "final contraction:");
        // final contraction
        auto logz_tmp = std::log((halfT*halfT).cplx().real());
        printfln( "%.16e", logz_tmp  );
        jlist.emplace_back( J1 );
        logz.emplace_back( logz_tmp );
    }
    std::ofstream fout("logz.out",std::ios::out);
    fout.precision(16);
    for(int i = 0; i <= maxpoints; ++i)
        fout << jlist[i] << " " << logz[i] << std::endl ;
}
