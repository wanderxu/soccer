//
// Writen by Xiao Yan Xu (wanderxu@gmail.com)
//
#ifndef __MEASURE_H_
#define __MEASURE_H_

//#include "latticeplaque.h"

namespace itensor {

typedef std::pair<std::string, int> opair;

bool cmp_by_value(const opair& lhs, const opair& rhs) {  
      return lhs.second < rhs.second;  
  }

template <class Tensor>
Real
mtwobody(MPSt<Tensor>& psi, 
     SpinHalf const& sites,
     std::vector<int> const& sites_tmp,
     std::string const& op1_label,
     std::string const& op2_label )
    {
    std::vector<opair> opstr = { std::make_pair(op1_label,sites_tmp[0]),
                                 std::make_pair(op2_label,sites_tmp[1])};
    // sort by site index
    std::sort(opstr.begin(), opstr.end(), cmp_by_value);

    auto opi = sites.op(opstr[0].first,opstr[0].second);
    auto opj = sites.op(opstr[1].first,opstr[1].second);

    psi.position(opstr[0].second);
    ITensor SStmp=psi.A(opstr[0].second);
    if( opstr[1].second != opstr[0].second ) {
        //println("i!=j");
        SStmp *= opi;
        auto ir1 = commonIndex(psi.A(opstr[0].second), psi.A(opstr[0].second+1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[0].second),Site),ir1));
        // propogate until opj
        for ( int i1 = opstr[0].second+1; i1<opstr[1].second; ++i1) { 
            SStmp *= psi.A(i1);
            SStmp *= dag(prime(psi.A(i1),Link));
        }
        SStmp *= psi.A( opstr[1].second );
        SStmp *= opj;
        auto ir2 = commonIndex(psi.A(opstr[1].second), psi.A(opstr[1].second-1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),ir2)); // return i!=j
    }
    else{ // opstr[1].second == opstr[0].second
        //println("i=j");
        SStmp *= opi;
        SStmp = noprime(SStmp,Site)*opj;
        SStmp *= dag(prime(psi.A(opstr[1].second),Site));
    }
    return (SStmp.cplx()).real();
  }

template <class Tensor>
Real
mthreebody(MPSt<Tensor>& psi, 
     SpinHalf const& sites,
     std::vector<int> const& sites_tmp,
     std::string const& op1_label,
     std::string const& op2_label,
     std::string const& op3_label )
    {
    std::vector<opair> opstr = { std::make_pair(op1_label,sites_tmp[0]),
                                 std::make_pair(op2_label,sites_tmp[1]),
                                 std::make_pair(op3_label,sites_tmp[2])};
    // sort by site index
    std::sort(opstr.begin(), opstr.end(), cmp_by_value);

    auto opi = sites.op(opstr[0].first,opstr[0].second);
    auto opj = sites.op(opstr[1].first,opstr[1].second);
    auto opk = sites.op(opstr[2].first,opstr[2].second);

    psi.position(opstr[0].second);
    ITensor SStmp=psi.A(opstr[0].second);
    SStmp *= opi;
    if( opstr[1].second != opstr[0].second ) {
        //println("i!=j");
        auto ir1 = commonIndex(psi.A(opstr[0].second), psi.A(opstr[0].second+1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[0].second),Site),ir1));
        // propogate until opj
        for ( int i1 = opstr[0].second+1; i1<opstr[1].second; ++i1) { 
            SStmp *= psi.A(i1);
            SStmp *= dag(prime(psi.A(i1),Link));
        }
        SStmp *= psi.A(opstr[1].second);
        SStmp *= opj;
    }
    else {
        //println("i=j");
        SStmp = noprime(SStmp,Site)*opj;
    }

    if ( opstr[2].second != opstr[1].second ){
        //println("j!=k");
        if ( opstr[1].second == opstr[0].second ) {
            auto ir2 = commonIndex(psi.A(opstr[1].second), psi.A(opstr[1].second+1),Link);
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),ir2));
        }
        else {
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),Link));
        }
        // propogate until opk
        for ( int i2 = opstr[1].second+1; i2<opstr[2].second; ++i2) { 
            SStmp *= psi.A(i2);
            SStmp *= dag(prime(psi.A(i2),Link));
        }
        SStmp *= psi.A(opstr[2].second);
        SStmp *= opk;
    }
    else{
        //println("j=k");
        SStmp = noprime(SStmp,Site)*opk;
    }

    if ( (opstr[1].second == opstr[0].second) && 
         (opstr[2].second == opstr[1].second) ) {
        SStmp *= dag(prime(psi.A(opstr[2].second),Site)); // return
    }
    else {
        auto ir2 = commonIndex(psi.A(opstr[2].second), psi.A(opstr[2].second-1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),ir2)); // return
    }
    return (SStmp.cplx()).real();
  }

template <class Tensor>
Cplx
mthreebodyC(MPSt<Tensor>& psi, 
     SpinHalf const& sites,
     std::vector<int> const& sites_tmp,
     std::string const& op1_label,
     std::string const& op2_label,
     std::string const& op3_label )
    {
    std::vector<opair> opstr = { std::make_pair(op1_label,sites_tmp[0]),
                                 std::make_pair(op2_label,sites_tmp[1]),
                                 std::make_pair(op3_label,sites_tmp[2])};
    // sort by site index
    std::sort(opstr.begin(), opstr.end(), cmp_by_value);

    auto opi = sites.op(opstr[0].first,opstr[0].second);
    auto opj = sites.op(opstr[1].first,opstr[1].second);
    auto opk = sites.op(opstr[2].first,opstr[2].second);

    psi.position(opstr[0].second);
    ITensor SStmp=psi.A(opstr[0].second);
    SStmp *= opi;
    if( opstr[1].second != opstr[0].second ) {
        //println("i!=j");
        auto ir1 = commonIndex(psi.A(opstr[0].second), psi.A(opstr[0].second+1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[0].second),Site),ir1));
        // propogate until opj
        for ( int i1 = opstr[0].second+1; i1<opstr[1].second; ++i1) { 
            SStmp *= psi.A(i1);
            SStmp *= dag(prime(psi.A(i1),Link));
        }
        SStmp *= psi.A(opstr[1].second);
        SStmp *= opj;
    }
    else {
        //println("i=j");
        SStmp = noprime(SStmp,Site)*opj;
    }

    if ( opstr[2].second != opstr[1].second ){
        //println("j!=k");
        if ( opstr[1].second == opstr[0].second ) {
            auto ir2 = commonIndex(psi.A(opstr[1].second), psi.A(opstr[1].second+1),Link);
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),ir2));
        }
        else {
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),Link));
        }
        // propogate until opk
        for ( int i2 = opstr[1].second+1; i2<opstr[2].second; ++i2) { 
            SStmp *= psi.A(i2);
            SStmp *= dag(prime(psi.A(i2),Link));
        }
        SStmp *= psi.A(opstr[2].second);
        SStmp *= opk;
    }
    else{
        //println("j=k");
        SStmp = noprime(SStmp,Site)*opk;
    }

    if ( (opstr[1].second == opstr[0].second) && 
         (opstr[2].second == opstr[1].second) ) {
        SStmp *= dag(prime(psi.A(opstr[2].second),Site)); // return
    }
    else {
        auto ir2 = commonIndex(psi.A(opstr[2].second), psi.A(opstr[2].second-1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),ir2)); // return
    }
    return SStmp.cplx();
  }

template <class Tensor>
Real
mfourbody(MPSt<Tensor>& psi, 
     SpinHalf const& sites,
     std::vector<int> const& sites_tmp,
     std::string const& op1_label,
     std::string const& op2_label,
     std::string const& op3_label,
     std::string const& op4_label )
    {
    std::vector<opair> opstr = { std::make_pair(op1_label,sites_tmp[0]),
                                 std::make_pair(op2_label,sites_tmp[1]),
                                 std::make_pair(op3_label,sites_tmp[2]),
                                 std::make_pair(op4_label,sites_tmp[3])};
    // sort by site index
    std::sort(opstr.begin(), opstr.end(), cmp_by_value);

    auto opi = sites.op(opstr[0].first,opstr[0].second);
    auto opj = sites.op(opstr[1].first,opstr[1].second);
    auto opk = sites.op(opstr[2].first,opstr[2].second);
    auto opl = sites.op(opstr[3].first,opstr[3].second);

    psi.position(opstr[0].second);
    ITensor SStmp=psi.A(opstr[0].second);
    SStmp *= opi;
    if( opstr[1].second != opstr[0].second ) {
        //println("i!=j");
        auto ir1 = commonIndex(psi.A(opstr[0].second), psi.A(opstr[0].second+1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[0].second),Site),ir1));
        // propogate until opj
        for ( int i1 = opstr[0].second+1; i1<opstr[1].second; ++i1) { 
            SStmp *= psi.A(i1);
            SStmp *= dag(prime(psi.A(i1),Link));
        }
        SStmp *= psi.A(opstr[1].second);
        SStmp *= opj;
    }
    else {
        //println("i=j");
        SStmp = noprime(SStmp,Site)*opj;
    }

    if ( opstr[2].second != opstr[1].second ){
        //println("j!=k");
        if ( opstr[1].second == opstr[0].second ) {
            auto ir2 = commonIndex(psi.A(opstr[1].second), psi.A(opstr[1].second+1),Link);
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),ir2));
        }
        else {
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),Link));
        }
        // propogate until opk
        for ( int i2 = opstr[1].second+1; i2<opstr[2].second; ++i2) { 
            SStmp *= psi.A(i2);
            SStmp *= dag(prime(psi.A(i2),Link));
        }
        SStmp *= psi.A(opstr[2].second);
        SStmp *= opk;
    }
    else{
        //println("j=k");
        SStmp = noprime(SStmp,Site)*opk;
    }

    if ( opstr[3].second != opstr[2].second ){
        //println("k!=l");
        if ( (opstr[1].second == opstr[0].second) && (opstr[2].second == opstr[1].second) ) {
            auto ir3 = commonIndex(psi.A(opstr[2].second), psi.A(opstr[2].second+1),Link);
            SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),ir3));
        }
        else{
            SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),Link));
        }
        // propogate until opk
        for ( int i3 = opstr[2].second+1; i3<opstr[3].second; ++i3) { 
            SStmp *= psi.A(i3);
            SStmp *= dag(prime(psi.A(i3),Link));
        }
        SStmp *= psi.A(opstr[3].second);
        SStmp *= opl;
    }
    else{
        //println("k=l");
        SStmp = noprime(SStmp,Site)*opl;
    }

    if ( (opstr[1].second == opstr[0].second) && 
         (opstr[2].second == opstr[1].second) &&
         (opstr[3].second == opstr[2].second) ) {
        SStmp *= dag(prime(psi.A(opstr[3].second),Site)); // return
    }
    else {
        auto ir3 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second-1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir3)); // return
    }
    return (SStmp.cplx()).real();
  }

struct opairstruct {
    std::vector<opair> op_pair;
    int ind; };

//bool cmp_by_value_of_1st_element_of_struct(const std::vector<opair>& lhs, const std::vector<opair>& rhs) {  
bool cmp_by_value_of_1st_element_of_struct(const opairstruct& lhs, const opairstruct& rhs) {  
      return lhs.op_pair[0].second < rhs.op_pair[0].second;  
}

// a serials of fourbody correlation, (op1,op2) is fixed, (op3,op4) in range op34array
// note the requirement of the input is: op1!=op2!=op3!=op4, and op3,op4 > op1,op2
template <class Tensor>
void
mfourbody_str(MPSt<Tensor>& psi,
     SpinHalf const& sites,
     std::vector<int> const& sites_12,
     std::string const& op1_label,
     std::string const& op2_label,
     std::vector< std::pair<int,int> > const& sites_34_array,
     std::string const& op3_label,
     std::string const& op4_label,
     std::vector<int>& corr_ind,
     std::vector<double>& corr_meas,
     double const& cfac )  // cfac is a constant factor
    {
    // op1, op2 is fixed
    std::vector<opair> opstr_12 = { std::make_pair(op1_label,sites_12[0]),
                                    std::make_pair(op2_label,sites_12[1])};
    // sort by site index
    std::sort(opstr_12.begin(), opstr_12.end(), cmp_by_value);

    auto opi = sites.op(opstr_12[0].first,opstr_12[0].second);
    auto opj = sites.op(opstr_12[1].first,opstr_12[1].second);
    //auto opk = sites.op(opstr[2].first,opstr[2].second);
    //auto opl = sites.op(opstr[3].first,opstr[3].second);

    std::vector< opairstruct > opstr_34_array;
    for ( int kli = 0; kli < sites_34_array.size(); ++kli ){
        std::vector<opair> opstr_34_tmp = { std::make_pair(op3_label, sites_34_array[kli].first),
                                            std::make_pair(op4_label, sites_34_array[kli].second) };
        // sort by site index
        std::sort(opstr_34_tmp.begin(), opstr_34_tmp.end(), cmp_by_value);
        opairstruct opairstruct_tmp;
        opairstruct_tmp.op_pair = opstr_34_tmp;
        opairstruct_tmp.ind = corr_ind[kli];
        opstr_34_array.emplace_back( opairstruct_tmp );
    }

    // sort by site index of first operator
    std::sort( opstr_34_array.begin(), opstr_34_array.end(), cmp_by_value_of_1st_element_of_struct);

    ////// for test
    ////println( "(i,j)=", opstr_12[0].second, opstr_12[1].second );
    ////for ( int kli = 0; kli < sites_34_array.size(); ++kli ){
    ////    println( "(k,l), ind =", opstr_34_array[kli].op_pair[0].second," ",
    ////                             opstr_34_array[kli].op_pair[1].second," ", 
    ////                             opstr_34_array[kli].ind );
    ////}

    psi.position(opstr_12[0].second);
    ITensor SStmp=psi.A(opstr_12[0].second);
    ITensor SStmp2=psi.A(opstr_12[0].second);
    SStmp *= opi;
    auto ir1 = commonIndex(psi.A(opstr_12[0].second), psi.A(opstr_12[0].second+1),Link);
    SStmp *= dag(prime(prime(psi.A(opstr_12[0].second),Site),ir1));
    // propogate until opj
    for ( int i1 = opstr_12[0].second+1; i1<opstr_12[1].second; ++i1) { 
        SStmp *= psi.A(i1);
        SStmp *= dag(prime(psi.A(i1),Link));
    }
    SStmp *= psi.A( opstr_12[1].second );
    SStmp *= opj;
    SStmp *= dag(prime(prime(psi.A(opstr_12[1].second),Site),Link));

    // propogate until first opk
    for ( int i1 = opstr_12[1].second+1; i1<opstr_34_array[0].op_pair[0].second; ++i1){
        SStmp *= psi.A(i1);
        SStmp *= dag(prime(psi.A(i1),Link));
    }

    for( int kli = 0; kli < sites_34_array.size(); ++kli ) {
        auto opk_label = opstr_34_array[kli].op_pair[0].first;
        auto k = opstr_34_array[kli].op_pair[0].second;
        auto opl_label = opstr_34_array[kli].op_pair[1].first;
        auto l = opstr_34_array[kli].op_pair[1].second;
        SStmp2 = SStmp * psi.A(k);
        auto opk = sites.op(opk_label,k);
        SStmp2 *= opk;
        SStmp2 *= dag(prime(prime(psi.A(k),Site),Link));
        // propogate until opl (pair with opk)
        for ( int i2 = k+1; i2<l; ++i2 ){
            SStmp2 *= psi.A(i2);
            SStmp2 *= dag(prime(psi.A(i2),Link));
        }

        SStmp2 *= psi.A( l );
        auto opl = sites.op(opl_label,l);
        SStmp2 *= opl;
        auto ir3 = commonIndex(psi.A(l), psi.A(l-1),Link);
        SStmp2 *= dag(prime(prime(psi.A(l),Site),ir3)); // i!=j!=k!=l
        corr_meas[ opstr_34_array[kli].ind ] += cfac*(SStmp2.cplx()).real(); // add measurement to corr_meas
        if( kli+1 < sites_34_array.size() ) {
            for ( int i3 = k; i3<opstr_34_array[kli+1].op_pair[0].second; ++i3 ){
                SStmp *= psi.A(i3);
                SStmp *= dag(prime(psi.A(i3),Link));
            }
        }
    }
  }

template <class Tensor>
Real
msixbody(MPSt<Tensor>& psi, 
     SpinHalf const& sites,
     std::vector<int> const& sites_tmp,
     std::string const& op1_label,
     std::string const& op2_label,
     std::string const& op3_label,
     std::string const& op4_label,
     std::string const& op5_label,
     std::string const& op6_label )
    {
    std::vector<opair> opstr = { std::make_pair(op1_label,sites_tmp[0]),
                                 std::make_pair(op2_label,sites_tmp[1]),
                                 std::make_pair(op3_label,sites_tmp[2]),
                                 std::make_pair(op4_label,sites_tmp[3]),
                                 std::make_pair(op5_label,sites_tmp[4]),
                                 std::make_pair(op6_label,sites_tmp[5])};
    // sort by site index
    std::sort(opstr.begin(), opstr.end(), cmp_by_value);

    auto opi = sites.op(opstr[0].first,opstr[0].second);
    auto opj = sites.op(opstr[1].first,opstr[1].second);
    auto opk = sites.op(opstr[2].first,opstr[2].second);
    auto opl = sites.op(opstr[3].first,opstr[3].second);
    auto opm = sites.op(opstr[4].first,opstr[4].second);
    auto opn = sites.op(opstr[5].first,opstr[5].second);

    psi.position(opstr[0].second);
    ITensor SStmp=psi.A(opstr[0].second);
    SStmp *= opi;
    if( opstr[1].second != opstr[0].second ) {
        //println("i!=j");
        auto ir1 = commonIndex(psi.A(opstr[0].second), psi.A(opstr[0].second+1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[0].second),Site),ir1));
        // propogate until opj
        for ( int i1 = opstr[0].second+1; i1<opstr[1].second; ++i1) { 
            SStmp *= psi.A(i1);
            SStmp *= dag(prime(psi.A(i1),Link));
        }
        SStmp *= psi.A(opstr[1].second);
        SStmp *= opj;
    }
    else {
        //println("i=j");
        SStmp = noprime(SStmp,Site)*opj;
    }

    if ( opstr[2].second != opstr[1].second ){
        //println("j!=k");
        if ( opstr[1].second == opstr[0].second ) {
            auto ir2 = commonIndex(psi.A(opstr[1].second), psi.A(opstr[1].second+1),Link);
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),ir2));
        }
        else {
            SStmp *= dag(prime(prime(psi.A(opstr[1].second),Site),Link));
        }
        // propogate until opk
        for ( int i2 = opstr[1].second+1; i2<opstr[2].second; ++i2) { 
            SStmp *= psi.A(i2);
            SStmp *= dag(prime(psi.A(i2),Link));
        }
        SStmp *= psi.A(opstr[2].second);
        SStmp *= opk;
    }
    else{
        //println("j=k");
        SStmp = noprime(SStmp,Site)*opk;
    }

    if ( opstr[3].second != opstr[2].second ){
        //println("k!=l");
        if ( (opstr[1].second == opstr[0].second) && (opstr[2].second == opstr[1].second) ) {
            auto ir3 = commonIndex(psi.A(opstr[2].second), psi.A(opstr[2].second+1),Link);
            SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),ir3));
        }
        else{
            SStmp *= dag(prime(prime(psi.A(opstr[2].second),Site),Link));
        }
        // propogate until opk
        for ( int i3 = opstr[2].second+1; i3<opstr[3].second; ++i3) { 
            SStmp *= psi.A(i3);
            SStmp *= dag(prime(psi.A(i3),Link));
        }
        SStmp *= psi.A(opstr[3].second);
        SStmp *= opl;
    }
    else{
        //println("k=l");
        SStmp = noprime(SStmp,Site)*opl;
    }

    if ( opstr[4].second != opstr[3].second ){
        //println("l!=m");
        if ( (opstr[1].second == opstr[0].second) && 
             (opstr[2].second == opstr[1].second) &&
             (opstr[3].second == opstr[2].second) ) {
            auto ir4 = commonIndex(psi.A(opstr[3].second), psi.A(opstr[3].second+1),Link);
            SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),ir4));
        }
        else{
            SStmp *= dag(prime(prime(psi.A(opstr[3].second),Site),Link));
        }
        // propogate until opk
        for ( int i4 = opstr[3].second+1; i4<opstr[4].second; ++i4) { 
            SStmp *= psi.A(i4);
            SStmp *= dag(prime(psi.A(i4),Link));
        }
        SStmp *= psi.A(opstr[4].second);
        SStmp *= opm;
    }
    else{
        //println("l=m");
        SStmp = noprime(SStmp,Site)*opm;
    }

    if ( opstr[5].second != opstr[4].second ){
        //println("m!=n");
        if ( (opstr[1].second == opstr[0].second) && 
             (opstr[2].second == opstr[1].second) &&
             (opstr[3].second == opstr[2].second) &&
             (opstr[4].second == opstr[3].second) ) {
            auto ir5 = commonIndex(psi.A(opstr[4].second), psi.A(opstr[4].second+1),Link);
            SStmp *= dag(prime(prime(psi.A(opstr[4].second),Site),ir5));
        }
        else {
            SStmp *= dag(prime(prime(psi.A(opstr[4].second),Site),Link));
        }
        // propogate until opk
        for ( int i5 = opstr[4].second+1; i5<opstr[5].second; ++i5) { 
            SStmp *= psi.A(i5);
            SStmp *= dag(prime(psi.A(i5),Link));
        }
        SStmp *= psi.A(opstr[5].second);
        SStmp *= opn;
    }
    else{
        //println("m=n");
        SStmp = noprime(SStmp,Site)*opn;
    }
    if ( (opstr[1].second == opstr[0].second) && 
         (opstr[2].second == opstr[1].second) &&
         (opstr[3].second == opstr[2].second) &&
         (opstr[4].second == opstr[3].second) &&
         (opstr[5].second == opstr[4].second) ) {
        SStmp *= dag(prime(psi.A(opstr[5].second),Site)); // return
    }
    else {
        auto ir5 = commonIndex(psi.A(opstr[5].second), psi.A(opstr[5].second-1),Link);
        SStmp *= dag(prime(prime(psi.A(opstr[5].second),Site),ir5)); // return
    }
    return (SStmp.cplx()).real();
  }

// a serials of sixbody correlation, (op1,op2,op3) is fixed, (op4,op5,op6) in range op456array
// note the requirement of the input is: op1!=op2!=op3!=op4=!op5=!op6, and op4,op5,op6 > op1,op2,op3
template <class Tensor>
void
msixbody_str(MPSt<Tensor>& psi,
     SpinHalf const& sites,
     std::vector<int> const& sites_123,
     std::string const& op1_label,
     std::string const& op2_label,
     std::string const& op3_label,
     std::vector< std::vector<int> > const& sites_456_array,
     std::string const& op4_label,
     std::string const& op5_label,
     std::string const& op6_label,
     std::vector<int>& corr_ind,
     std::vector<double>& corr_meas,
     double const& cfac )  // cfac is a constant factor
    {
    // op1, op2 is fixed
    std::vector<opair> opstr_123 = { std::make_pair(op1_label,sites_123[0]),
                                     std::make_pair(op2_label,sites_123[1]),
                                     std::make_pair(op3_label,sites_123[2])};
    // sort by site index
    std::sort(opstr_123.begin(), opstr_123.end(), cmp_by_value);

    auto opi = sites.op(opstr_123[0].first,opstr_123[0].second);
    auto opj = sites.op(opstr_123[1].first,opstr_123[1].second);
    auto opk = sites.op(opstr_123[2].first,opstr_123[2].second);
    //auto opk = sites.op(opstr[2].first,opstr[2].second);
    //auto opl = sites.op(opstr[3].first,opstr[3].second);

    std::vector< opairstruct > opstr_456_array;
    for ( int kli = 0; kli < sites_456_array.size(); ++kli ){
        std::vector<opair> opstr_456_tmp = { std::make_pair(op4_label, sites_456_array[kli][0]),
                                             std::make_pair(op5_label, sites_456_array[kli][1]),
                                             std::make_pair(op6_label, sites_456_array[kli][2]) };
        // sort by site index
        std::sort(opstr_456_tmp.begin(), opstr_456_tmp.end(), cmp_by_value);
        opairstruct opairstruct_tmp;
        opairstruct_tmp.op_pair = opstr_456_tmp;
        opairstruct_tmp.ind = corr_ind[kli];
        opstr_456_array.emplace_back( opairstruct_tmp );
    }

    // sort by site index of first operator
    std::sort( opstr_456_array.begin(), opstr_456_array.end(), cmp_by_value_of_1st_element_of_struct);

    ////// for test
    ////println( "(i,j,k)=", opstr_123[0].second, opstr_123[1].second, opstr_123[2].second );
    ////for ( int kli = 0; kli < sites_456_array.size(); ++kli ){
    ////    println( "(l,m,n), ind =", opstr_456_array[kli].op_pair[0].second," ",
    ////                               opstr_456_array[kli].op_pair[1].second," ", 
    ////                               opstr_456_array[kli].op_pair[2].second," ", 
    ////                               opstr_456_array[kli].ind );
    ////}

    psi.position(opstr_123[0].second);
    ITensor SStmp=psi.A(opstr_123[0].second);
    ITensor SStmp2=psi.A(opstr_123[0].second);
    SStmp *= opi;
    auto ir1 = commonIndex(psi.A(opstr_123[0].second), psi.A(opstr_123[0].second+1),Link);
    SStmp *= dag(prime(prime(psi.A(opstr_123[0].second),Site),ir1));
    // propogate until opj
    for ( int i1 = opstr_123[0].second+1; i1<opstr_123[1].second; ++i1) { 
        SStmp *= psi.A(i1);
        SStmp *= dag(prime(psi.A(i1),Link));
    }
    SStmp *= psi.A( opstr_123[1].second );
    SStmp *= opj;
    SStmp *= dag(prime(prime(psi.A(opstr_123[1].second),Site),Link));

    // propogate until opk
    for ( int i1 = opstr_123[1].second+1; i1<opstr_123[2].second; ++i1) { 
        SStmp *= psi.A(i1);
        SStmp *= dag(prime(psi.A(i1),Link));
    }
    SStmp *= psi.A( opstr_123[2].second );
    SStmp *= opk;
    SStmp *= dag(prime(prime(psi.A(opstr_123[2].second),Site),Link));

    // propogate until first opl
    for ( int i1 = opstr_123[2].second+1; i1<opstr_456_array[0].op_pair[0].second; ++i1){
        SStmp *= psi.A(i1);
        SStmp *= dag(prime(psi.A(i1),Link));
    }

    for( int kli = 0; kli < sites_456_array.size(); ++kli ) {
        auto opl_label = opstr_456_array[kli].op_pair[0].first;
        auto l = opstr_456_array[kli].op_pair[0].second;
        auto opm_label = opstr_456_array[kli].op_pair[1].first;
        auto m = opstr_456_array[kli].op_pair[1].second;
        auto opn_label = opstr_456_array[kli].op_pair[2].first;
        auto n = opstr_456_array[kli].op_pair[2].second;

        // continue things on SStmp
        SStmp2 = SStmp * psi.A(l);
        auto opl = sites.op(opl_label,l);
        SStmp2 *= opl;
        SStmp2 *= dag(prime(prime(psi.A(l),Site),Link));
        // propogate until opm (l->m)
        for ( int i2 = l+1; i2<m; ++i2 ){
            SStmp2 *= psi.A(i2);
            SStmp2 *= dag(prime(psi.A(i2),Link));
        }

        SStmp2 *=  psi.A(m);
        auto opm = sites.op(opm_label,m);
        SStmp2 *= opm;
        SStmp2 *= dag(prime(prime(psi.A(m),Site),Link));
        // propogate until opn (m->n)
        for ( int i2 = m+1; i2<n; ++i2 ){
            SStmp2 *= psi.A(i2);
            SStmp2 *= dag(prime(psi.A(i2),Link));
        }

        SStmp2 *= psi.A( n );
        auto opn = sites.op(opn_label,n);
        SStmp2 *= opn;
        auto ir3 = commonIndex(psi.A(n), psi.A(n-1),Link);
        SStmp2 *= dag(prime(prime(psi.A(n),Site),ir3)); // i!=j!=k!=l!=m!=n
        corr_meas[ opstr_456_array[kli].ind ] += cfac*(SStmp2.cplx()).real(); // add measurement to corr_meas

        if( kli+1 < sites_456_array.size() ) {
            for ( int i3 = l; i3<opstr_456_array[kli+1].op_pair[0].second; ++i3 ){
                SStmp *= psi.A(i3);
                SStmp *= dag(prime(psi.A(i3),Link));
            }
        }
    }
  }

} //namespace itensor

#endif
