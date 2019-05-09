#include "itensor/all.h"
#include "soccer.h"
#include "measure.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    println("//////////////////////////");
    println("Reading input file ......\n");
    //Parse the input file
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");
    //Read in individual parameters from the input file
    //second argument to getXXX methods is a default
    //in case parameter not provided in input file
    double J1 = input.getReal("J1");
    auto readmps = input.getYesNo("readmps",false);
    auto eneropt = input.getYesNo("eneropt",true);
    auto domeas = input.getYesNo("domeas",false);
    auto meas_spincorr = input.getYesNo("meas_spincorr",false);
    auto quiet = input.getYesNo("quiet",true);

    // Read the sweeps parameters
    auto nsweeps = input.getInt("nsweeps");
    auto table = InputGroup(input,"sweeps");

    // number of sites
    auto N = 60;

    // suggested output file name from model parameters
    std::string runlogfile="runlog_";
    runlogfile += "_J1="; runlogfile += std::to_string(J1);
    runlogfile += ".out";

    println("output file name suggested: ", runlogfile);

    //Create the sweeps class & print
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    SpinHalf sites;
    MPS psi;
    MPO H;

    ////// check whether "psi_file" exists
    ////{
    ////    std::ifstream ifile("psi_file");
    ////    if(ifile) {
    ////        readmps=true;
    ////        println("psi_file is found, readmps is set to true"); }
    ////    else{
    ////        readmps=false;
    ////        println("psi_file is NOT found, readmps is set to false"); }
    ////}

    if(readmps) {
        println("\n//////////////////////////////////////////////////");
        println("Reading basis, wavefunction and H from files ......");
        //SpinHalf sites;
        readFromFile("sites_file", sites);
        //MPS psi(sites);
        psi=MPS(sites);
        readFromFile("psi_file", psi);
        //MPO H(sites);
        H=MPO(sites);
        readFromFile("H_file", H);
        //auto psiHpsi = overlap(psi,H,psi);
        ////auto psiHHpsi = overlap(psi,H,H,psi);
        //printfln(" Intial energy information from input file: ");
        //printfln("\n<psi|H|psi> = %.10f", psiHpsi );
        ////printfln("\n<psi|H^2|psi> = %.10f", psiHHpsi );
        ////printfln("\n<psi|H^2|psi> - <psi|H|psi>^2 = %.10f", psiHHpsi-psiHpsi*psiHpsi );
        //println("\nTotal QN of Ground State = ",totalQN(psi));
    }
    else {
        println("\n//////////////////////////////////////////////////");
        println("Building basis, wavefunction and H from scratch ......\n");
        //
        // Initialize the site degrees of freedom.
        //
        //auto sites = SpinHalf(N);
        sites = SpinHalf(N);

        //
        // Use the AutoMPO feature to create the 
        // nearest-neighbor XXZ model with ring exchange.
        //

        auto ampo = AutoMPO(sites);
        auto lattice = soccer();

        println("H is made up of ");
        println("\nBound:\n", lattice);
        println("Total number of nn bound: ", lattice.size());

        // two-body term, nearest neighbor
        for(auto bnd : lattice)
        {
            ampo += 4.0*J1,"Sz",bnd.s1,"Sz",bnd.s2;
        }

        //auto H = MPO(ampo);
        H = MPO(ampo);

        // Set the initial wavefunction matrix product state
        // to be a Neel state.
        //
        auto state = InitState(sites);
        if(J1 > 0 ) {
            for(int i = 1; i <= N; ++i)
                state.set(i,"Up");
        } else {
            for(int i = 1; i <= N; ++i) {
                if(i%2 == 1)
                    state.set(i,"Up");
                else
                    state.set(i,"Dn");
            }
        }

        //auto psi = MPS(state);
        psi = MPS(state);

        //
        // overlap calculates matrix elements of MPO's with respect to MPS's
        // overlap(psi,H,psi) = <psi|H|psi>
        //
        printfln("\nInitial energy = %.5f\n", overlap(psi,H,psi));
    }


    if(eneropt){
        println("\n//////////////////////////////////////////////////////////////////");
        println("Beigin energy optimization, to get ground state wavefunction ......");
        //
        // Begin the DMRG calculation
        //
        auto energy = dmrg(psi,H,sweeps,{"Quiet",quiet,"WriteM",3100});

        // after the MPS converged, write basis, psi, and H to disk
        writeToFile("sites_file", sites);
        writeToFile("psi_file", psi);
        writeToFile("H_file", H);

        //
        // Calculate entanglement entropy
        for ( int b = 1; b <= N-1; ++b ) {
            psi.position(b);
            auto wf = psi.A(b)*psi.A(b+1);
            auto U = psi.A(b);
            ITensor S, V;
            auto spectrum = svd(wf,U,S,V);
            Real SvN = 0.;
            for ( auto p : spectrum.eigs() ) {
                if ( p > 1E-12 ) SvN += -p*log(p);
            }
            printfln("Across bond b=%d, SvN = %.10f", b, SvN);
        }

        //
        // Print the final energy reported by DMRG
        //
        //println("\nTotal QN of Ground State = ",totalQN(psi));
        printfln("\nGround State Energy = %.10f", energy);
        printfln("\n<psi|H|psi> / N = %.10f", energy/N );

        // since calculate overlap(psi,H,H,psi) is extremely memory expensive, we only calculate it when maxm <= max(640, 2560*16/N)
        if( sweeps.maxm( sweeps.nsweep() ) <= std::max(160, 640*16/N) ) { 
            auto psiHHpsi = overlap(psi,H,H,psi);
            printfln("\n<psi|H^2|psi> = %.10f", psiHHpsi );
            printfln("\n<psi|H^2|psi> - <psi|H|psi>^2 = %.10f", psiHHpsi-energy*energy);
            printfln("\n<psi|H^2|psi> / N^2 = %.10f", psiHHpsi/(N*N) );
            printfln("\n( <psi|H^2|psi> - <psi|H|psi>^2 ) / N^2 = %.10f", (psiHHpsi-energy*energy)/(N*N) );
            printfln("\nsqrt( | <psi|H^2|psi> - <psi|H|psi>^2 | ) / N = %.10f", sqrt(abs(psiHHpsi-energy*energy))/N );
        }

    }

    if(domeas && meas_spincorr) {
        println("\n////////////////////////////");
        println("Start to perform measurement of spin correlation\n");
        //
        // Measure Si.Sj of every {i,j}, and total M
        //
        auto totalM = 0.0;
        auto Msquare = 0.0;
        auto Mzsquare = 0.0;
        std::vector<double> SiSj_meas={};
        std::vector<double> SiSjzz_meas={};
        std::vector<double> SiSjpm_meas={};
        std::vector<double> Sz_meas={};
        std::vector<Cplx> Sp_meas={};
        std::vector<Cplx> Sm_meas={};
        for ( int i = 1; i <= N; ++i ) {
            //'gauge' the MPS to site i
            psi.position(i); 
            
            //psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation

            // i == j part

            // magnetization
            auto ket = psi.A(i);
            auto bra = dag(prime(ket,Site));
            auto sz_tmp = ((bra*sites.op("Sz",i)*ket).cplx()).real();
            Sz_meas.emplace_back(sz_tmp);
            totalM +=  sz_tmp;

            auto sp_tmp = (bra*sites.op("S+",i)*ket).cplx();
            Sp_meas.emplace_back(sp_tmp);
            auto sm_tmp = (bra*sites.op("S-",i)*ket).cplx();
            Sm_meas.emplace_back(sm_tmp);

            auto szsz_tmp = 0.0;
            szsz_tmp = (( prime(bra*sites.op("Sz",i),Site)*sites.op("Sz",i)*ket).cplx()).real();
            SiSjzz_meas.emplace_back(szsz_tmp);
            Mzsquare += szsz_tmp;
            auto spsm_tmp = 0.0;
            spsm_tmp = (( prime(bra*sites.op("S+",i),Site)*sites.op("S-",i)*ket).cplx()).real();
            SiSjpm_meas.emplace_back(spsm_tmp);

            auto ss_tmp = szsz_tmp + spsm_tmp;
            Msquare += ss_tmp;
            SiSj_meas.emplace_back(ss_tmp);
            println( i, " ", i, " ", ss_tmp );
            
            if ( i < N ) {
                // i != j part
                //index linking i to i+1:
                auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);
   
                auto op_ip = sites.op("S+",i);
                auto op_im = sites.op("S-",i);
                auto op_iz = sites.op("Sz",i);
                auto Cpm = psi.A(i)*op_ip*dag(prime(psi.A(i),Site,ir));
                auto Cmp = psi.A(i)*op_im*dag(prime(psi.A(i),Site,ir));
                auto Czz = psi.A(i)*op_iz*dag(prime(psi.A(i),Site,ir));
                for(int j = i+1; j <= N; ++j) {
                    Cpm *= psi.A(j);
                    Cmp *= psi.A(j);
                    Czz *= psi.A(j);

                    auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);

                    auto op_jm = sites.op("S-",j);
                    auto ss_tmp = 0.0;
                    ss_tmp += 0.5*(( (Cpm*op_jm)*dag(prime(psi.A(j),jl,Site)) ).cplx()).real();

                    auto op_jp = sites.op("S+",j);
                    ss_tmp += 0.5*(( (Cmp*op_jp)*dag(prime(psi.A(j),jl,Site)) ).cplx()).real();

                    auto spm_tmp = ss_tmp;
                    SiSjpm_meas.emplace_back(ss_tmp);

                    auto op_jz = sites.op("Sz",j);
                    ss_tmp += (( (Czz*op_jz)*dag(prime(psi.A(j),jl,Site)) ).cplx()).real();

                    SiSjzz_meas.emplace_back(ss_tmp-spm_tmp);
                    Mzsquare += (ss_tmp - spm_tmp)*2.0;
                    
                    SiSj_meas.emplace_back(ss_tmp);
                    Msquare += ss_tmp*2.0;
                    println( i, " ", j, " ", ss_tmp ); 

                    if(j < N) {
                        Cpm *= dag(prime(psi.A(j),Link));
                        Cmp *= dag(prime(psi.A(j),Link));
                        Czz *= dag(prime(psi.A(j),Link));
                    }
                }
            }
        }
        printfln("Total M = %.10e", totalM );
        printfln("Msquare = %.10e", Msquare );
        printfln("Mzsquare = %.10e", Mzsquare );

        std::ofstream fSzout("Siz.out",std::ios::out);
        fSzout.precision(16);
        for (std::vector<double>::const_iterator i = Sz_meas.begin(); i != Sz_meas.end(); ++i)
                fSzout << *i << ' ';

        std::ofstream fSpout("Sip.out",std::ios::out);
        fSpout.precision(16);
        for (std::vector<Cplx>::const_iterator i = Sp_meas.begin(); i != Sp_meas.end(); ++i)
                fSpout << *i << ' ';

        std::ofstream fSmout("Sim.out",std::ios::out);
        fSmout.precision(16);
        for (std::vector<Cplx>::const_iterator i = Sm_meas.begin(); i != Sm_meas.end(); ++i)
                fSmout << *i << ' ';

        std::ofstream fSiSjout("SiSj.out",std::ios::out);
        fSiSjout.precision(16);
        for (std::vector<double>::const_iterator i = SiSj_meas.begin(); i != SiSj_meas.end(); ++i)
                fSiSjout << *i << ' ';

        std::ofstream fSiSjzzout("SiSjzz.out",std::ios::out);
        fSiSjzzout.precision(16);
        for (std::vector<double>::const_iterator i = SiSjzz_meas.begin(); i != SiSjzz_meas.end(); ++i)
                fSiSjzzout << *i << ' ';

        std::ofstream fSiSjpmout("SiSjpm.out",std::ios::out);
        fSiSjpmout.precision(16);
        for (std::vector<double>::const_iterator i = SiSjpm_meas.begin(); i != SiSjpm_meas.end(); ++i)
                fSiSjpmout << *i << ' ';
    }

    println( "\nRUNNING FINISHED ^_^ !!! " );
    return 0;
}
