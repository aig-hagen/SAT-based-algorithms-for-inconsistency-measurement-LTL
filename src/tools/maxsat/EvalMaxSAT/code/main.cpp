#include <iostream>
#include <cassert>
#include <csignal>
#include <zlib.h>

#include "EvalMaxSAT.h"
#include "lib/CLI11.hpp"


using namespace MaLib;

EvalMaxSAT* monMaxSat = nullptr;

void signalHandler( int signum ) {
    std::cout << "c Interrupt signal (" << signum << ") received"<< std::endl;
    std::cout << "c o >=" << monMaxSat->getCost() << std::endl;
    std::cout << "s UNKNOWN" << std::endl;

   delete monMaxSat;

   exit(signum);
}



int main(int argc, char *argv[])
{
    Chrono chrono("c Total time");
    signal(SIGINT, signalHandler);
    signal(SIGTERM, signalHandler);

    /////// PARSE ARG //////////////////////
    CLI::App app{"EvalMaxSAT Solver"};

    std::string file;
    app.add_option("file", file, "File with the formula to be solved (wcnf format)")->check(CLI::ExistingFile)->required();

    unsigned int paralleleThread=0;
    app.add_option("-p", paralleleThread, toString("Number of minimization threads (default = ",paralleleThread,")"));

    unsigned int timeOutFastMinimize=60;
    app.add_option("--timeout_fast", timeOutFastMinimize, toString("Timeout in second for fast minimize (default = ",timeOutFastMinimize,")"));

    unsigned int coefMinimizeTime=2;
    app.add_option("--coef_minimize", coefMinimizeTime, toString("Multiplying coefficient of the time spent to minimize cores (default = ",coefMinimizeTime,")"));

    bool oldOutputFormat = false;
    app.add_flag("--old", oldOutputFormat, "Use old output format.");

    CLI11_PARSE(app, argc, argv);
    ////////////////////////////////////////

    if(paralleleThread != 0) {
        std::cerr << "Multi thread not yet compatible with CaDiCal" << std::endl;
        return -1;
    }
    auto monMaxSat = new EvalMaxSAT(paralleleThread);
    monMaxSat->setTimeOutFast(timeOutFastMinimize);
    monMaxSat->setCoefMinimize(coefMinimizeTime);

    if(!monMaxSat->parse(file)) {
        std::cerr << "Error parsing : " << file << std::endl;
        return -1;
    }

    if(!monMaxSat->solve()) {
        std::cout << "s UNSATISFIABLE" << std::endl;
        return 0;
    }

    if(oldOutputFormat) {
        ////// PRINT SOLUTION OLD FORMAT //////////////////
        std::cout << "s OPTIMUM FOUND" << std::endl;
        std::cout << "o " << monMaxSat->getCost() << std::endl;
        std::cout << "v";
        for(unsigned int i=1; i<=monMaxSat->nInputVars; i++) {
            if(monMaxSat->getValue(i))
                std::cout << " " << i;
            else
                std::cout << " -" << i;
        }
        std::cout << std::endl;
        ///////////////////////////////////////
    } else {
        ////// PRINT SOLUTION NEW FORMAT //////////////////
        std::cout << "s OPTIMUM FOUND" << std::endl;
        std::cout << "o " << monMaxSat->getCost() << std::endl;
        std::cout << "v ";
        for(unsigned int i=1; i<=monMaxSat->nInputVars; i++) {
            std::cout << monMaxSat->getValue(i);
        }
        std::cout << std::endl;
        ///////////////////////////////////////
    }


    delete monMaxSat;
    return 0;
}



