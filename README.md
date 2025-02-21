# SAT-based-algorithms-for-inconsistency-measurement-LTL

### MaxSAT-based implementations of the LTL time inconsistency measure and the LTL contension inconsistency measure

Usage: ``` python3 main.py <measure> --ltl-m <m> [--solver=<solver_name>]```
* Available measures:
    * ```t-ltl```: LTL time inconsistency measure
    * ```c-ltl```: LTL contension inconsistency measure
* ```<m>``` is the index of the final temporal state (i.e., any number >= 0)

### CEGAR-based implementations of the LTL problematic inconsistency measure and the LTL MUS-variable-based inconsistency measure
Usage: ``` python3 main.py <measure> --ltl-m <m> --solver=cegar [--trim-autarky] [--disjoint-cores] [--cegar-maximize]```
* Available measures:
    * ```t-ltl```: LTL time inconsistency measure
    * ```c-ltl```: LTL contension inconsistency measure
* ```<m>``` is the index of the final temporal state (i.e., any number >= 0)

### Acknowledgements: 
This code is based on the implementations by Andreas Niskanen: ```https://bitbucket.org/coreo-group/sat4im/```
See also the corresponding papers: 
* _MaxSAT-Based Inconsistency Measurement._ Andreas Niskanen, Isabelle Kuhlmann, Matthias Thimm, and Matti Järvisalo. Proceedings of the 26th European Conference on Artificial Intelligence (ECAI 2023).
* _Computing MUS-Based Inconsistency Measures._ Isabelle Kuhlmann, Andreas Niskanen, and Matti Järvisalo. Proceedings of the 18th European Conference on Logics in Artificial Intelligence (JELIA 2023). 
