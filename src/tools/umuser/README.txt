
-------------------------------------------------------------------------------
On Computing the Union of MUSes


This folder contains the tools proposed in [1] for the task of computing 
the union of MUSes of a given unsatisfiable CNF formula.

It contains the following files:

  1) umuser: Linux executable with the main tool

  2) mcscacheUMU: Linux executable with the MCS enumerator mcscache

  3) umuser-star.py: Python 2.7 script implementing the bootstrapped version
      of umuser, termed umuser*, running first mcscacheUMU and then umuser.

  4) ex.cnf: The CNF formula used in Example 1 in [1]
       ------------
	p cnf 4 7
	1 0
	-1 0
	-2 0
	-1 2 0
	-2 1 0
	-3 4 0
	-4 3 0
       -------------


 Usage instructions:
---------------------

1) umuser

usage: umuser [ <option> ... ] <input> 
where <input> is a CNF formula and <option> is one of the following:
  -h          prints this help and exits
  -T   TTT    specify timeout [default: -T 3600]
  -bt <file>  bootstrap with clauses in <file> [default: off]

The option -bt serves for initializing umuser with some clauses known to
be in the union of MUSes. The specified file (<file>), must contain a single
line with a list of numbers (specifying the clauses), and ending with "0".


As an example:

> ./umuser -T 10 ex.cnf

This command specifies a timeout of 10 seconds, and the ouput is:

  c *** umuser: Computing the Union of MUSes ***
  c *** instance: ex.cnf ***
  c 
  c Parsing ...
  c Parsing CPU Time: 0
  c Running umuser ...
  c CUMU: 1 2 3 4 
  c 
  c CUMU size: 4
  c Instance solved. Terminating umuser ...
  c CPU Time: 0

This indicates that the instance was solved, and the union of MUSes are
clauses 1, 2, 3 and 4.

If umuser times out, CUMU would represent an under-approximation of the
union of MUSes.

As another example: for bootstrapping umuser with clauses 2 and 3, it would
be necessary to create a file, e.g. f.txt with the following content:

----------------- 
2 3 0
-----------------

Then umuser should be invoked with the following command:

> ./umuser -bt f.txt ex.cnf



2) mcscacheUMU

usage: mcscacheUMU [ <option> ... ] <input> 
where <option> is one of the following:
  -h        prints this help and exits
  -T   TTT  specify timeout [default: -T 3600]

The output follows the same structure as in umuser, indicating the clauses
in CUMU when it terminates.



3) umuser-star.py

usage: python umuser-star.py <ident> <to_mcscache> <to_umuser> <input>
where: 
  <ident>:       unique identifier for the execution
  <to_mcscache>: Timeout in seconds for mcscache (integer >= 1)
  <to_umuser>:   Timeout in seconds for umuser (integer >= 1)
  <input>:       Input CNF file


This script implements umuser*, running mcscacheUMU for at most <to_mcscache>
seconds and, then initializing umuser with the clauses found and running it for
at most <to_umuser> seconds.

If the instance is solved in the first phase umuser is not invoked.

The first argument, <ident>, is a unique identifier for each run. This identifier
is necessary since umuser-star.py writes some auxiliary files during its execution.
These files contain the provided identifier in their names, and are deleted when the
tool terminates.

Note: The script should be run with Python 2.7.


As an example:

> python umuser-star.py id1 3 10 ex.cnf

With this command, umuser* would run mcscacheUMU for 3 seconds and, if it times 
out, it would run the bootstrapped umuser for 10 seconds.

The tool will write the following files in the same directory: id1.mcscache, id1.mcs 
and id1.umuser, which will be removed at the end.

The output would be the following:

  c *** umuser*: Computing the Union of MUSes ***
  c *** instance: ex.cnf ***
  c
  c Running mcscache ...
  c CUMU (mcscache): 1 2 3 4 
  c CUMU size  (mcscache): 4
  c CPU Time  (mcscache): 0
  c Instance solved by mcscache
  c Terminating ...

In this case, the instance is solved in the first phase, so umuser 
is not invoked.

For a harder instance the output would be something as:

  c *** umuser*: Computing the Union of MUSes ***
  c *** instance: fml.cnf ***
  c
  c Running mcscache ...
  c CUMU (mcscache): 2 195 3 4 969 5 6 1164 7 8 9 10 11 12 13 14 ... 
  c CUMU size  (mcscache): 146
  c CPU Time  (mcscache): 2.944
  c mcscache timed out ...
  c
  c Running umuser ...
  c CUMU (umuser): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
  c CUMU size  (umuser): 1218
  c CPU Time  (umuser): 4.328
  c Instance solved by umuser
  c Terminating ...



 References:
-------------

[1] Carlos Mencía, Oliver Kullmann, Alexey Ignatiev, Joao Marques-Silva.
    On Computing the Union of MUSes. SAT 2019: pp. 211-221.

-------------------------------------------------------------------------------
