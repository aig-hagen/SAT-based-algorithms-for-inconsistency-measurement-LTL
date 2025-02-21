"""
Copyright (c) <2023> <Andreas Niskanen, University of Helsinki>
Copyright (c) <2025> <Isabelle Kuhlmann, University of Hagen>


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:



The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.



THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import core, encoder, solver, cegar
import argparse, os, sys, itertools, resource

def get_utime():
	self_utime = resource.getrusage(resource.RUSAGE_SELF).ru_utime
	children_utime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
	return self_utime + children_utime

measures = ["t-ltl", "c-ltl", "p-ltl", "mv-ltl"]

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="SAT-based Inconsistency Measurement")
	parser.add_argument("filename", help="Input filename for KB instance in Tweety format.")
	parser.add_argument("measure", choices=measures, help="Inconsistency measure. Allowed values are " + ", ".join(measures), metavar="measure")
	parser.add_argument("--format", choices=["pl", "cnf", "gcnf"], default="pl", help="Input file format.")
	parser.add_argument("--output-gcnf", default=None, help="Output a PL file to GCNF and exit.")
	parser.add_argument("--output-wcnf", default=None, help="Output the MaxSAT encoding to WCNF and exit.")
	parser.add_argument("--solver", help="Override default solver.")
	parser.add_argument("--lp-solver", default="GLOP", help="LP solver (column generation for the eta measure).")
	parser.add_argument("--trim-autarky", default=False, action="store_true", help="Compute the maximal autarky and remove it from softs (p, mv).")
	parser.add_argument("--disjoint-cores", default=False, action="store_true", help="Compute a set of disjoint MUSes and remove them from softs (p, mv).")
	parser.add_argument("--cegar-maximize", default=False, action="store_true", help="Always perform subset-maximization in CEGAR (p, mv).")
	parser.add_argument("--ltl-m", default=5, help="Index of final state (for LTLff inconsistency measures)")
	args = parser.parse_args()
	start_time = get_utime()
	filename = args.filename
	measure = args.measure
	contents = [line for line in open(filename).read().split("\n") if len(line) > 0]
	kb = core.KnowledgeBase(contents, args.format)
	if args.format == "pl" and args.output_gcnf is not None:
		hard, soft = kb.to_group_cnf("global")
		with open(args.output_gcnf, "w") as file:
			solver.gcnf_to_file(kb, hard, soft, file)
		sys.exit(0)

	kb_time = get_utime()
	print("c Parsing: ", round(kb_time-start_time, 2), "seconds", flush=True)

	if measure in ["t-ltl", "c-ltl"]:

		if measure == "t-ltl":
			hard, soft = encoder.drastic_ltl_encoding(kb,args.ltl_m)

		elif measure == "c-ltl":
			hard, soft = encoder.contension_ltl_encoding(kb,args.ltl_m)

		enc_time = get_utime()
		print("c Encoding:", round(enc_time-kb_time, 2), "seconds", flush=True)

		if args.output_wcnf is not None:
			with open(args.output_wcnf, "w") as file:
				solver.new_wcnf_to_file(hard, soft, [1]*len(soft), file)
			sys.exit(0)
		cost, _ = solver.maxsat_solve(hard, soft, len(soft)*[1], "rc2" if args.solver is None else args.solver)
		end_time = get_utime()
		print("c Solving: ", round(end_time-enc_time, 2), "seconds")
		print("o", cost)

	elif measure in ["p-ltl", "mv-ltl"]:
		original_sig_size = len(kb.atoms)
		kb.ltl_kb_to_pl(args.ltl_m)
		hard, soft = kb.to_group_cnf("global")
		enc_time = get_utime()
		print("c Encoding:", round(enc_time-kb_time, 2), "seconds", flush=True)
		if measure == "p-ltl":
			if args.solver == "cegar":
				val = cegar.p_measure_ltl(kb, hard, soft, args.trim_autarky, args.disjoint_cores, args.cegar_maximize)
			elif args.solver in ["umuser", "mcscacheUMU"]:
				umus = solver.union_of_muses(kb, hard, soft, args.solver)
				val = len(umus)
			elif args.solver in ["marco", "unimus", "optux"]:
				umus = set()
				for mus in solver.enumerate_muses(kb, hard, soft, args.solver):
					umus.update(mus)
					if len(umus) == len(kb.formulas):
						break
				val = len(umus)
			else:
				umcs = set()
				for mcs in solver.enumerate_mcses(kb, hard, soft, "lbx" if args.solver is None else args.solver):
					umcs.update(mcs)
					if len(umcs) == len(kb.formulas):
						break
				val = len(umcs)

		elif measure == "mv-ltl":
			val = cegar.mv_ltl_measure(kb, hard, soft, original_sig_size, args.trim_autarky, args.disjoint_cores, args.cegar_maximize)

		end_time = get_utime()
		print("c Solving: ", round(end_time-enc_time, 2), "seconds")
		print("o", val)
