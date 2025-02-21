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

from pysat.solvers import Solver
from pysat.formula import WCNF
from pysat.examples.lbx import LBX

import itertools
import sys
import core

def lit_is_true(model, lit):
	return model[abs(lit)-1] == lit

def disjoint_phase(solver, soft_lits):
	assumptions = soft_lits[:]
	while not solver.solve(assumptions):
		core = solver.get_core()
		i = 0
		while i < len(core):
			to_test = core[:i] + core[(i+1):]
			if solver.solve(to_test):
				i += 1
			else:
				core = to_test
		core = set(core)
		assumptions = [lit for lit in assumptions if lit not in core]
	return assumptions

def mcs_overlap(solver, hard_clauses, soft_lits, query, maximize=True):
	top = solver.nof_vars()
	top = max(top, max([abs(lit) for lit in soft_lits]))

	q = top+1
	solver.add_clause([-q] + [-lit for lit in query])
	for lit in query:
		solver.add_clause([q, lit])

	iters = 0

	while solver.solve([q]):
		iters += 1
		model = solver.get_model()

		if maximize:
			solver.add_clause([lit for lit in soft_lits if not lit_is_true(model, lit)])
			while solver.solve([q] + [lit for lit in soft_lits if lit_is_true(model, lit)]):
				model = solver.get_model()
				solver.add_clause([lit for lit in soft_lits if not lit_is_true(model, lit)])

		if not solver.solve([-q] + [lit for lit in soft_lits if lit_is_true(model, lit)]):
			mcs = [lit for lit in soft_lits if not lit_is_true(model, lit)]

			if not maximize:
				solver.add_clause(mcs)
				while solver.solve([lit for lit in soft_lits if lit_is_true(model,lit)]):
					model = solver.get_model()
					mcs = [lit for lit in soft_lits if not lit_is_true(model, lit)]
					solver.add_clause(mcs)

			print("c CEGAR:", iters, "iterations")
			return mcs

		model = solver.get_model()
		solver.add_clause([lit for lit in soft_lits if not lit_is_true(model, lit)])

	print("c CEGAR:", iters, "iterations")
	return None

def compute_autarky(hard_clauses, soft_lits):
	n = 0
	if len(soft_lits) > 0:
		n = max(n, max([abs(lit) for lit in soft_lits]))
	for clause in hard_clauses:
		n = max(n, max([abs(lit) for lit in clause]))
	var_counter = itertools.count(1)
	tvars = { var : next(var_counter) for var in range(1,n+1) }
	fvars = { var : next(var_counter) for var in range(1,n+1) }
	wcnf = WCNF()
	for clause in hard_clauses:
		for l1 in clause:
			if l1 > 0:
				wcnf.append([-fvars[abs(l1)]] + [tvars[abs(l2)] for l2 in clause if l2 > 0 and l2 != l1] + [fvars[abs(l2)] for l2 in clause if l2 < 0])
			else:
				wcnf.append([-tvars[abs(l1)]] + [tvars[abs(l2)] for l2 in clause if l2 > 0] + [fvars[abs(l2)] for l2 in clause if l2 < 0 and l2 != l1])
	for lit in soft_lits:
		if lit > 0:
			wcnf.append([-fvars[abs(lit)]])
		else:
			wcnf.append([-tvars[abs(lit)]])
	for var in range(1,n+1):
		wcnf.append([-tvars[var], -fvars[var]])
	for lit in soft_lits:
		if lit > 0:
			wcnf.append([tvars[abs(lit)]], weight=1)
		else:
			wcnf.append([fvars[abs(lit)]], weight=1)
	solver = LBX(wcnf, use_cld=True)
	model = set(solver.compute())
	return [lit for i, lit in enumerate(soft_lits) if i+1 not in model]


def p_measure_ltl(kb, hard, soft, trim_autarky=True, disjoint_cores=True, maximize=True):
	autarky = set()
	softs_remaining = soft[:]
	print("c Original:", len(softs_remaining), "softs remaining")
	if trim_autarky:
		autarky = compute_autarky(hard, soft)
		softs_remaining = [lit for lit in soft if lit not in autarky]
		print("c Trim autarky:", len(softs_remaining), "softs remaining")
	solver = Solver(name="cadical153")
	for clause in hard:
		solver.add_clause(clause)
	if disjoint_cores:
		softs_remaining = disjoint_phase(solver, softs_remaining)
		print("c Disjoint phase:", len(softs_remaining), "softs remaining")
	iters = 0
	while True:
		iters += 1
		query = softs_remaining
		mcs = mcs_overlap(solver, hard, soft, query, maximize)
		if mcs is None:
			break
		mcs = set(mcs)
		softs_remaining = [lit for lit in softs_remaining if lit not in mcs]
		print("c Iteration", iters, ":", len(softs_remaining), "softs remaining")
	return len(soft)-len(softs_remaining)-len(autarky)


def mv_ltl_measure(kb, hard, soft, original_sig_size, trim_autarky=True, disjoint_cores=True, maximize=True):
	original_atoms_covered = set()
	autarky = set()
	softs_remaining = soft[:]
	print("c Original:", len(softs_remaining), "softs remaining")

	if trim_autarky:
		autarky = set(compute_autarky(hard, soft))
		softs_remaining = [lit for lit in soft if lit not in autarky]
		print("c Trim autarky:", len(softs_remaining), "softs remaining")

	solver = Solver(name="cadical153")
	for clause in hard:
		solver.add_clause(clause)
	atoms_remaining = [atom for atom in kb.atoms]

	print("c Original:", len(atoms_remaining), "atoms remaining")
	if disjoint_cores:
		softs_remaining = disjoint_phase(solver, softs_remaining)
		print("c Disjoint phase:", len(softs_remaining), "softs remaining")
		softs_remaining = set(softs_remaining)
		atoms_covered = set()
		for i, lit in enumerate(soft):
			if lit not in softs_remaining and lit not in autarky:
				atoms_covered.update(kb.formulas[i].calculate_atoms())
				original_atoms_covered.update(kb.formulas[i].atoms)
		atoms_remaining = [atom for atom in atoms_remaining if atom not in atoms_covered]
		print("c Disjoint phase:", len(atoms_remaining), "atoms remaining")

	atom_to_occurrences = { atom : set(lit for i, lit in enumerate(soft) if atom in kb.formulas[i].calculate_atoms()) for atom in kb.atoms }

	iters = 0
	while True:
		iters += 1
		print(iters)
		query = set()
		for atom in atoms_remaining:
			query.update(atom_to_occurrences[atom])
		query = query.difference(autarky)
		query = list(query)
		mcs = mcs_overlap(solver, hard, soft, query, maximize)
		if mcs is None:
			break
		original_atoms_covered = set()
		for i, lit in enumerate(soft):
			if lit in mcs:
				original_atoms_covered.update(kb.formulas[i].atoms)

		tmp_atoms_remaining = set()
		for at in atoms_remaining:
			curr_original_atom = str(at).rsplit('_',1)[0]
			if curr_original_atom not in original_atoms_covered:
				tmp_atoms_remaining.add(at)
		atoms_remaining = tmp_atoms_remaining
		print("c Iteration", iters, ":", len(atoms_remaining), "atoms remaining")

	return (len(original_atoms_covered)/original_sig_size)


if __name__ == "__main__":
	filename = sys.argv[1]
	kb = core.KnowledgeBase(open(filename).read().split("\n"))
	print("v", p_measure(kb))
