"""
Copyright (c) <2023> <Andreas Niskanen, University of Helsinki>



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


import os, subprocess, sys, tempfile

from pysat.solvers import Solver
from pysat.formula import WCNF
from pysat.examples.rc2 import RC2
from pysat.examples.mcsls import MCSls
from pysat.examples.lbx import LBX
from pysat.examples.optux import OptUx

def lit_is_true(model, lit):
	assert(abs(model[abs(lit)-1]) == abs(lit))
	return model[abs(lit)-1] == lit

def gcnf_to_file(kb, hard, soft, file):
	n_vars = next(kb.var_counter)-1
	n_clauses = len(hard) + len(soft)
	n_groups = len(soft)
	for atom in kb.atoms:
		file.write("c " + atom + " = " + str(kb.atom_vars[atom]) + "\n")
	for i, lit in enumerate(soft):
		file.write("c F" + str(i+1) + " = " + str(lit) + "\n")
	file.write("p gcnf " + str(n_vars) + " " + str(n_clauses) + " " + str(n_groups) + "\n")
	for clause in hard:
		file.write("{0} " + " ".join(map(str, clause)) + " 0" + "\n")
	for i, lit in enumerate(soft):
		file.write("{" + str(i+1) + "} " + str(lit) + " 0" + "\n")
	file.flush()

def wcnf_to_file(kb, hard, soft, file):
	n_vars = next(kb.var_counter)-1
	n_clauses = len(hard) + len(soft)
	top = len(soft)+1
	for atom in kb.atoms:
		file.write("c " + atom + " = " + str(kb.atom_vars[atom]) + "\n")
	for i, lit in enumerate(soft):
		file.write("c F" + str(i+1) + " = " + str(lit) + "\n")
	file.write("p wcnf " + str(n_vars) + " " + str(n_clauses) + " " + str(top) + "\n")
	for i, lit in enumerate(soft):
		file.write("1 " + str(lit) + " 0" + "\n")
	for clause in hard:
		file.write(str(top) + " " + " ".join(map(str, clause)) + " 0" + "\n")
	file.flush()

def cnf_to_file(clauses, file):
    n_vars = 0
    for clause in clauses:
        n_vars = max(n_vars, max([abs(lit) for lit in clause]))
    n_clauses = len(clauses)
    file.write("p cnf " + str(n_vars) + " " + str(n_clauses) + "\n")
    for clause in clauses:
        file.write(" ".join(map(str, clause)) + " 0" + "\n")
    file.flush()

def old_wcnf_to_file(hard_clauses, soft_clauses, weights, file):
	n_vars = 0
	for clause in hard_clauses:
		n_vars = max(n_vars, max([abs(lit) for lit in clause]))
	for clause in soft_clauses:
		n_vars = max(n_vars, max([abs(lit) for lit in clause]))
	n_clauses = len(hard_clauses) + len(soft_clauses) - len([w for w in weights if w == 0])
	top = sum(weights)+1
	file.write("p wcnf " + str(n_vars) + " " + str(n_clauses) + " " + str(top) + "\n")
	for clause in hard_clauses:
		#print(str(top) + " " + " ".join(map(str, clause)) + " 0")
		file.write(str(top) + " " + " ".join(map(str, clause)) + " 0" + "\n")
	for clause, w in zip(soft_clauses, weights):
		if w == 0:
			continue
		#print(str(w) + " " + " ".join(map(str, clause)) + " 0")
		file.write(str(w) + " " + " ".join(map(str, clause)) + " 0" + "\n")
	file.flush()

def new_wcnf_to_file(hard_clauses, soft_clauses, weights, file):
	for clause in hard_clauses:
		#print("h " + " ".join(map(str, clause)) + " 0")
		file.write("h " + " ".join(map(str, clause)) + " 0" + "\n")
	for clause, w in zip(soft_clauses, weights):
		if w == 0:
			continue
		#print(str(w) + " " + " ".join(map(str, clause)) + " 0")
		file.write(str(w) + " " + " ".join(map(str, clause)) + " 0" + "\n")
	file.flush()

def maxsat_solve(hard_clauses, soft_clauses, weights, solver="rc2"):
	if solver == "rc2":
		rc2 = RC2(WCNF())
		for clause in hard_clauses:
			#print(clause)
			rc2.add_clause(clause)
		for clause, w in zip(soft_clauses, weights):
			#print(clause, w)
			if w == 0:
				continue
			rc2.add_clause(clause, weight=w)
		model = rc2.compute()
		if model is not None:
			return rc2.cost, model
		return float("inf"), None
	elif solver in ["maxhs", "uwrmaxsat", "evalmaxsat"]:
		with tempfile.NamedTemporaryFile(mode="w", suffix=".wcnf") as tmp:
			new_wcnf_to_file(hard_clauses, soft_clauses, weights, tmp)
			if solver == "maxhs":
				output = subprocess.check_output([os.path.join(os.path.dirname(__file__), "tools/maxsat/MaxHS/bin/maxhs"), "-no-printOptions", "-printSoln", "-verb=0", tmp.name]).decode("utf-8").split("\n")
			elif solver == "uwrmaxsat":
				output = subprocess.check_output([os.path.join(os.path.dirname(__file__), "tools/maxsat/UWrMaxSat/bin/uwrmaxsat"), "-v0", "-no-sat", "-no-bin", "-m", "-bm", tmp.name]).decode("utf-8").split("\n")
			elif solver == "evalmaxsat":
				output = subprocess.check_output([os.path.join(os.path.dirname(__file__), "tools/maxsat/EvalMaxSAT/bin/EvalMaxSAT_bin"), tmp.name]).decode("utf-8").split("\n")
		if "s UNSATISFIABLE" in output:
			return float("inf"), None
		elif "s OPTIMUM FOUND" in output:
			cost_line = [line for line in output if line.startswith("o ")][-1]
			cost = int(cost_line.replace("o ", "").replace("-oo", "0"))
			model_line = [line for line in output if line.startswith("v ")][-1]
			model_line = model_line.replace("v ", "")
			model = [i if model_line[i-1] == "1" else -i for i in range(1, len(model_line)+1)]
			return cost, model
		else:
			print("ERROR: MaxSAT solver did not find optimal solution. Output dump:")
			for line in output:
				print(line)
			return None, None
	elif solver in ["pacose"]:
		with tempfile.NamedTemporaryFile(mode="w", suffix=".wcnf") as tmp:
			old_wcnf_to_file(hard_clauses, soft_clauses, weights, tmp)
			if solver == "pacose":
				try:
					output = subprocess.check_output([os.path.join(os.path.dirname(__file__), "tools/maxsat/pacose/bin/Pacose"), tmp.name]).decode("utf-8").split("\n")
				except subprocess.CalledProcessError as error:
					output = error.output.decode("utf-8").split("\n")
		if "UNSATISFIABLE" in output:
			return float("inf"), None
		elif "s OPTIMUM FOUND" in output:
			cost_line = [line for line in output if line.startswith("o ")][-1]
			cost = int(cost_line.replace("o ", ""))
			model_line = [line for line in output if line.startswith("v ")][-1]
			model_line = model_line.replace("v ", "")
			model = [i if model_line[i-1] == "1" else -i for i in range(1, len(model_line)+1)]
			return cost, model
		else:
			print("ERROR: MaxSAT solver did not find optimal solution. Output dump:")
			for line in output:
				print(line)
			return None, None
	else:
		print("ERROR: Unknown solver", solver)
		sys.exit(1)

def count_assignments(hard_clauses, soft_lits, solver="glucose3"):
	sat_solver = Solver(name=solver)
	for clause in hard_clauses:
		sat_solver.add_clause(clause)
	unsat_soft_lits = soft_lits[:]

	count = 0

	while len(unsat_soft_lits) > 0:
		sat_solver.add_clause(unsat_soft_lits)
		sat_solver.set_phases(unsat_soft_lits)
		sat = sat_solver.solve()
		if not sat:
			return 0
		count += 1
		model = sat_solver.get_model()
		unsat_soft_lits = [lit for lit in unsat_soft_lits if not lit_is_true(model, lit)]

	return count


def smallest_mus(hard_clauses, soft_lits, solver="optux"): # TODO: forqes
	if solver == "optux":
		wcnf = WCNF()
		for clause in hard_clauses:
			wcnf.append(clause)
		for lit in soft_lits:
			wcnf.append([lit], weight=1)
		with OptUx(wcnf) as optux:
			smus = optux.compute()
			return smus
	elif solver == "forqes":
		with tempfile.NamedTemporaryFile(mode="w", suffix=".cnf") as tmp_hard:
			cnf_to_file(hard_clauses, tmp_hard)
			with tempfile.NamedTemporaryFile(mode="w", suffix=".cnf") as tmp_soft:
				cnf_to_file([[lit] for lit in soft_lits], tmp_soft)
				try:
					output = subprocess.check_output([os.path.join(os.path.dirname(__file__), "tools/forqes/forqes-linux-x86-64"), "-vv", "-n", tmp_hard.name, tmp_soft.name]).decode("utf-8").split("\n")
				except subprocess.CalledProcessError as error:
					output = error.output.decode("utf-8").split("\n")
		smus_line = [line for line in output if line.startswith("v ")][0]
		smus = list(map(int, smus_line.replace("v ", "").split()[:-1]))
		return smus
	else:
		print("ERROR: Unknown solver", solver)
		sys.exit(1)

def union_of_muses(kb, hard_clauses, soft_lits, solver):
	if solver == "umuser":
		with tempfile.NamedTemporaryFile(mode="w", suffix=".wcnf") as tmp:
			wcnf_to_file(kb, hard_clauses, soft_lits, tmp)
			output = subprocess.check_output([os.path.join(os.path.dirname(__file__), "tools/umuser/umuser"), tmp.name]).decode("utf-8").split("\n")
		umus_line = [line for line in output if line.startswith("c CUMU: ")][0]
		umus = list(map(int, umus_line.replace("c CUMU: ", "").split()))
		return umus
	elif solver == "mcscacheUMU":
		with tempfile.NamedTemporaryFile(mode="w", suffix=".wcnf") as tmp:
			wcnf_to_file(kb, hard_clauses, soft_lits, tmp)
			output = subprocess.check_output([os.path.join(os.path.dirname(__file__), "tools/umuser/mcscacheUMU"), tmp.name]).decode("utf-8").split("\n")
		umus_line = [line for line in output if line.startswith("c CUMU: ")][0]
		umus = list(map(int, umus_line.replace("c CUMU: ", "").split()))
		return umus
	else:
		print("ERROR: Unknown solver", solver)
		sys.exit(1)

def enumerate_muses(kb, hard_clauses, soft_lits, solver):
	if solver == "optux":
		wcnf = WCNF()
		for clause in hard_clauses:
			wcnf.append(clause)
		for lit in soft_lits:
			wcnf.append([lit], weight=1)
		with OptUx(wcnf, unsorted=True) as optux:
			for mus in optux.enumerate():
				yield mus
	elif solver == "marco":
		with tempfile.NamedTemporaryFile(mode="w", suffix=".gcnf") as tmp:
			gcnf_to_file(kb, hard_clauses, soft_lits, tmp)
			proc = subprocess.Popen(["python3", os.path.join(os.path.dirname(__file__), "tools/MARCO/marco.py"), "--parallel", "MUS", "-v", tmp.name], stdout=subprocess.PIPE)
			while True:
				line = proc.stdout.readline()
				if not line:
					break
				line = line.rstrip().decode("utf-8")
				if line.startswith("U "):
					mus = list(map(int, line.replace("U ", "").split()))
					yield mus
	elif solver == "unimus": # TODO: unimus seems to report incorrect numbers
		with tempfile.NamedTemporaryFile(mode="w", suffix=".gcnf") as tmp:
			gcnf_to_file(kb, hard_clauses, soft_lits, tmp)
			proc = subprocess.Popen([os.path.join(os.path.dirname(__file__), "tools/unimus/unimus"), "-v", "1", tmp.name], stdout=subprocess.PIPE)
			while True:
				line = proc.stdout.readline()
				if not line:
					break
				line = line.rstrip().decode("utf-8")
				if line.startswith("MUS "):
					mus = list(map(int, line.replace("MUS ", "").split()))
					mus = [i+1 for i in mus]
					yield mus
	else:
		print("ERROR: Unknown solver", solver)
		sys.exit(1)

def enumerate_mcses(kb, hard_clauses, soft_lits, solver):
	if solver == "mcsls": # TODO: mcsls seems to report incorrect numbers
		wcnf = WCNF()
		for clause in hard_clauses:
			wcnf.append(clause)
		for lit in soft_lits:
			wcnf.append([lit], weight=1)
		with MCSls(wcnf) as mcsls:
			for mcs in mcsls.enumerate():
				mcsls.block(mcs)
				yield mcs
	elif solver == "lbx":
		wcnf = WCNF()
		for clause in hard_clauses:
			wcnf.append(clause)
		for lit in soft_lits:
			wcnf.append([lit], weight=1)
		with LBX(wcnf) as lbx:
			for mcs in lbx.enumerate():
				lbx.block(mcs)
				yield mcs
	elif solver == "marcoMCS":
		with tempfile.NamedTemporaryFile(mode="w", suffix=".gcnf") as tmp:
			gcnf_to_file(kb, hard_clauses, soft_lits, tmp)
			proc = subprocess.Popen(["python3", os.path.join(os.path.dirname(__file__), "tools/MARCO/marco.py"), "--parallel", "MCS", "--print-mcses", "-v", tmp.name], stdout=subprocess.PIPE)
			while True:
				line = proc.stdout.readline()
				if not line:
					break
				line = line.rstrip().decode("utf-8")
				if line.startswith("C "):
					mcs = list(map(int, line.replace("C ", "").split()))
					yield mcs
	elif solver == "mcscache":
		with tempfile.NamedTemporaryFile(mode="w", suffix=".wcnf") as tmp:
			wcnf_to_file(kb, hard_clauses, soft_lits, tmp)
			proc = subprocess.Popen([os.path.join(os.path.dirname(__file__), "tools/mcscache/mcscache"), tmp.name], stdout=subprocess.PIPE)
			while True:
				line = proc.stdout.readline()
				if not line:
					break
				line = line.rstrip().decode("utf-8")
				if line.startswith("c MCS: "):
					mcs = list(map(int, line.replace("c MCS: ", "").split()))
					yield mcs
	elif solver == "rime":
		with tempfile.NamedTemporaryFile(mode="w", suffix=".wcnf") as tmp:
			wcnf_to_file(kb, hard_clauses, soft_lits, tmp)
			proc = subprocess.Popen([os.path.join(os.path.dirname(__file__), "tools/rime/rime"), "-v", "1", tmp.name], stdout=subprocess.PIPE)
			while True:
				line = proc.stdout.readline()
				if not line:
					break
				line = line.rstrip().decode("utf-8")
				if line.startswith("MCS "):
					mcs = list(map(int, line.replace("MCS ", "").split()))
					mcs = [i+1 for i in mcs]
					yield mcs
	else:
		print("ERROR: Unknown solver", solver)
		sys.exit(1)

def count_muses(kb, hard_clauses, soft_lits):
	with tempfile.NamedTemporaryFile(mode="w", suffix=".gcnf") as tmp:
		gcnf_to_file(kb, hard_clauses, soft_lits, tmp)
		output = subprocess.check_output(["python3", "counter.py", tmp.name], cwd=os.path.join(os.path.dirname(__file__), "tools/exactMUSCounter")).decode("utf-8").split("\n")
	mus_line = [line for line in output if line.startswith("MUS count:")][0]
	return int(mus_line.replace("MUS count: ", ""))

def count_mcses(kb, hard_clauses, soft_lits):
	with tempfile.NamedTemporaryFile(mode="w", suffix=".gcnf") as tmp:
		gcnf_to_file(kb, hard_clauses, soft_lits, tmp)
		output = subprocess.check_output(["python3", "counter.py", tmp.name], cwd=os.path.join(os.path.dirname(__file__), "tools/MSSCounting")).decode("utf-8").split("\n")
	mss_line = [line for line in output if line.startswith("MSS count:")][0]
	return int(mss_line.replace("MSS count: ", ""))

def exists_sc_formula(hard_clauses, soft_lits, solver="glucose3"):
	sat_solver = Solver(name=solver)
	for clause in hard_clauses:
		sat_solver.add_clause(clause)
	for lit in soft_lits:
		if not sat_solver.solve([lit]):
			return True
	return False

def count_sc_formulas(hard_clauses, soft_lits, solver="glucose3"):
	count = 0
	sat_solver = Solver(name=solver)
	for clause in hard_clauses:
		sat_solver.add_clause(clause)
	for lit in soft_lits:
		if not sat_solver.solve([lit]):
			count += 1
	return count

def get_entailed(hard_clauses, soft_lits, solver="glucose3"):
	sat_solver = Solver(name=solver)
	for clause in hard_clauses:
		sat_solver.add_clause(clause)
	entailed = set()
	for i, ilit in enumerate(soft_lits):
		if i in entailed:
			continue
		for j, jlit in enumerate(soft_lits):
			if i == j or j in entailed:
				continue
			if not sat_solver.solve([ilit, -jlit]):
				entailed.add(j)
	return entailed
