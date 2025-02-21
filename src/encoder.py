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


from core import *
from solver import exists_sc_formula
from pysat.card import ITotalizer
import itertools

def contension_ltl_encoding(kb,m):
	kb.ltl_kb_to_pl(m)

	hard, soft_lits = kb.to_group_cnf("unique")

	for lit in soft_lits:
		hard.append([lit])
	soft = []

	for i in range(0, int(m)+1):

		for atom in kb.occs.keys():

			equiv_var = next(kb.var_counter)

			atom_state = f"{atom}_{i}"
			if atom_state in kb.atoms:
				for occ_var in kb.occurrences[atom_state]:
					hard.append([-equiv_var, -kb.atom_vars[atom_state], occ_var])
					hard.append([-equiv_var, kb.atom_vars[atom_state], -occ_var])

			soft.append([equiv_var])

	return hard, soft


def drastic_ltl_encoding(kb,m):

	kb.ltl_kb_to_pl(m)

	hard, soft_lits = kb.to_group_cnf("unique")

	for lit in soft_lits:
		hard.append([lit])
	soft = []

	for i in range(0, int(m)+1):

		equiv_var = next(kb.var_counter)

		for atom in kb.occs.keys():
			atom_state = f"{atom}_{i}"
			if atom_state in kb.atoms:
				for occ_var in kb.occurrences[atom_state]:
					hard.append([-equiv_var, -kb.atom_vars[atom_state], occ_var])
					hard.append([-equiv_var, kb.atom_vars[atom_state], -occ_var])

		soft.append([equiv_var])

	return hard, soft
