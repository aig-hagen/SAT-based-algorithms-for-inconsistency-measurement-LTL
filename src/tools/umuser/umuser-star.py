# umuser-star.py
# Author: Carlos Mencia

import sys
import os;

class umuserStar:

	def __init__(self, ident, tmcs, tumuser, fname):
		self.fname = fname
		self.ident = ident
		self.tmcs = int(tmcs);
		self.tumuser = int(tumuser);

		self.cumu_mcs = []
		self.solved_mcs = False;
		self.cputime_mcs = "";

		self.cumu_umuser = []
		self.solved_umuser = False;
		self.cputime_umuser = "";

	
	def print_header(self):
		print "c *** umuser*: Computing the Union of MUSes ***"
		print "c *** instance: " + self.fname + " ***"
		print "c"

			
	def solve(self):
		self.print_header();
		self.invoke_mcscache();
		self.process_file_mcs();
		self.report_cumu(self.cumu_mcs, self.cputime_mcs, " (mcscache)");
		if self.solved_mcs == True:
			print "c Instance solved by mcscache"
			
		else:
			print "c mcscache timed out ..."
			self.create_mcsfile();
			self.invoke_umuser();
			self.process_file_umuser();
			self.report_cumu(self.cumu_umuser, self.cputime_umuser, " (umuser)")
			if self.solved_umuser == True:
				print "c Instance solved by umuser"
			else:
				print "c umuser timed out ..."

		print "c Terminating ..."
		


	def invoke_mcscache(self):
		auxname = self.ident + ".mcscache"
		cmd = "./mcscacheUMU -T " + str(self.tmcs) + " " + self.fname + " > " + auxname;
		print "c Running mcscache ...";
		os.system(cmd);

	def process_file_mcs(self):
		auxname = self.ident + ".mcscache"
		if not os.path.isfile(auxname):
			print "c Error running mcscache. Terminating ..."
			exit();

		ok = False;
		for line in open(auxname):
			if "c CUMU:" in line:
				self.cumu_mcs = line.split()[2:]
				ok = True
			elif "c Instance solved" in line:
				self.solved_mcs = True;
			elif "c CPU Time" in line:
				self.cputime_mcs = line.split()[3];

		if not ok:
			print "c Error running mcscache. Terminating ..."
			exit();

	def process_file_umuser(self):
		auxname = self.ident + ".umuser"
		if not os.path.isfile(auxname):
			print "c Error running umuser. Terminating ..."
			exit();

		ok = False;
		for line in open(auxname):
			if "c CUMU:" in line:
				self.cumu_umuser = line.split()[2:]
				ok = True;
			elif "c Instance solved" in line:
				self.solved_umuser = True;
			elif "c CPU Time" in line:
				self.cputime_umuser = line.split()[3];

		if not ok:
			print "c Error running umuser. Terminating ..."
			exit();

	def report_cumu(self, cumu, cputime, string):
		st = "c CUMU" + string + ": "
		for c in cumu:
			st += c + " "
		print st;
		
		print "c CUMU size " + string + ": " + str(len(cumu));
		print "c CPU Time " + string + ": " + cputime;


	def create_mcsfile(self):
		st = ""
		for c in self.cumu_mcs:
			st += c + " "
		st += "0";

		auxfname = self.ident + ".mcs"
		f = open(auxfname, "w+")
		f.write(st);
		f.close();
		
	
	def invoke_umuser(self):
		mcsfile = self.ident + ".mcs"
		umuserfile = self.ident + ".umuser"
		cmd = "./umuser -bt " + mcsfile + " -T " + str(self.tumuser) + " " + self.fname + " > " + umuserfile
		print "c"
		print "c Running umuser ..."
		os.system(cmd);

	def finishup(self):
		if os.path.isfile(self.ident + ".mcscache"):
			os.system("rm " + self.ident + ".mcscache");
		if os.path.isfile(self.ident + ".umuser"):
			os.system("rm " + self.ident + ".umuser");
		if os.path.isfile(self.ident + ".mcs"):
			os.system("rm " + self.ident + ".mcs");




def print_help():
	print "umuser*: Computing the Union of MUSes"
	print ""
	print "usage: python umuser-star.py <ident> <to_mcscache> <to_umuser> <input>"
	print "where: "
	print "  <ident>:       unique identifier for the execution"
	print "  <to_mcscache>: Timeout in seconds for mcscache (integer >= 1)"
	print "  <to_umuser>:   Timeout in seconds for umuser (integer >= 1)"
	print "  <input>:       Input CNF file"
	print ""
	print "additional info: see distribution README or contact main author"
	print "author: Carlos Mencia (cmencia@gmail.com)"
	

def check_cmdline():
	if len(sys.argv) != 5:
		print_help();
		exit();

	ok = True;
	try:
		t1 = int(sys.argv[2])
		t2 = int(sys.argv[3])
		if(t1 <= 0 or t2 <= 0):
			ok = False;
	except:
		ok = False

	if not os.path.isfile(sys.argv[4]):
		ok = False

	if not ok:
		print_help();
		exit();
		


def main():
	check_cmdline()
	us = umuserStar(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]);
	us.solve();
	us.finishup();


main();



