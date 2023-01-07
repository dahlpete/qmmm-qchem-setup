import numpy as np
import math
import MDAnalysis as mda
import sys
import time

class residue:
	def __init__(self,rname,rid,chn,aname):
		self.rname = rname
		self.rid = rid
		self.chn = chn
		self.aname = aname

	def __eq__(self,other):
		if isinstance(other, self.__class__):
			return self.__dict__==other.__dict__
		else:
			return False

	def __ne__(self,other):
		return not self.__eq__(other)

	def __repr__(self):
		try:
			return repr((self.rname, self.rid, self.chn, self.aname))
		except:
			return repr((self.rname, self.rid, self.chn))
		finally:
			return repr((self.rid, self.chn))


def text2list(text_filename):
	'''
	This function reads a text file in which each residue of the QM
	region is listed on its own line. The line has the following format:
	(resname)        (resid)        (chain)      (atomname)

	OUTPUT: The function returns a list whose elements are composed of 
	  tuples containing the residue information for each qm member
	'''

	fileIN = open(text_filename).readlines()
	qm_list = [[] for i in range(len(fileIN))]
	for lnum in range(len(fileIN)):
		line = fileIN[lnum]
		values = line.split()
		qm_list[lnum] = residue(values[0],values[1],values[2],values[3])
		
	return qm_list


def selection_text(qm_list,mode='bb',cutmode='cc',check_bb=False):
	'''
	This function takes the list of qm residues and generates a string
	to be used as MDAnalysis selection text. When bb is set to True, the
	function will perform qm/mm cuts according to the cutmode setting.
	When bb is set to False, all cuts will be between alpha and beta 
	carbon atoms.

	OUTPUT: selection text string for select_atoms mda function
	'''

	non_prot = ['HER','HEO']
	bb_atoms_charmm = ['CA','N','C','O','HN']

	def uniquity(starting_list):
		ulist = [starting_list[0]]
		for i in range(1,len(starting_list)):
			if not (starting_list[i] in ulist):
				ulist.append(starting_list[i])
		return ulist

	def sort_list(qm_list):
		if check_bb:
			print('check bb: this is a work in progress')

		else:
			# to reduce the size to unique residues we must eliminate the atom names
			for q in qm_list:
				del q.aname

			qm_list.sort(key = lambda resd: int(resd.rid))
			urid = uniquity(qm_list) # gets the unique residues

			cnt = 0; not_prot = [False for i in range(len(urid))]
			for u in urid:
				if (u.rname in non_prot):
					not_prot[cnt] = True
				del u.rname
				cnt += 1
			#urid.sort(key = lambda resd: int(resd.rid))

		return urid,not_prot

	
	unique_res,not_prot = sort_list(qm_list)
	if mode == 'bb':
		if cutmode == 'cc':
			caplist = []; count = 0
			for res in unique_res:

				resm1 = str(int(res.rid) - 1)
				resp1 = str(int(res.rid) + 1)
				if count == 0:
					seltext = '(resid %s and segid %s and name C O or (resid %s and segid %s and not name C O)) ' % (resm1,res.chn,res.rid,res.chn)
				else:
					seltext += 'or (resid %s and segid %s and name C O or (resid %s and segid %s and not name C O)) ' % (resm1,res.chn,res.rid,res.chn)	
				
				# Check if resm1 is in the reslist
				m1_class = residue('dummy',resm1,res.chn,'dummy')
				del m1_class.aname; del m1_class.rname
				if not (m1_class in unique_res) and not not_prot[count]:
					# if it is not there, create the captext for the prior residue
					# had it been there, the prior iteration of the loop would have already accounted for those atoms
					#print(not_prot[count])
					captext_mm = 'resid %s and segid %s and name CA' % (resm1,res.chn)
					captext_qm = 'resid %s and segid %s and name C' % (resm1,res.chn)
					captext = [captext_qm,captext_mm]
					caplist.append(captext)

				# Check if resp1 is in the reslist
				p1_class = residue('dummy',resp1,res.chn,'dummy')
				del p1_class.aname; del p1_class.rname
				if (p1_class in unique_res) and not not_prot[count]:
					# if it is there, append C and O atoms to seltext
					seltext += 'or (resid %s and segid %s and name C O) ' % (res.rid,res.chn)					

				elif not (p1_class in unique_res) and not not_prot[count]:
					# if it is not there, create the captext for the current residue
					captext_mm = 'resid %s and segid %s and name C' % (res.rid,res.chn)
					captext_qm = 'resid %s and segid %s and name CA' % (res.rid,res.chn)
					captext = [captext_qm,captext_mm]
					caplist.append(captext)

				count += 1

		if cutmode == 'nc':
			caplist = []; count = 0
			for res in unique_res:

				resm1 = str(int(res.rid) - 1)
				resp1 = str(int(res.rid) + 1)
				if count == 0:
					seltext = '(resid %s and segid %s) ' % (res.rid,res.chn)
				else:
					seltext += 'or (resid %s and segid %s) ' % (res.rid,res.chn)

				# Check if resm1 is in the reslist
				m1_class = residue('dummy',resm1,res.chn,'dummy')
				del m1_class.aname; del m1_class.rname
				if not (m1_class in unique_res) and not not_prot[count]:
					# if it is not there, create the captext bridging the prior residue
					captext_mm = 'resid %s and segid %s and name C' % (resm1,res.chn)
					captext_qm = 'resid %s and segid %s and name N' % (res.rid,res.chn)
					captext = [captext_qm,captext_mm]
					caplist.append(captext)

					try:
						zerocharge += ' or (resid %s and segid %s and name N HN HT1 HT2 HT3 HA C O CA)' % (resm1,res.chn)
					except:
						zerocharge = '(resid %s and segid %s and name N HN HT1 HT2 HT3 HA C O CA)' % (resm1,res.chn)

				#Check if resp1 is in the reslist
				p1_class = residue('dummy',resp1,res.chn,'dummy')
				del p1_class.aname; del p1_class.rname
				if not (p1_class in unique_res) and not not_prot[count]:
					#if it is not there, create the captext bridging the next residue
					captext_mm = 'resid %s and segid %s and name N' % (resp1,res.chn)
					captext_qm = 'resid %s and segid %s and name C' % (res.rid,res.chn)
					captext = [captext_qm,captext_mm]
					caplist.append(captext)

					try:
						zerocharge += ' or (resid %s and segid %s and name N HN HT1 HT2 HT3 HA C O CA)' % (resp1,res.chn)
					except:
						zercharge = '(resid %s and segid %s and name N HN HT1 HT2 HT3 HA C O CA)' % (resp1,res.chn)

				count += 1
		return seltext,caplist,zerocharge


	if mode == 'sc':
		caplist = []; count = 0
		for res in unique_res:
			if count == 0:
				if not_prot[count]:
					seltext = '(resid %s and segid %s and not name C O N CA HN HT1 HT2 HT3)' % (res.rid,res.chn)
				else: 
					seltext = '(resid %s and segid %s and not name C O N CA HN HA HT1 HT2 HT3)' % (res.rid,res.chn)
			else:
				if not_prot[count]:
					seltext += 'or (resid %s and segid %s and not name C O N CA HN HT1 HT2 HT3) ' % (res.rid,res.chn)
				else:
					seltext += 'or (resid %s and segid %s and not name C O N CA HA HN HT1 HT2 HT3) ' % (res.rid,res.chn)

			if not not_prot[count]:
				captext_mm = 'resid %s and segid %s and name CA' % (res.rid,res.chn)
				captext_qm = 'resid %s and segid %s and name CB' % (res.rid,res.chn)
				captext = [captext_qm,captext_mm]
				caplist.append(captext)
				
				try:
					zerocharge += ' or (resid %s and segid %s and name N HN HT1 HT2 HT3 HA C O CA)' % (res.rid,res.chn)
				except:
					zerocharge = '(resid %s and segid %s and name N HN HT1 HT2 HT3 HA C O CA)' % (res.rid,res.chn)

			count += 1

		return seltext,caplist,zerocharge



def struct_generator(input_psf,input_pdb,seltext,caps,ftype='xyz',ofile='qm_coords',rem=None,tail_text=None,charge=None,mult=None,mm=False,
										    qmmm_mode=None,zerocharge=None,ext_coords=None,wat_thresh=None,
										    wat_layer=None,params=None,bonded=False,read_types=False):

	'''
	This function creates a pdb structure of the qm region specified by seltext
	The caps are used to add hydrogens to satisfy the valency of atoms at the
	QM/MM border

	OUTPUT: xyz coordinate file written to the current directory
	'''

	def h_coords(qm_cap,mm_cap,bond_length=1.10):

		'''
		This function will determine the cartesian coordinates of the hydrogen
		atoms used to cap the QM region of the system
		'''

		qm_coords = u.select_atoms(qm_cap).positions.ravel()

		mm_coords = u.select_atoms(mm_cap).positions.ravel()

		veclength = np.sqrt(np.sum((mm_coords - qm_coords)**2))
		h_xyz = qm_coords + (mm_coords - qm_coords) * (bond_length/veclength)
		
		return h_xyz

	def vdW_params(atype,par_file):
		'''
		This function will read through a list of parameter files to find the
		Lennard Jones van der Waals parameters for the atoms in your system. 

		The current implementation of this function may require a minor adjustment
		of your parameter files: Please be sure to remove any empty lines within
		the NONBONDED section of your CHARMM parameter files, leaving an empty line 
		only at the end of the section.
		'''

		global LJparams
		import re

		try:
			repsilon,rsigma = LJparams[atype]

		except:
			# Check if LJparams dictionary exists
			try: LJparams
			except NameError: LJparams={}

			if type(par_file) != list:
				par_file == list([par_file])

			# Read through the parameter files, looking for the LJ parameter 
			# for the specified atom type
			fcount = 0
			for filename in par_file:
				fileIN = open(filename,'r')
				file_lines = fileIN.readlines()

				activate = False
				for line in file_lines:
					if re.search(r'NONBONDED',line):
						activate = True
						acnt = 0
				
					if activate and acnt > 2:
						values = line.split()

						# Determine if NONBONDED section is done
						if len(values) == 0:
							break

						chars = values[0][0].split()
						if chars[0] != '!':
							epsilon = float(values[2])
							sigma = float(values[3])

							# add the parameters to the dictionary
							dict_val = [epsilon,sigma]
							LJparams[values[0]] = dict_val

							if values[0] == atype:
								repsilon = float(values[2])
								rsigma = float(values[3])

					if activate:
						acnt += 1

				fileIN.close()
				fcount += 1

		try:
			return repsilon,rsigma
		except UnboundLocalError: 
			print('ERROR: Van der Waals parameter for atom type %s could not be found' % atype)
			return None


	def bond_params(par_file,print_all=False,btype=None):
		global bparams
		import re

		try:
			rforce_const,requil_dist = bparams[btype]

		except:
			try:
				rforce_const,requil_dist = bparams[btype.sort(reverse=True)]
			except:
				# Check if bparams dictionary exists
				try: bparams
				except NameError: bparams={}

				if type(par_file) != list:
					par_file == list([par_file])

				for filename in par_file:
					fileIN = open(filename,'r')
					file_lines = fileIN.readlines()

					activate = False
					for line in file_lines:
						if re.search(r'BONDS',line):
							activate = True
							acnt = 0

						if activate and acnt > 2:
							# Read the atom pairs and associated bond parameters
							values = line.split()

							# Determine if the BONDS section is done
							if len(values) == 0:
								break

							chars = values[0][0].split()
							if chars[0] != '!':
								apair = values[0:2]
								par_values = values[2:4]

								if print_all:
									determination = True
									for atype in apair:
										try: type_mapper[atype]
										except KeyError: determination = False
									
									if determination:
										for i in type_mapper[apair[0]]:
											for j in type_mapper[apair[1]]:
												print('Bond %s %s %s %s' % (i,j,par_values[0],
															    par_values[1]),file=fileOUT)
										
								else:
									# add the parameters to the dictionary
									bparams[apair] = par_values
									
									# Check if the atom pair matches the desired bond type
									if apair == btype or apair.sort(reverse=True) == btype:
										rforce_const = par_values[0]
										requil_dist = par_values[1]

						if activate:
							acnt += 1

					fileIN.close()

		if print_all:
			return None
		else:
			return rforce_const,requil_dist


	def angle_params(par_file,print_all=False,theta_type=None):
		global theta_params
		import re

		try:
			rthetak,requil_theta = theta_params[theta_type]
		
		except:
			try:
				rthetak,requil_theta = theta_params[theta_type.sort(reverse=True)]
			except:
				# Check if theta_params dictionary exists
				try: theta_params
				except NameError: theta_params={}

				if type(par_file) != list:
					par_file == list([par_file])

				for filename in par_file:
					fileIN = open(filename,'r')
					file_lines = fileIN.readlines()

					activate = False
					for line in file_lines:
						if re.search(r'ANGLES',line):
							activate = True
							acnt = 0

						if activate and acnt > 2:
							# Read the angle parameters
							values = line.split()
							
							# Determine if ANGLES section is done
							if len(values) == 0:
								break

							chars = values[0][0].split()
							if chars[0] != '!':
								theta_set = values[0:3]
								par_set = values[3:5]
								#par_set[0] = str(float(par_set[0])*57.2958**-2)

								if print_all:
									determination = True
									for atype in theta_set:
										try: type_mapper[atype]
										except KeyError: determination = False
									
									if determination:
										for i in type_mapper[theta_set[0]]:
											for j in type_mapper[theta_set[1]]:
												for k in type_mapper[theta_set[2]]:
													print('Angle %s %s %s %s %s' % (i,j,k,par_set[0],
																	par_set[1]),file=fileOUT)

								else:
									# add the parameters to the dictionary
									theta_params[theta_set] = par_set

									# Check if the atom set matches the desired angle type
									if theta_set == theta_type or theta_set.sort(reverse=True) == theta_type:
										rthetak,requil_theta = par_set
						
						if activate:
							acnt += 1

					fileIN.close()

		if print_all:
			return None
		else:
			return rthetak,requil_theta

		
	def torsion_params(par_file,print_all=False,phipsi_type=None):
		global phipsi_params
		import re

		try:
			rphipsik,rphipsi_ph,rphipsi_n = phipsi_params[phipsi_type]

		except:
			try:
				rphipsik,rphipsi_ph,rphipsi_n = phipsi_params[phipsi_type.sort(reverse=True)]
			except:
				# Check if the phipsi dictionary exists
				try: phipsi_params
				except NameError: phipsi_params={}

				if type(par_file) != list:
					par_file == list([par_file])

				for filename in par_file:
					fileIN = open(filename,'r')
					file_lines = fileIN.readlines()

					activate = False
					for line in file_lines:
						if re.search(r'DIHEDRALS',line):
							activate = True
							acnt = 0

						if activate and acnt > 2:
							# Read the dihedral parameters
							values = line.split()

							# Determine if DIHEDRALS section is done
							if len(values) == 0:
								break

							chars = values[0][0].split()
							if chars[0] != '!':
								phipsi_set = values[0:4]
								par_set = [values[4],values[6],values[5]]

								if print_all:
									determination = True
									for atype in phipsi_set:
										try: type_mapper[atype]
										except KeyError: determination = False

									if determination:
										for i in type_mapper[phipsi_set[0]]:
											for j in type_mapper[phipsi_set[1]]:
												for k in type_mapper[phipsi_set[2]]:
													for l in type_mapper[phipsi_set[3]]:
														print('Torsion %s %s %s %s %s %s %s' % (i,j,k,l,par_set[0],
																			par_set[1],par_set[2]),file=fileOUT)

								else:
									# add the parameters to the dictionary
									phipsi_params[phipsi_set] = par_set

									# Check if the atom set matches the desired dihedral type
									if phipsi_set == phipsi_type or phipsi_set.sort(reverse=True) == phipsi_type:
										rphipsik,rphipsi_ph,rphipsi_n = par_set

						if activate:
							acnt += 1

					fileIN.close()
		
		if print_all:
			return None
		else:
			return rphipsik,rphipsi_ph,rphipsi_n


	def type_reader(filename):
		fileIN = open(filename).readlines()	
		#print(fileIN)
		types = [[] for i in range(len(fileIN))]
		cnt = 0
		for t in fileIN:
			#print(t[:-1])
			types[cnt] = t[:-1]
			cnt += 1
		#print(types)

		return np.array(types)


	def create_qchem_dict(asel,qchem_dict={}):
		count = len(qchem_dict)+1
		for atom in asel.atoms:
			qchem_dict[(atom.name,atom.resid,atom.segid)] = count
			count += 1
			
		return qchem_dict


	def get_h2o(Nwat,qm_seltext,init_thresh=7.0):
		from MDAnalysis.analysis import distances

		init_text = 'name OH2 and around %.1f %s' % (float(init_thresh),qm_seltext)
		thresh_wat = u.select_atoms(init_text)
		index_array = thresh_wat.indices

		# compute the distances between the selected oxygen atoms 
		dist_arr = distances.distance_array(qm_region.positions,
						    thresh_wat.positions,
						    box=u.dimensions)
		min_dists = np.amin(dist_arr,axis=0)	
		sort_instructions = np.argsort(min_dists)
		index_array_sort = index_array[sort_instructions[::-1]]
		
		sel_ix = index_array_sort[-1*int(Nwat):]
		wat_seltext1 = ''
		ix_cnt = 0
		for ix in sel_ix:
			if ix_cnt == 0:
				wat_seltext1 += 'index %s ' % ix
			else:
				wat_seltext1 += 'or index %s ' % ix
			ix_cnt += 1
			
		wat_seltext = ' or ('+wat_seltext1+') '+'or (name H1 H2 and bonded ('+wat_seltext1+'))'
		
		return wat_seltext

	if ftype == 'xyz':
		ofilename = ofile+'.xyz'
		fileOUT = open(ofilename,'w')

	def write_xyz(seltext,h2o_thresh):
		global qm_region

		qm_region = u.select_atoms(seltext)
		
		try:
			h2o_seltext = get_h2o(h2o_thresh,seltext)
			seltext = seltext+h2o_seltext
			qm_region = u.select_atoms(seltext)

		except:
			qm_region = u.select_atoms(seltext)

		qm_names = qm_region.names
		qm_coords = qm_region.positions

		natoms = len(qm_region.atoms) + len(caps)
		if ftype == 'xyz':
			print('%s\n' % natoms,file=fileOUT)

		for i in range(len(qm_region.atoms)):
			atom_name = list(qm_names[i])[0]
			if atom_name == 'F':
				atom_name = 'Fe'
			atom_xyz = qm_coords[i]
			
			print('%s\t%.3f\t%.3f\t%.3f' % (atom_name,atom_xyz[0],atom_xyz[1],atom_xyz[2]),file=fileOUT)

		if len(caps) > 0:
			for cap in caps:
				h_xyz = h_coords(cap[0],cap[1])	
				#print(cap[0],cap[1])
				#print(h_xyz)
				print('H\t%.3f\t%.3f\t%.3f' % (h_xyz[0],h_xyz[1],h_xyz[2]),file=fileOUT)

		return 0

	def xyz_janus(seltext,ext_mol=None,h2o_thresh=None,h2o_layer=None,bonds=False):

		'''

		'''

		# Define the global variables created within this function
		global charge_list; global atom_types
		global tot_qm;      global h2o_seltext
		global qm_region;   global type_dict
		global type_mapper

		# Initialize variables
		charge_list = []
		atom_types = []
		type_dict = dict()
		type_mapper = dict()


		# Partition the QM and MM regions into MDAnalysis atomselections
		# QM First
		qm_region = u.select_atoms(seltext)

		if h2o_layer == 'QM':
			h2o_seltext = get_h2o(h2o_thresh,seltext)
			seltext = seltext+h2o_seltext
			qm_region = u.select_atoms(seltext)
		else:
			qm_region = u.select_atoms(seltext)


		# Now the MM
		# With the MM we must write this section for calculations both with and without a zerocharge selection
		if zerocharge == None:
			if h2o_layer == 'MM':
				mm_seltext = 'not ('+seltext+' or name H1 H2 OH2 SOD CLA)'
				h2o_seltext = get_h2o(h2o_thresh,seltext)
				tot_seltext = mm_seltext+h2o_seltext
			
			else:
				tot_seltext = 'not ('+seltext+')'

		else:
			if h2o_layer == 'MM':
				mm_seltext = 'not ('+seltext+'or '+zerocharge+' or name H1 H2 OH2 SOD CLA)'
				h2o_seltext = get_h2o(h2o_thresh,seltext)
				tot_seltext = mm_seltext+h2o_seltext

			else:
				tot_seltext = 'not ('+seltext+'or '+zerocharge+')'

		mm_region = u.select_atoms(tot_seltext)


		# If bonded parameters are required, generate the dictionary mapping
		# the atom identities to their Q-Chem indices (1-based)
		# Identity has form (name,resid,segid)
		if bonds:
			qchem_dict = create_qchem_dict(qm_region)
			qchem_dict = create_qchem_dict(mm_region,qchem_dict)


		# Write the QM region to the $molecule section 
		qm_names = qm_region.names
		qm_coords = qm_region.positions
		qm_charge = qm_region.charges
		if read_types:
			qm_type = type_reader('qm_types.txt')
			#print(qm_type)
		else:
			qm_type = qm_region.types

		if ext_mol != None:
			# check that the external molecule has the same number of atoms as your qm selection
			len_metric = len(ext_mol) - len(caps)
			if len_metric == len(qm_region.atoms):
				for k in range(len(qm_coords)):
					ext_coords = list(map(float,ext_mol[k][1:]))
					qm_coords[k,:] = ext_coords
			else:
				print('External molecule has a different number of atoms than your qm atom selection.')

		natoms = len(qm_region.atoms) + len(caps)
		count = 1; type_count = 1
		for i in range(len(qm_region.atoms)):
			atom_name = list(qm_names[i])[0]
			if atom_name == 'F':
				atom_name = 'Fe'
			atom_xyz = qm_coords[i]
			atom_charge = qm_charge[i]

			try:
				typeID = type_dict[(qm_type[i],atom_charge)]
				#typeID = type_dict[qm_type[i]]
			except:
				typeID = type_count * -1
				atom_types.append([typeID,(qm_type[i],atom_charge)])
				#type_dict = {atom_types[k][1]:str(atom_types[k][0]) for k in range(len(atom_types))}
				type_dict[(qm_type[i],atom_charge)] = typeID
				type_count += 1
			
			bd_ix = np.zeros(4)
			if bonds:
				# Determine the atom indices for the bonded atoms
				cur_atom = qm_region.atoms[i]
				bdseltext = 'bonded (name %s and resid %s and segid %s)' % (cur_atom.name,cur_atom.resid,cur_atom.segid)
				bdsel = u.select_atoms(bdseltext)

				cnt = 0
				for atom in bdsel.atoms:
					bd_ix[cnt] = qchem_dict[(atom.name,atom.resid,atom.segid)]
					cnt += 1

			print('%s\t%.8f\t%.8f\t%.8f\t%s\t%s\t%s\t%s\t%s' % (atom_name,atom_xyz[0],atom_xyz[1],atom_xyz[2],typeID,
									    int(bd_ix[0]),int(bd_ix[1]),int(bd_ix[2]),int(bd_ix[3])),file=fileOUT)

			charge_list.append(atom_charge)

			count += 1

		tot_qm = (count - 1) + len(caps)
		#print(tot_qm)
		
		# mod 11/13/22 
		if len(caps) > 0:
			typeID = type_count * -1
			charge_list.append(0.09)
			for cap in caps:
				h_xyz = h_coords(cap[0],cap[1])
				if ext_mol != None:
					k += 1
					h_xyz = list(map(float,ext_mol[k][1:]))
				print('H\t%.8f\t%.8f\t%.8f\t%s\t0\t0\t0\t0' % (h_xyz[0],h_xyz[1],h_xyz[2],typeID),file=fileOUT)
			atom_types.append([typeID,('HA3',0.09)])
			#type_dict = {atom_types[k][1]:str(atom_types[k][0]) for k in range(len(atom_types))}     
			type_dict[('HA3',0.09)] = typeID

			count += 1; type_count += 1

		# Write the MM region to the $molecule section
		mm_names = mm_region.names
		mm_coords = mm_region.positions
		mm_charge = mm_region.charges
		mm_type = mm_region.types

		natoms = len(mm_region.atoms)
		for i in range(len(mm_region.atoms)):
			atom_name = list(mm_names[i])[0]
			if mm_names[i] == 'SOD':
				atom_name = 'Na'
			elif atom_name == 'F':
				atom_name = 'Fe'
			elif mm_names[i] == 'CLA':
				atom_name = 'Cl'

			atom_xyz = mm_coords[i]
			atom_charge = mm_charge[i]
			try:
				typeID = type_dict[(mm_type[i],atom_charge)]
			except:
				typeID = type_count * -1
				atom_types.append([typeID,(mm_type[i],atom_charge)])
				#type_dict = {atom_types[k][1]:str(atom_types[k][0]) for k in range(len(atom_types))}     
				type_dict[(mm_type[i],atom_charge)] = typeID
				type_count += 1

			duplicity = float(len(atom_types)) - float(len(type_dict))
			if duplicity != 0.0:
				print(duplicity)

			bd_ix = np.zeros(4)
			if bonds:
				# Determine the atom indices for the bonded atoms (same as we did for the QM region above)
				cur_atom = mm_region.atoms[i]
				bdseltext = 'bonded (name %s and resid %s and segid %s)' % (cur_atom.name,cur_atom.resid,cur_atom.segid)
				bdsel = u.select_atoms(bdseltext)

				cnt = 0
				for atom in bdsel.atoms:
					bd_ix[cnt] = qchem_dict[(atom.name,atom.resid,atom.segid)]
					cnt += 1

			print('%s\t%.8f\t%.8f\t%.8f\t%s\t%s\t%s\t%s\t%s' % (atom_name,atom_xyz[0],atom_xyz[1],atom_xyz[2],typeID,
									    int(bd_ix[0]),int(bd_ix[1]),int(bd_ix[2]),int(bd_ix[3])),file=fileOUT)

			charge_list.append(atom_charge)

			count += 1

		atom_types = {str(atom_types[i][0]):atom_types[i][1] for i in range(len(atom_types))}

		if bonds:
			for i in range(len(atom_types)):
				qtype = -1*(i+1)
				ctype,q = atom_types[str(qtype)]
				try: 
					dict_val = type_mapper[ctype] 
					type_mapper[ctype] = dict_val+[qtype]
				except KeyError: type_mapper[ctype] = [qtype]
		
		return 0	
	
	def write_external_charges(default=True,seltext=None):
		if default:
			mm_seltext = 'not ('+seltext+'or '+zerocharge+')'
		else:
			mm_seltext = seltext
		mm_region = u.select_atoms(mm_seltext)

		mm_coords = mm_region.positions
		mm_charge = mm_region.charges

		natoms = len(mm_region.atoms)
		for i in range(len(mm_region.atoms)):
			atom_xyz = mm_coords[i]
			atom_charge = mm_charge[i]
			
			print('%.3f %.3f %.3f %.3f' % (atom_xyz[0],atom_xyz[1],atom_xyz[2],atom_charge),file=fileOUT)

		# Write zero charge atoms with a charge of 0.000
		#zc_region = u.select_atoms(zerocharge)
		#zc_coords = zc_region.positions
		
		#for i in range(len(zc_region.atoms)):
		#	atom_xyz = zc_coords[i]
		#	zc = 0.000

		#	print('%.3f %.3f %.3f %.3f' % (atom_xyz[0],atom_xyz[1],atom_xyz[2],zc),file=fileOUT)

		return 0
		

	u = mda.Universe(input_psf,input_pdb)

	if ftype == 'xyz':
		ofilename = ofile+'.xyz'
		fileOUT = open(ofilename,'w')

		write_xyz(h2o_thresh=wat_thresh)

	if ftype == 'qchem':
		ofilename = ofile+'.inp'
		fileOUT = open(ofilename,'w')
		
		remfile = open(rem,'r').readlines()
		for line in remfile:
			print(line[:-1],file=fileOUT)

		if qmmm_mode == 'janus':
			print('$molecule\n%s %s' % (charge,mult),file=fileOUT)

			xyz_janus(seltext,ext_coords,h2o_thresh=wat_thresh,h2o_layer=wat_layer,bonds=bonded)

			print('$end',file=fileOUT)

			print('\n$qm_atoms',file=fileOUT)
			print('1:%s' % tot_qm, file=fileOUT)
			print('$end',file=fileOUT)

			# Write the force field parameters required for the atoms in the simulation
			print('\n$force_field_params', file=fileOUT)

			try:
				print('NumAtomTypes %s' % len(atom_types),file=fileOUT)
				for i in range(len(atom_types)):
					idx = -1*(i+1)
					atm_type,chg = atom_types[str(idx)]
					#print(atm_type)
					eps,sig = vdW_params(atm_type,params)
					#print(eps)
					print('AtomType %s %.3f %s %s' % (idx,chg,sig,eps),file=fileOUT)
			except:
				print('NumAtomTypes %s' % len(charge_list),file=fileOUT)
				for i in range(len(charge_list)):
					idx = -1*(i+1)
					print('AtomType %s %.3f 1 0' % (idx,charge_list[i]),file=fileOUT)

			if bonded:
				bond_params(params,print_all=True)
				angle_params(params,print_all=True)
				torsion_params(params,print_all=True)

			print('$end',file=fileOUT)

		#	if wat_thresh != None:
		#		print('\n$external_charges',file=fileOUT)
		#		external_seltext = 'name H1 H2 OH2 and not (name H1 H2 OH2 and around %.1f protein)' % float(wat_thresh)
		#		write_external_charges(default=False,seltext=external_seltext)
		#		print('$end',file=fileOUT)	
		
		else:
			#print(seltext)
			print('$molecule\n%s %s' % (charge,mult),file=fileOUT)
			write_xyz(seltext,h2o_thresh=wat_thresh)	
			print('$end',file=fileOUT)
			
			if mm:
				print('\n$external_charges',file=fileOUT)
				write_external_charges()
				print('$end',file=fileOUT)
		try:
			tailfile = open(tail_text,'r').readlines()
			for line in tailfile:
				print(line[:-1],file=fileOUT)
		except:
			last_action = None

	fileOUT.close()

	return 0

