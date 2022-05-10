#!/gpfs/loomis/project/batista/pd455/conda_envs/mdanalysis/bin/python

import sys
import build_qmmm_struct_v5 as build
import subprocess

#hnum = sys.argv[1]
base = sys.argv[1]
n_snaps = int(sys.argv[2])

mypsf = '../../snapshots/omcs_2chains_ox_wb_ionized.psf'
text_filename = 'qm.txt'

print('Setting up Q-Chem calculations')

for i in range(n_snaps):
	mypdb = '../../snapshots/%s_%s.pdb' % (base,i)

	for state in ['oxidized','reduced']:
		if state == 'oxidized':
			qm_charge = -1
			qm_multiplicity = 2
		elif state == 'reduced':
			qm_charge = -2
			qm_multiplicity = 1

		if i == 0:
			shellcommand = "mkdir %s" % (state)
			subprocess.run([shellcommand],shell=True)
			
		
		ofilename = '%s/%s_on_%s_%s' % (state,state,base,i)
		
		qlist = build.text2list(text_filename)
		seltext,caplist,zcharge = build.selection_text(qlist,mode='sc')

		build.struct_generator(mypsf,
				       mypdb,
				       seltext,
				       caplist,
				       ofile = ofilename,
				       ftype='qchem',
				       rem='myrem.rem',
				       charge=qm_charge,
				       mult=qm_multiplicity,
				       mm=True,
				       zerocharge=zcharge)
