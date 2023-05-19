#!{path2python}/python

import sys
import build_qmmm_struct_v11 as build
import subprocess

base = sys.argv[1]
nsnaps = int(sys.argv[2])

mypsf = '{psf_filename}.psf'
text_filename = 'qm.txt' 

print('Setting up Q-Chem calculations')

	
qlist = build.text2list(text_filename)
my_seltext,caplist,zcharge = build.selection_text(qlist,mode='sc')

state = 'qmmm_dir'

for i in range(nsnaps):
	mypdb = '{path2pdbfiles}/%s_%s.pdb' % (base,i)

	if i == 0:
		shellcommand = "mkdir %s" % (state)
		subprocess.run([shellcommand],shell=True)

	ofilename = '%s/%s_on_%s_%s' % (state,state,base,i)

	build.struct_generator(input_psf  = mypsf,
			       input_pdb  = mypdb,
			       seltext    = my_seltext,
			       caps       = caplist,
			       ofile      = ofilename,
			       ftype      = 'qchem',
			       rem        = 'myrem.rem',
			       charge     = 0,
			       mult       = 1,
			       mm         = True,
			       qmmm_mode  = 'janus',
			       zerocharge = zcharge,
			       params     = ['{path2params}/par_all36_prot.prm',
					     '{path2params}/toppar_water_ions_modified.prm'])
