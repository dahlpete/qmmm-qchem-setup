#!/Users/peterdahl/opt/anaconda3/bin/python

import sys
import vmdpath2reslist as list_gen
import build_qmmm_struct_v4 as build

path_number = sys.argv[1]
mypsf = 'omcs_au_ox_wb_ions.psf'
mypdb = 'final_conformation_pose1.pdb'

# This line only needs to be run once
list_gen.read_vmdlog('path_identities_pose1.txt',dest='path_lists',heme_substituents=True)

text_filename = 'path_lists/path%s.txt' % path_number
qlist = build.text2list(text_filename)

seltext,caplist = build.selection_text(qlist,mode='bb',cutmode='nc')
ofilename = 'coordinate_files/qm_coords_path%s.xyz' % path_number
build.struct_generator(mypsf,mypdb,seltext,caplist,ofile=ofilename)
