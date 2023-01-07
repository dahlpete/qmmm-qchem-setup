# qmmm-qchem-setup
Collection of scripts for the automated setup of QM/MM simulations using Q-Chem.


The scripts in this directory allow for generation of .xyz coordinate files
or Q-Chem input files corresponding to a user defined QM region. The program 
requires a text file containing a list of atoms formatted as follows:

  resname       resid   chain   atom_name


An accompanying script is included which allows for the generation of QM
coordinate files from the residues identified in electron transfer path
calculations, using the VMD Pathways plugin. This script explicitly allows
for the inclusion of heme c axial ligands from OmcS (PDB ID: 6ef8). Their
inclusion can be toggled via the "heme_substituents" variable.

The run_qmmm_setup.py script then sequentially calls three functions from
the build_qmmm_struct_v{version}.py script:

         The first (text2list) generates a list of python class objects
whose attributes correspond to the atom information in the input text file.

        The second (selection_text) takes this list as input and generates
two outputs: an MDAnalysis selection text string (atom names use the CHARMM
naming scheme) and a list of paired MDAnalysis selection text strings.
The list of paired strings define the caps or boundaries between the QM and
MM regions of the model. This function allows for inclusion as well as
exclusion of peptide backbone atoms. If backbone atoms are included (toggled
via the "mode" variable which takes values of bb or sc) two modes
of partitioning are available. The first is C-C, in reference to defining
the boundary between the alpha carbon and carbonyl carbons of the peptide
backbone. The second is N-C, in reference to defining the boundary at the
peptide bond.

        The third (struct_generator) takes the psf and pdb files as well as
the selection text output to generate an xyz coordinate file for the QM
region. This function uses the MDAnalysis module, so make sure you have that
downloaded.

          
VMD Pathways to QM USAGE: ./run_path2qmmm.py {path_num}
          
        path_num refers to the n-th path in the vmd log file
        

QChem QMMM setup USAGE: ./run_qmmm_setup.py {textfile} {basename} {num snapshots}
      
        The run_qmmm_setup.py script needs to be updated with the user's
        file paths. 
        
        For the qchem example and qchem janus example, use omcs_ox_snap for the basename
    
** Note: The header of the run scripts should be updated with the path to
         the user's version of Python3 (must have MDAnalysis downloaded)


VERSION HISTORY

Version 5 of build_qmmm_struct introduced the inclusion of the MM atoms via
the electrostatic embedding method when using 'sc' mode. Future versions of
this program will also allow for electrostatic embedding when in 'bb' mode.

Version 8 


Version 9


Version 10 of build_qmmm_struct introduced the ability to include specific water 
molecules into either the QM layer or the MM layer of the calculation. Two modes
are available. You can either include water based on a distance threshold or you 
can include a specific number of waters, selected by proximity to the QM selection.
When in the QM layer, the rest of the waters are included in the MM layer. When in 
the MM layer, the remainder of the waters are excluded.


Version 11 of build_qmmm_struct introduced the ability to include bonded 
parameters in a Janus QM/MM calculation.
