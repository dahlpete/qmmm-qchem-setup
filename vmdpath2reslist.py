import numpy as np
import subprocess
import sys
import re

def heme_sub(resnum,chain,fileOUT):
	class heme1:
		h1 = 13
		h2 = 147
		c1 = 9
		c2 = 12

	class heme2:
		h1 = 2
		h2 = 51
		c1 = 47
		c2 = 50

	class heme3:
		h1 = 144
		h2 = 335
		c1 = 140
		c2 = 143

	class heme4:
		h1 = 110
		h2 = 243
		c1 = 239
		c2 = 242

	class heme5:
		h1 = 16
		h2 = 332
		c1 = 328
		c2 = 331

	class heme6:
		h1 = 257
		h2 = 404
		c1 = 400
		c2 = 403

	def prot_chain(chain,resnum):
		if chain == 'Y':
			pchn = 'A'; npchn = 'B'
		elif chain == 'Z':
			pchn = 'B'; npchn = 'A'
		return pchn,npchn


	pchn,npchn = prot_chain(chain,resnum)

	if int(resnum) == 412:
		h = heme1
	elif int(resnum) == 410:
		h = heme2
	elif int(resnum) == 411:
		h = heme3
	elif int(resnum) == 409:
		h = heme4
	elif int(resnum) == 408:
		h = heme5
	elif int(resnum) == 413:
		h = heme6

	if resnum == 408:
		new_text = '%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s' % ('His',h.h1,npchn,'A','His',h.h2,pchn,'A','Cys',h.c1,pchn,'A','Cys',h.c2,pchn,'A')
	else:
		new_text = '%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s' % ('His',h.h1,pchn,'A','His',h.h2,pchn,'A','Cys',h.c1,pchn,'A','Cys',h.c2,pchn,'A')
	
	print(new_text,file=fileOUT)
	
	

def read_vmdlog(filename,dest,heme_substituents=False):
	fileIN = open(filename,'r').readlines()

	command = 'mkdir %s' % dest
	subprocess.run([command],shell=True)
	
	read_lines = False
	pcount = 0; lcount = 0
	for line in fileIN:
		if re.search(r'PATH 1:',line):
			pcount += 1
			read_lines = True
			lthresh = lcount
			
			outfilename = '%s/path%s.txt' % (dest,pcount)
			fileOUT = open(outfilename,'w')
		
		if read_lines:
			data = line.split()
	
			if lcount > lthresh and len(data) > 1:
				chain = data[1]; restype = data[2]; resnum = data[3]; atmname = data[4]
				print('%s\t%s\t%s\t%s' % (restype,resnum,chain,atmname),file=fileOUT)		
				if (int(resnum) in [412,410,411,409,408,413]):
					heme_sub(resnum,chain,fileOUT)
					
	
			elif lcount > lthresh and len(data) <= 1:
				data = fileIN[lcount-1].split()
				chain = data[7]; restype = data[8]; resnum = data[9]; atmname = data[10]
				print('%s\t%s\t%s\t%s' % (restype,resnum,chain,atmname),file=fileOUT)
				if (int(resnum) in [412,410,411,409,408,413]):
					heme_sub(resnum,chain,fileOUT)
	
				read_lines = False
				fileOUT.close()
	
		lcount += 1
