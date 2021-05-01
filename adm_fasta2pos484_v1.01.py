#!/usr/bin/env python3

import sys
import datetime
import re

start = datetime.datetime.now()
timer0 = str(start)
sys.stderr.write("Start time:			%s\n\nPython computations running...\n" % timer0)

inf = sys.argv[1]
outf = open(inf[:-6] + "_single.fasta", "w+")

with open(inf) as f:
	i = 0
	header = ''
	seq = ''
	newline= 0
	for line in f:
		if line.startswith('>'):
			i += 1
			outf.write(seq)
			if newline == 1:
	    			outf.write('\n')
			newline = 1
			header = line.strip('\n') + '\n'
			outf.write(header)
			seq = ''
		else:
			seq += line.strip('\n')
#		if i == 10:
#			break
outf.write(seq)
outf.close()

outf2 = open(inf[:-6] + "_pos484-res.tsv", "w+")
with open(inf[:-6] + "_single.fasta") as f2:
	outf2.write('lab_id\tres_484')
	header = "no"
	num_a = 0
	num_c = 0
	num_g = 0
	num_t = 0
	for line in f2:
		if line.startswith('>'):
			header = "yes"
			fasta_title = line[14:21]
		else:
			header = "no"
		if header == "no":
			mut = ''	
			vr_string = ''	
			q =   'ACACCTTGTAATGGTGTT'
			m = re.search(r'(?<=' + q + ')\w+', line)
			t = 'G'
			l = len(t)
			if m:
				tempseq = m.group(0)
				VR = tempseq[0:l]
				#Removing consensus residues			
				n = 0
				for res in VR:
					if res == t[n]:
						vr_string += '.'
						num_g += 1
						mut = 'E484'
						#print(line)
					else:
						vr_string += res
						if res == "C":
							num_c += 1
							mut = 'E484Q'
							#print(line)
						if res == "A":
							num_a += 1
							#print(line[position -10: position +2])
							mut = 'E484K'
						if res == "T":
							num_t += 1
							mut = 'STOP'
							#print(line)
					n = n + 1
				outf2.write('\n' + fasta_title + '\t' + mut)
				#print(vr_string)
			
	print('Number of A:	' + str(num_a))
	print('Number of C:	' + str(num_c))
	print('Number of G:	' + str(num_g))
	print('Number of T:	' + str(num_t))
outf2.close()

end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("Python Completed:		%s\n\n" % timer_0)

"""
PURPOSE:

PROCEDURE:
1.
2.
3.
4.




Sample input file:


Sample output file:

"""
