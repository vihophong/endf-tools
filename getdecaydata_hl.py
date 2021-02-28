
import struct
import numpy as np

endf = []

def gettype(input):
	if input == 0:
		return "g"
	elif input == 1:
		return "b"
	elif input == 2:
		return "b_p"
	elif input == 4:
		return "a"
	elif input == 5:
		return "n"
	elif input == 6:
		return "sf"
	elif input == 7:
		return "p"
	elif input == 8:
		return "e_m"
	elif input == 9:
		return "x"
	elif input == 10:
		return "a_nue"
	elif input == 11:
		return "nue"
	else:
		return "n/a"

def getval(input):
	expsign=input.rfind("E")
	if (expsign != -1):
		right = int(input[expsign+1:])
		left = float(input[:expsign])
		return left*(10**(right))

	plus = input.rfind("+")
	minus = input.rfind("-")
	right = -9999
	left = -9999
	if (plus != -1 and plus!=0):
		right = int(input[plus:])
		left = float(input[:plus])
	if (minus != -1 and minus!=0):
		right = int(input[minus:])
		left = float(input[:minus])
	if (right == -9999 or left == -9999):
		return float(input)
	return left*(10**(right))
	
def getdata(infile):
	file1 = open(infile)
	count = 0


	#general header
	za = 0; awr = 0; lis = 0; liso = 0; nst = -1; nsp = 0; t12 = 0
	dt12 = 0; nc = 0; ex1 = 0; dex1 = 0; ex2 = 0; dex2 = 0; ex3 = 0; dex3 = 0
	spi = 0; par = 0; ndk = -1
	nstart = -1

	specall = []


	nnext = -1
	ncurr = -1

	flag_spec = False
	flag_read = False
	nstart = -1

	ispec = 0 
	ival = 0

	stype = ""

	rad = []

	spec = []


	while True:
		count+=1
		line = file1.readline()
		#break
		if not line :
			break
		(hl,mat,mf,mt,ns) = struct.unpack("66s4s2s3s5sxx",line)
		(mat,mf,mt,ns) = map(int,(mat.strip(),mf.strip(),mt.strip(),ns.strip()))
		if (mt==457):
			(word1,word2,word3,word4,word5,word6) = struct.unpack("11s11s11s11s11s11s",hl)

			if ns == 1:
				za = getval(word1.strip())
				awr = getval(word2.strip())#nuclear mass ratio to neutron
				lis = int(word3.strip())#excited state
				liso = int(word4.strip())#isomeric state
				nst = int(word5.strip())#stablility flag
				nsp = int(word6.strip())# number of radiation type
			if (nst==0):#if radioactive nucleus
				# SECTION HEADER
				if ns == 2:
					t12 = getval(word1.strip())#half life
					dt12 = getval(word2.strip())#d half life
					nc = int(word5.strip())/2#number of decay energy given
					#if nc != 3:
						#print (za)
				if ns == 3:
					(ex1,dex1,ex2,dex2,ex3,dex3) = map(lambda x: getval(x),(word1.strip(),word2.strip(),word3.strip(),word4.strip(),word5.strip(),word6.strip()))
				if ns == 4:
					spi = getval(word1.strip())
					par = getval(word2.strip())
					ndk = int(word6.strip())
				if (ns >= 5 and ns < 5+ndk):
					rad.append({"rtyp": getval(word1.strip()),"rfs": getval(word2.strip()),"q": getval(word3.strip()),
						"dq": getval(word4.strip()),"br": getval(word5.strip()),"dbr": getval(word6.strip())})

	Z = int(round(za))/1000
	A = int(round(za))-Z*1000
	endf.append({"input": infile, "Z":Z, "A":A, "za": za, "awr": awr, "lis": lis, "liso": liso, "nst": nst,
				 "nsp": nsp,"t12": t12,"dt12": dt12,"nc": nc,"ex1": ex1,"dex1": dex1,"ex2": dex2,"ex3": dex3,
				 "spi": spi,"par": par,"ndk": ndk,"radia": rad})

#getdata("data/decay_1565_50-Sn-134.dat")
#print endf[0]["specs"][0]["stype"]

def getneuspec(listfiles):
	f = open (listfiles)
	while True:
		line = f.readline()
		if not line :
			break
		getdata(line.strip())
		#print line.strip()
	
	endf_neuspecs = []
	for index in range(len(endf)):
		if endf[index]["lis"]==0 and endf[index]["liso"]==0:#only ground state
			endf_neuspecs.append(endf[index])
	return endf_neuspecs

def writebin(input, output):
	endf_neuspecs = getneuspec(input)
	print len(endf_neuspecs),"entries saved"
	np.save(output,endf_neuspecs)

def load_bin(infile):
	return np.load(infile,allow_pickle='TRUE')


#writebin("listfiles.txt","all_beta_specs_hlonly")
#endf_neuspecs = load_bin("all_beta_specs.npy")

#endf_neuspecs = getneuspec("listfiles.txt")
#count = 0
# for index in range(len(endf_neuspecs)):
# 	if (endf_neuspecs[index]["nst"]==0 and endf_neuspecs[index]["dt12"]>0):
# 		print endf_neuspecs[index]["Z"],endf_neuspecs[index]["A"],endf_neuspecs[index]["t12"],endf_neuspecs[index]["dt12"]
# 		count+=1

# print count,len(endf_neuspecs)

endf_neuspecs = load_bin("all_beta_specs_hlonly.npy")

from ROOT import TTree, TFile, TH2F
x1 = 0; x2 = 240
y1 = 0; y2 = 115
nx = x2-x1
ny = y2-y1

output_file = TFile.Open("halflives.root", 'recreate')
h2 = TH2F("halflives","halflives",nx,x1,x2,ny,y1,y2)
for index in range(len(endf_neuspecs)):
	if (endf_neuspecs[index]["nst"]==0 and endf_neuspecs[index]["dt12"]>0):
		h2.Fill(endf_neuspecs[index]["A"]-endf_neuspecs[index]["Z"],endf_neuspecs[index]["Z"],endf_neuspecs[index]["t12"])
		h2.SetBinError(endf_neuspecs[index]["A"]-endf_neuspecs[index]["Z"]+1,endf_neuspecs[index]["Z"]+1,endf_neuspecs[index]["dt12"])

h2.Write()
output_file.Close()