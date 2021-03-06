
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

				# SUB SESSION
				if (ndk >= 0 and ns >= 5+ndk):
					if (word1.strip()=="0.000000+0" and word5.strip()=="6"):
						flag_spec = True
						ispec += 1
						if (ispec!=1):
							specall.append({"stype": stype, "spec":spec})

						stype = gettype(int(getval(word2.strip())))
						spec = []

					if (flag_spec):
						if (word4.strip()=="" and word5.strip()=="" and word6.strip()==""):
							#print hl
							flag_spec = False
							nstart = ns+1
							ival = 0					

					if (nstart>0 and  (not flag_spec) and (not (word4.strip()=="" and word5.strip()=="" and word6.strip()==""))):
						#print ispec-1,getval(word1.strip())
						if word1.strip() != "":
							val1 = getval(word1.strip())
							dval1 = getval(word2.strip())
							spec.append({"val": val1, "dval": dval1})
							ival+=1
						if word3.strip() != "":
							val2 = getval(word3.strip())
							dval2 = getval(word4.strip())
							spec.append({"val": val2, "dval": dval2})
							ival+=1
						if word5.strip() != "":
							val3 = getval(word5.strip())
							dval3 = getval(word6.strip())
							spec.append({"val": val3, "dval": dval3})
							ival+=1

	if ispec > 0:
		specall.append({"stype": stype,"spec":spec})
		Z = int(round(za))/1000
		A = int(round(za))-Z*1000
		endf.append({"input": infile, "Z":Z, "A":A, "za": za, "awr": awr, "lis": lis, "liso": liso, "nst": nst,
				 "nsp": nsp,"t12": t12,"dt12": dt12,"nc": nc,"ex1": ex1,"dex1": dex1,"ex2": dex2,"ex3": dex3,
				 "spi": spi,"par": par,"ndk": ndk,"radia": rad,"qbn": [-1,-1,-1,-1,-1], "dqbn": [-1,-1,-1,-1,-1],"pn": [0,0,0,0,0],"dpn": [0,0,0,0,0],
				 "specs": specall})

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
			for ip in range(len(endf[index]["specs"])):
				if endf[index]["specs"][ip]["stype"] == "n":
					if endf[index]["lis"]==0 and endf[index]["liso"]==0:#only ground state
						flag_bdn = False
						for irad in range(len(endf[index]["radia"])):
							if endf[index]["radia"][irad]["rtyp"] == 1.0:
								endf[index]["qbn"][0] = endf[index]["radia"][irad]["q"]
								endf[index]["dqbn"][0] = endf[index]["radia"][irad]["dq"]
								endf[index]["pn"][0] = endf[index]["radia"][irad]["br"]
								endf[index]["dpn"][0] = endf[index]["radia"][irad]["dbr"]
								flag_bdn = True
							if endf[index]["radia"][irad]["rtyp"] == 1.5:
								endf[index]["qbn"][1] = endf[index]["radia"][irad]["q"]
								endf[index]["dqbn"][1] = endf[index]["radia"][irad]["dq"]
								endf[index]["pn"][1] = endf[index]["radia"][irad]["br"]
								endf[index]["dpn"][1] = endf[index]["radia"][irad]["dbr"]
								flag_bdn = True
							if endf[index]["radia"][irad]["rtyp"] == 1.55:
								endf[index]["qbn"][2] = endf[index]["radia"][irad]["q"]
								endf[index]["dqbn"][2] = endf[index]["radia"][irad]["dq"]
								endf[index]["pn"][2] = endf[index]["radia"][irad]["br"]
								endf[index]["dpn"][2] = endf[index]["radia"][irad]["dbr"]
							if endf[index]["radia"][irad]["rtyp"] == 1.555:
								endf[index]["qbn"][3] = endf[index]["radia"][irad]["q"]
								endf[index]["dqbn"][3] = endf[index]["radia"][irad]["dq"]
								endf[index]["pn"][3] = endf[index]["radia"][irad]["br"]
								endf[index]["dpn"][3] = endf[index]["radia"][irad]["dbr"]
							if endf[index]["radia"][irad]["rtyp"] == 1.5555:
								endf[index]["qbn"][4] = endf[index]["radia"][irad]["q"]
								endf[index]["dqbn"][4] = endf[index]["radia"][irad]["dq"]
								endf[index]["pn"][4] = endf[index]["radia"][irad]["br"]
								endf[index]["dpn"][4] = endf[index]["radia"][irad]["dbr"]
						if flag_bdn:
							endf_neuspecs.append(endf[index])
	return endf_neuspecs

def writebin(input, output):
	endf_neuspecs = getneuspec(input)
	print len(endf_neuspecs),"entries saved"
	np.save(output,endf_neuspecs)

def load_bin(infile):
	return np.load(infile,allow_pickle='TRUE')

elements={"h": 1, "he": 2, "li": 3, "be": 4, "b": 5, "c": 6, "n": 7, "o": 8, "f": 9, "ne": 10, "na": 11, "mg": 12, "al": 13, 
"si": 14, "p": 15, "s": 16, "cl": 17, "ar": 18, "k": 19, "ca": 20, "sc": 21, "ti": 22, "v": 23, "cr": 24, "mn": 25, "fe": 26,
 "co": 27, "ni": 28, "cu": 29, "zn": 30, "ga": 31, "ge": 32, "as": 33, "se": 34, "br": 35, "kr": 36, "rb": 37, "sr": 38, "y": 39,
  "zr": 40, "nb": 41, "mo": 42, "tc": 43, "ru": 44, "rh": 45, "pd": 46, "ag": 47, "cd": 48, "in": 49, "sn": 50, "sb": 51, "te": 52,
   "i": 53, "xe": 54, "cs": 55, "ba": 56, "la": 57, "ce": 58, "pr": 59, "nd": 60, "pm": 61, "sm": 62, "eu": 63, "gd": 64, "tb": 65,
    "dy": 66, "ho": 67, "er": 68, "tm": 69, "yb": 70, "lu": 71, "hf": 72, "ta": 73, "w": 74, "re": 75, "os": 76, "ir": 77, "pt": 78,
     "au": 79, "hg": 80, "tl": 81, "pb": 82, "bi": 83, "po": 84, "at": 85, "rn": 86, "fr": 87, "ra": 88, "ac": 89, "th": 90, "pa": 91,
      "u": 92, "np": 93, "pu": 94, "am": 95, "cm": 96, "bk": 97, "cf": 98, "es": 99, "fm": 100, "md": 101, "no": 102, "lr": 103, "rf": 104,
       "db": 105, "sg": 106, "bh": 107, "hs": 108, "mt": 109, "ds": 110, "rg": 111, "cn": 112, "nh": 113, "fl": 114, "mc": 115, "lv": 116, "ts": 117, "og": 118,
       "al-":119, "al*": 120}

def getnamebyz(z):
	"""
	Get element name by atomic number Z
	
	Parameters:
	   z ( int ): Atomic number Z
	"""
	keys=elements.keys()
	values=elements.values()
	return keys[values.index(z)]


#writebin("listfiles.txt","all_beta_specs")
endf_neuspecs = load_bin("all_beta_specs.npy")

from ROOT import TTree, TFile, TH1F

for index in range(len(endf_neuspecs)):
			for ip in range(len(endf_neuspecs[index]["specs"])):
				if endf_neuspecs[index]["specs"][ip]["stype"] == "n":
					val=np.zeros(len(endf_neuspecs[index]["specs"][ip]["spec"]))
					dval=np.zeros(len(endf_neuspecs[index]["specs"][ip]["spec"]))
					for ips in range(len(endf_neuspecs[index]["specs"][ip]["spec"])):
						val[ips]=(float(endf_neuspecs[index]["specs"][ip]["spec"][ips]["val"]))/1000000
						dval[ips]=(float(endf_neuspecs[index]["specs"][ip]["spec"][ips]["dval"]))
						
					if len(endf_neuspecs[index]["specs"][ip]["spec"])>0:
						filename="specs/"+getnamebyz(int(endf_neuspecs[index]["Z"])).capitalize()+str(endf_neuspecs[index]["A"])+".root"
						output_file = TFile.Open(filename, 'recreate')
						h1 = TH1F("hSpecRebin","hSpecRebin",len(endf_neuspecs[index]["specs"][ip]["spec"])-1,val)
						for ii in range(len(endf_neuspecs[index]["specs"][ip]["spec"])-1):
							h1.SetBinContent(ii+1,dval[ii])
						h1.Scale(1/h1.Integral())
						h1.Write()
						output_file.Close()
						print "./neueff_from_spec.sh endf-tools/"+filename+" hSpecRebin upc_brikenV69_wClover.txt "+str(endf_neuspecs[index]["qbn"][1]/1000000)+" "+str(endf_neuspecs[index]["qbn"][0]/1000000)
