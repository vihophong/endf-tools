
import struct
import numpy as np
import re

ensdf = []

findees = ['Y','D','H','M','S','MS','US','NS','PS','FS','AS']
units = {'Y': 365.2422*24*60*60, 'D': 24*60*60, 'H': 60*60, 'M': 60, 'S':1, 'MS': 1e-3, 'US': 1e-6, 'NS': 1e-9, 'PS': 1e-12, 'FS': 1e-15, 'AS': 1e-18}	

elements={"h": 1, "he": 2, "li": 3, "be": 4, "b": 5, "c": 6, "n": 7, "o": 8, "f": 9, "ne": 10, "na": 11, "mg": 12, "al": 13, 
"si": 14, "p": 15, "s": 16, "cl": 17, "ar": 18, "k": 19, "ca": 20, "sc": 21, "ti": 22, "v": 23, "cr": 24, "mn": 25, "fe": 26,
 "co": 27, "ni": 28, "cu": 29, "zn": 30, "ga": 31, "ge": 32, "as": 33, "se": 34, "br": 35, "kr": 36, "rb": 37, "sr": 38, "y": 39,
  "zr": 40, "nb": 41, "mo": 42, "tc": 43, "ru": 44, "rh": 45, "pd": 46, "ag": 47, "cd": 48, "in": 49, "sn": 50, "sb": 51, "te": 52,
   "i": 53, "xe": 54, "cs": 55, "ba": 56, "la": 57, "ce": 58, "pr": 59, "nd": 60, "pm": 61, "sm": 62, "eu": 63, "gd": 64, "tb": 65,
    "dy": 66, "ho": 67, "er": 68, "tm": 69, "yb": 70, "lu": 71, "hf": 72, "ta": 73, "w": 74, "re": 75, "os": 76, "ir": 77, "pt": 78,
     "au": 79, "hg": 80, "tl": 81, "pb": 82, "bi": 83, "po": 84, "at": 85, "rn": 86, "fr": 87, "ra": 88, "ac": 89, "th": 90, "pa": 91,
      "u": 92, "np": 93, "pu": 94, "am": 95, "cm": 96, "bk": 97, "cf": 98, "es": 99, "fm": 100, "md": 101, "no": 102, "lr": 103, "rf": 104,
       "db": 105, "sg": 106, "bh": 107, "hs": 108, "mt": 109, "ds": 110, "rg": 111, "cn": 112, "nh": 113, "fl": 114, "mc": 115, "lv": 116, "ts": 117, "og": 118,"nn" : 400}

def getZ(input):
	"""
	Get atomic number Z by element name
	
	Parameters:
	   input ( str ): Element name
	"""
	if (input==""):
		return -8888
	else:
		sep=re.split('(\d+)',input)
		if sep[2] == "":
			return 500
		if len(sep)==1:
			if sep[2]=="n":
				return int(0)
			elif (sep[2]=="p" or sep[2]=="d" or sep[2]=="t"):
				return int(1)			
			else:
				print "Something wrong! ",input
		else:
			return int(elements[sep[2]])

def getA(input):
	"""
	Get mass number A by element name
	
	Parameters:
	   input ( str ): Element name
	"""
	if (input==""):
		return -9999
	else:
		sep=re.split('(\d+)',input)
		if len(sep)==1:
			if sep[2]=="n":
				return 1
			elif sep[2]=="p":
				return 1
			elif sep[2]=="d":
				return 2
			elif sep[2]=="t":
				return 3
			else:
				print "Something wrong! ",input
		else:
			return int(sep[1])

def geterr(val,err):
	pos = val.rfind(".")
	Epos = val.rfind("E")
	if (Epos != -1):
		ep = float(val[Epos+1:])
		if pos !=-1 :
			errnoep = float(err)*(10**(-(Epos-pos-1)))
			return errnoep*(10**ep)
		else:
			return float(err)*(10**ep)
	else:
		if pos !=-1 :
			return float(err)*(10**(-(len(val)-pos-1)))
		else:
			return float(err)

def getdata(infile):
	file1 = open(infile)
	count = 0
	while True:
		count+=1
		line = file1.readline()
		if not line:
			break
		#print line
		(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16) = struct.unpack("5s1s1s1s1s10s2s18s10s6s9s10s2s1s2s1sx",line)
		if (w3 == " " and w4 == "L" and w5 ==" " and w6.strip() != ""  and w9.strip() != ""):
			w6float = -1
			try:
				w6float = float(w6.strip())
			except:
				w6float = -1
			if (w6float == 0 or w6.strip() == "0.0+X" or w6.strip() == "0+X" or w6.strip() == "0+V" or w6.strip() == "0.0+V" or w6.strip() == "X"):
				isstable = False
				if (w9.strip()=="STABLE"):
					isstable = True
					#print "Stable",w1.strip()
				else:
					if (w10.strip()!=""):
						tmpt12 = w9.strip().split()
						if len(tmpt12) == 2:
							t12str = tmpt12[0]
							t12unitstr = tmpt12[1]
							if (t12unitstr != "MEV" and t12unitstr != "KEV" and t12unitstr != "EV"):
								t12 = float(t12str)
								t12unit = units[t12unitstr]
								t12unc_type = ""
								A = getA(w1.strip().lower())
								Z = getZ(w1.strip().lower())
								nuc = w1.strip().lower()
								dt12p = 0
								dt12m = 0
								if (not (w10.strip().rfind("+") !=-1 and w10.strip().rfind("-") !=-1)):
									if (w10.strip() == "LT" or w10.strip() == "AP" or w10.strip() == "GT" or w10.strip() == "GE"):
										dt12p = 0
										dt12m = 0
										t12unc_type = w10.strip()
									else:
										dt12p = geterr(t12str,w10.strip())
										dt12m = geterr(t12str,w10.strip())
								else:
									if w10.strip().find("+")<w10.strip().find("-"):
										dt12p= geterr(t12str,w10.strip()[w10.strip().find("+")+1:w10.strip().find("-")])
										dt12m= geterr(t12str,w10.strip()[w10.strip().find("-")+1:])
									else:
										dt12m= geterr(t12str,w10.strip()[w10.strip().find("-")+1:w10.strip().find("+")])
										dt12p= geterr(t12str,w10.strip()[w10.strip().find("+")+1:])
								#print A,Z,nuc,t12,dt12p,dt12m,t12unc_type
								ensdf.append({"input": infile, "Z":Z, "A":A, "nuc": nuc, "t12": t12, "dt12p": dt12p, "dt12m": dt12m, "t12unctype": t12unc_type})


		
def get_mult_data(listfile):
	f = open (listfile)
	while True:
		line = f.readline()
		if not line :
			break
		getdata(line.strip())



def writebin(output,data):
	np.save(output,data)

def load_bin(infile):
	return np.load(infile,allow_pickle='TRUE')

def duplicateout():
	curr_data  = []
	for index in range(len(ensdf)):
		flag_fill = True
		for index2 in range(len(curr_data)):
			if (ensdf[index]["Z"] == curr_data[index2]["Z"] and ensdf[index]["A"] == curr_data[index2]["A"]):
				flag_fill = False
		if flag_fill:
			curr_data.append(ensdf[index])

	return curr_data

from ROOT import TTree, TFile, TH2F
def writerootfile(inp,outp):
	x1 = 0; x2 = 240
	y1 = 0; y2 = 115
	nx = x2-x1
	ny = y2-y1
	ensdf_clean = load_bin(inp)
	output_file = TFile.Open(outp, 'recreate')
	h2 = TH2F("halflives","halflives",nx,x1,x2,ny,y1,y2)
	for index in range(len(ensdf_clean)):
		h2.Fill(ensdf_clean[index]["A"]-ensdf_clean[index]["Z"],ensdf_clean[index]["Z"],ensdf_clean[index]["t12"])
		h2.SetBinError(ensdf_clean[index]["A"]-ensdf_clean[index]["Z"]+1,ensdf_clean[index]["Z"]+1,ensdf_clean[index]["dt12p"])
	h2.Write()
	output_file.Close()


# get_mult_data("list1.txt")
# get_mult_data("list2.txt")
# get_mult_data("list3.txt")
# print len(ensdf)
# ensdf_clean = duplicateout()
# print len(ensdf_clean)
# writebin("ensdfdata_t12.npy",ensdf_clean)

writerootfile("ensdfdata_t12.npy","halflives.root")