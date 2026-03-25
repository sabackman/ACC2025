mapping = {"7504":"P1",
	"7510":"P2",
	"8745":"M1",
	"8748":"M2",
	"8749":"M3",
	"10817-2":"P",
	"10820-2":"M",
	"5589":"P",
	"6031":"R",
	"2269":"P",
	"4465":"R",
	"18381-nr3":"P",
	"11529":"R",
	"11259":"P1",
	"11260":"P2",
	"8553-nr4":"M1",
	"A7551-20AT":"M2",
	"2177-nr5":"M",
	"11599":"P",
	"11691":"P1",
	"11694":"P2",
	"WC-3639-11695":"Thrombus",
	"11712":"M1",
	"11912":"M3",
	"A9008-19":"M2",
	"11504":"M",
	"11913":"P1",
	"11913n2":"P2",
	"11914":"R"
}


import os
basedir = "concordance_acc 2"
chromosomes = ["chr"+str(i) for i in range(1, 23)]
for sample in os.listdir(basedir):
	with open(os.path.join(basedir, sample), "r") as infile:
		pairs = []
		line = infile.readline()
		pair = {"name": line.strip()}
		if "vs" in line.strip():
			s1 = mapping[line.strip().split(" vs ")[0]]
			s2 = mapping[line.strip().split(" vs ")[1]]
			pair = {"name": s1+" vs "+s2}
		for line in infile:
			if "chromosome" in line: continue
			if "vs" in line:
				pairs.append(pair)
				s1 = mapping[line.strip().split(" vs ")[0]]
				s2 = mapping[line.strip().split(" vs ")[1]]
				pair = {"name": s1+" vs "+s2}
				continue
			if line.strip() == "": #we're done
				pairs.append(pair)
				continue
			words = line.strip().split()
			pair[words[0]]	= words[-1]
		pairs.append(pair)
	with open(sample.split("_")[0]+"_matrix.txt", "w+") as outfile:
		outline = ["pair"]
		for c in chromosomes:
			outline.append(c)
		outfile.write("\t".join(outline))
		for pair in pairs:
			outline = [pair["name"]]
			for c in chromosomes:
				outline.append(pair[c])
			outfile.write("\n"+"\t".join(outline))
			
