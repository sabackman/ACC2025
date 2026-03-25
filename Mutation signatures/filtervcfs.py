import os, gzip

name_map = {"Mutect2_filtered_SJ-2335-6031_vs_SJ-2335-6028.vcf.gz":"ACC3_6031",
	"Mutect2_filtered_SJ-2335-10820-2_vs_SJ-2335-1851.vcf.gz":"ACC2_10820-2",
	"Mutect2_filtered_SJ-2335-10817-2_vs_SJ-2335-1851.vcf.gz":"ACC2_10817-2",
	"Mutect2_filtered_SJ-2335-11599_vs_SJ-2335-11600.vcf.gz":"ACC7_11599",
	"Mutect2_filtered_SJ-2335-7504_vs_SJ-2335-1039.vcf.gz":"ACC1_7504",
	"Mutect2_filtered_UI-3077-11912_vs_SJ-2335-2813.vcf.gz":"ACC8_11912",
	"Mutect2_filtered_UI-3077-11914_vs_UI-3077-2625.vcf.gz":"ACC9_11914",
	"Mutect2_filtered_UI-3077-11260_vs_SJ-2335-2414.vcf.gz":"ACC6_11260",
	"Mutect2_filtered_WC-3639-11695_vs_SJ-2335-2813.vcf.gz":"ACC8_11695",
	"Mutect2_filtered_UI-3077-11913_vs_UI-3077-2625.vcf.gz":"ACC9_11913",
	"Mutect2_filtered_UI-3077-A7551-20AT_vs_SJ-2335-2414.vcf.gz":"ACC6_A7551-20AT",
	"Mutect2_filtered_SJ-2335-2177-nr5_vs_SJ-2335-11600.vcf.gz":"ACC7_2177-nr5",
	"Mutect2_filtered_SJ-2335-11259_vs_SJ-2335-2414.vcf.gz":"ACC6_11259",
	"Mutect2_filtered_SJ-2335-11529_vs_SJ-2335-2647.vcf.gz":"ACC5_11529",
	"Mutect2_filtered_SJ-2335-18381-nr3_vs_SJ-2335-2647.vcf.gz":"ACC5_18381-nr3",
	"Mutect2_filtered_UI-3077-11913n2_vs_UI-3077-2625.vcf.gz":"ACC9_11913n2",
	"Mutect2_filtered_SJ-2335-8553-nr4_vs_SJ-2335-2414.vcf.gz":"ACC6_8553-nr4",
	"Mutect2_filtered_SJ-2335-4465_vs_SJ-2335-4471.vcf.gz":"ACC4_4465",
	"Mutect2_filtered_SJ-2335-8749_vs_SJ-2335-1039.vcf.gz":"ACC1_8749",
	"Mutect2_filtered_SJ-2335-5589_vs_SJ-2335-6028.vcf.gz":"ACC3_5589",
	"Mutect2_filtered_UI-3077-A9008-19_vs_SJ-2335-2813.vcf.gz":"ACC8_A9008-19",
	"Mutect2_filtered_UI-3077-11504_vs_UI-3077-2625.vcf.gz":"ACC9_11504",
	"Mutect2_filtered_SJ-2335-8745_vs_SJ-2335-1039.vcf.gz":"ACC1_8745",
	"Mutect2_filtered_UI-3077-7510_vs_SJ-2335-1039.vcf.gz":"ACC1_7510",
	"Mutect2_filtered_UI-3077-11694_vs_SJ-2335-2813.vcf.gz":"ACC8_11694",
	"Mutect2_filtered_SJ-2335-8748_vs_SJ-2335-1039.vcf.gz":"ACC1_8748",
	"Mutect2_filtered_SJ-2335-11691_vs_SJ-2335-2813.vcf.gz":"ACC8_11691",
	"Mutect2_filtered_SJ-2335-11712_vs_SJ-2335-2813.vcf.gz":"ACC8_11712",
	"Mutect2_filtered_SJ-2335-2269_vs_SJ-2335-4471.vcf.gz":"ACC4_2269" }

for f in os.listdir("2025ACCmutpatternvcfs"):
	if not ".vcf.gz" in f: continue
	
	with gzip.open(os.path.join("2025ACCmutpatternvcfs", f), "rt") as infile:
		with open(os.path.join("filteredvcfs", name_map[f]+".vcf"), "w+") as outfile:
			for line in infile:
				if line.startswith("#") or line.split("\t")[6]=="PASS": outfile.write(line)
