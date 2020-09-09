file="filelist_hzd.txt"
fo=open("hzd_m0p25_ct100.txt","w")

with open(file,"r") as f:
    for line in f:
        newline = line
        if "mzd0p25_ctau100mm" in line:
            fo.write(newline)
fo.close
