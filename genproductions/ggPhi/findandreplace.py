file="unweighted_events_100000_mass6.lhe"
fo=open("unweighted_events_6000211_100000events_mass6.lhe","w")
with open(file,"r") as f:
    for line in f:
        newline = line
        spline=line.split()
        if(len(spline)==13):
            if(spline[0]=="25" and spline[1]=="1" and spline[4] == "0"): #and spline[2]=="1"):
#                print(spline)
                #spline[0] = "6000211"
                #newline = line[1*(len(spline[0])+1):] + '\n'
                newline = line.replace("     25", "6000211")
#                print(newline)
                
        fo.write(newline)
fo.close


