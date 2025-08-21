import os
from os import path

# list the directory contents
#directory = "U:\\PLAN\\AMC\\GC-Agilent+Markes\\Data\\MDVR\\2021\\202105 MDL\\"
#directory = "U:\\PLAN\\Rkaullman\\GC Data Exp\\"
#directory = "D:\\Testing\\"
directory = "C:\\Users\\rkaullman\\AutoGC Data\\RB\\Result\\RB\\"
if (path.exists(directory)):
    files = os.listdir(directory)

    
    for filename in files:
        
        if filename.endswith(".dat.tx1"):
            newname = filename.replace(".dat", "", 1)
            print(filename + " => " + newname)
            # Change hwqe20h.dattx1 to hwqe20h.tx1
            os.rename(directory + filename, directory + newname)
            
