#!usr/bin/env python3
"""
pimp_my_network.py:

python3 pimp_my_network.py 

"""
## IMPORTS

from sys import argv

import re

## functions

def split(word): 
    return [char for char in word]

## main

if __name__ == "__main__":

    if len(argv) == 2:

        input_names = argv[1]
      
        infile1 = open(input_names)

        while True:

            line1 = infile1.readline()
   
            if len(line1) > 0:

                file_name = ("{}".format(line1.split("\n")[0]))
        
                infile2= open(file_name)

                outfile_name = ("{}".format(("out" + file_name)))
                outfile = open(outfile_name, "w")
        
                while True:

                    line2 = infile2.readline()

                    if len(line2) > 0:

                         test = line2.split(" ")               

                         if split(test[0])[1] == "S":
                     
                             col_rep3 = ("SSU")

                         elif split(test[0])[1] == "L":

                             col_rep3 = ("LSU")

                         else:
                             col_rep3 = "NA"


                         if split(test[1])[1] == "S":

                             col_rep4 = ("SSU")

                         elif split(test[1])[1] == "L":

                             col_rep4 = ("LSU")

                         else:
                             col_rep4 = "NA"

                         outfile.write(test[0] + " " + test[1] + " " + test[2].split("\n")[0] + " " + col_rep3 + " " + col_rep4 + "\n")

                    else:
                        break

                infile2.close()
                outfile.close()

            else:
                break


        infile1.close()

        print("done")

    elif len(argv) == 4:
        input_names = argv[1]
        input_IntCryOmics = argv[2]
        protein_name = argv[3]

        infile1 = open(input_names)
        infile3 = open(input_IntCryOmics)

        count = 0

        while True:

            count = count + 1

            #print(count)

            line3 = infile3.readline()
        
            if len(line3) > 0:
                
                print(line3)

                file_name3 = ("{}".format(line3.split("\n")[0]))

                infile4 = open(file_name3)

                count2 = 0

                Regions = str()

                while True:

                    count2 = count2 + 1

                    #print(count2)

                    line4 = infile4.readline()    

                    if len(line4) > 0:

                        if re.search("Region", line4):

                            if re.search(("{}".format(protein_name)), line4):

                                Region = str()

                                for i in range(len(line4.split("'"))):

                                    test1 = ((2*i)+1)

                                    if test1 < len(line4.split("'")):

                                        runner = line4.split("'")[test1]

                                    else:
                                        runner = str()

                                    Region = Region + " " + runner

                                    #print(Region)

                                Regions = Regions + " " + Region  
                              
                    else:
                        break
                        
                ## remove empty spaces from regions

                Regions_list = list(Regions.split(" "))

                while("" in Regions_list):
                    Regions_list.remove("")

            ## mapping region belingings

            line1 = infile1.readline()
   
            if len(line1) > 0:

                file_name = ("{}".format(line1.split("\n")[0]))
        
                infile2= open(file_name)

                outfile_name = ("{}".format(("out2" + file_name)))
                outfile = open(outfile_name, "w")
        
                while True:

                    line2 = infile2.readline()

                    if len(line2) > 0:

                        test = line2.split(" ")                         

                        #print(test[0] in Regions_list)

                        if test[0] in Regions_list:

                            col_rep3 = (protein_name)

                        elif split(test[0])[1] == "S":
        
                            col_rep3 = ("SSU")

                        elif split(test[0])[1] == "L":

                            col_rep3 = ("LSU")

                        else:
                            col_rep3 = "NA"

                        if test[1] in Regions_list:

                            col_rep4 = (protein_name)

                        elif split(test[1])[1] == "S":

                            col_rep4 = ("SSU")

                        elif split(test[1])[1] == "L":

                            col_rep4 = ("LSU")

                        else:
                            col_rep4 = "NA"

                        outfile.write(test[0] + " " + test[1] + " " + test[2].split("\n")[0] + " " + col_rep3 + " " + col_rep4 + "\n")

                    else:
                        break

                infile4.close()
                infile2.close()
                outfile.close()

            else:
                 break


        infile1.close()
        infile3.close()

    else:
        print('################################################################################################\n')
        print('    Usage: python3 pimp_my_network.py [Names_file] [Intcryomics_file] [protein_ID]    \n')
        print('################################################################################################\n')


