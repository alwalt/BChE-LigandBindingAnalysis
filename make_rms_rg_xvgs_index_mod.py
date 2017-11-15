#!/usr/bin/python

import os
import commands

def find_xtc_file(working_folder):
    """ we should only find one file """
    n = 0
    for file in os.listdir(working_folder):
        if file[len(file)-4:len(file)] == ".xtc" and 'frame' not in file:
            result = "%s/%s" % (working_folder, file)
            n += 1
    if n != 1:
        result = ""
    return result

def find_tpr_file(working_folder):
    """ we should only find one file """
    n = 0
    for file in os.listdir(working_folder):
        if file[len(file)-4:len(file)] == ".tpr":
            result = "%s/%s" % (working_folder, file)
            n += 1
    if n != 1:
        result = ""
    return result

def find_edr_file(working_folder):
    """ we should only find one file """
    n = 0
    for file in os.listdir(working_folder):
        if file[len(file)-4:len(file)] == ".edr":
            result = "%s/%s" % (working_folder, file)
            n += 1
    if n != 1:
        result = ""
    return result

def find_index_file(working_folder):
    """ we should only find one file """
    n = 0
    for file in os.listdir(working_folder):
        if file[len(file)-4:len(file)] == ".ndx":
            result = "%s/%s" % (working_folder, file)
            n += 1
    if n != 1:
	result = ""
    return result

current_folder = os.getcwd()
folders = os.listdir(current_folder)
folders.sort()
for working_folder in folders:
    if os.path.isdir(working_folder):
            xtc_file = find_xtc_file(working_folder)
            tpr_file = find_tpr_file(working_folder)
	    edr_file = find_edr_file(working_folder)
            index_file = find_index_file(working_folder)
            gyrate_file = "%s/gyrate.xvg" % (working_folder)
            rmsd_file = "%s/rmsd.xvg" % (working_folder)
	    coul_file = "%s/coul.xvg" % (working_folder)
	    lj_file = "%s/lj.xvg" % (working_folder)

#            command = "echo 24 0|trjconv -s %s -f %s -n %s -o %s -center -pbc nojump -ur tric" % (tpr_file, xtc_file, index_file, xtc_file)
#            print command
 #           status, output = commands.getstatusoutput(command)
  #          if status != 0:
   #             print output
    #            sys.exit()

     #       command = "echo 1 0|trjconv -s %s -f %s -n %s -o %s -fit rot+trans" % (tpr_file, xtc_file, index_file, xtc_file)
      #      print command
       #     status, output = commands.getstatusoutput(command)
        #    if status != 0:
         #       print output
          #      sys.exit()

#            command = "echo 24 24|g_rms -s %s -f %s -n %s -o %s" % (tpr_file, xtc_file, index_file, rmsd_file)
 #           print command
  #          status, output = commands.getstatusoutput(command)
   #         if status != 0:
    #            print output
     #           sys.exit()
         
            command = "echo 51 '\n' |g_energy -s %s -f %s -o %s" % (tpr_file, edr_file, lj_file)
            print command
            status, output = commands.getstatusoutput(command)
            if status != 0:
                print output
                sys.exit()

            command = "echo 52 '\n' |g_energy -s %s -f %s -o %s" % (tpr_file, edr_file, coul_file)
            print command
	    status, output = commands.getstatusoutput(command)
            if status != 0:
                print output
                sys.exit()
