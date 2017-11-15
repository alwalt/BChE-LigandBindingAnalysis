#!/usr/bin/python

import matplotlib
matplotlib.use( 'Agg' )
import pylab
import subprocess
import os

PWD = subprocess.Popen('pwd', stdout = subprocess.PIPE)  # create a subprocess called PWD
WorkingDir = PWD.stdout.read() # Run the command
WorkingDir = WorkingDir.rstrip('\r\n') # remove endline characters from the output string
PlotDir = WorkingDir.split("choline_set_1") # SHOULD SAY APO-data, get rid of extraneous foldernames
# print PlotDir
Plot_title = PlotDir[1]
print Plot_title # something like /1LS4/300/A94/rf/L1.2-C1.2-V1.2-S1.0/1LS4-A94/
PlotFolders = Plot_title.split("/") # split on each / to separate folder names
s = "_"
FileName = s.join( PlotFolders ) # joins the items together w/ specified char
# following subprocess finds all rmsd.xvg files in current folder tree
findrmsd = subprocess.Popen("locate rmsd.xvg | grep %s" % WorkingDir, stdout = subprocess.PIPE, shell = True)  
rmsdfilepaths = findrmsd.stdout.read() # makes a list of all rmsd files found
findrg = subprocess.Popen("locate gyrate.xvg | grep %s" % WorkingDir, stdout = subprocess.PIPE, shell = True)
rgfilepaths = findrg.stdout.read() # lists all gyrate.xvg

if rmsdfilepaths and rgfilepaths: # quick loop to ensure that the search was successful
	print "The rmsdfilepaths and rgfilepaths were successful!" 
	# output is something like ~/001/rmsd.xvg, ~/002/rmsd.xvg, ~/003/rmsd.xvg, etc.
else:
	print 'No .xvg files were found...\n' ### if the .xvg files were recently generated, try running su updatedb FIRST! ### 
	exit()
rmsdfiles = rmsdfilepaths.split("\n") # make each filepath a separate item
rgfiles = rgfilepaths.split("\n")

finaltime = '0'
finalrg = '0'
for rmsdfile in rmsdfiles:
	if "rmsd" in rmsdfile:
		#print rmsdfile # it checks out	
		rmsdtail = subprocess.Popen('tail %s -n 1' % rmsdfile, stdout = subprocess.PIPE, shell = True)
		rmsdline = rmsdtail.stdout.read() # looks like 'time    rmsdvalue\n'
		rmsdline = rmsdline.rstrip('\r\n')
		rmsdsplit = rmsdline.split("   ") # "[time, rmsdvalue]"
		if "rmsdsplit[0]" > finaltime:
			finaltime = rmsdsplit[0] # will eventually give the longest runtime
			#print finaltime, "is finaltime"
	else:
		rmsdfiles.remove(rmsdfile)# removes non-empty useless items in list (this was an issue)
for rgfile in rgfiles:
	if "gyrate" in rgfile:
		#print rgfile # it checks out	
		rgtail = subprocess.Popen('tail %s -n 1' % rgfile, stdout = subprocess.PIPE, shell = True)
		rgline = rgtail.stdout.read() # looks like 'time    rgvalue\n'
		rgline = rgline.rstrip('\r\n')
		rgsplit = rgline.split(" ") # "[time, rgx, rgy, rgz]"
		for rgsplitvalues in rgsplit:	# only problem is the runtime is in that line and that is always bigger than rg
			if "." in rgsplitvalues:				
				if rgsplitvalues > finalrg:
					finalrg = rgsplitvalues # will eventually give the longest runtime
	else:
		rgfiles.remove(rgfile)# removes non-empty useless items in list (this was an issue)

splittime = finaltime.split(".")
time = splittime[0]
time = time.rstrip('\r\n')
scalartime = int(time)
timebuffer = 100
newfinaltime = scalartime + timebuffer
#print newfinaltime, "is x-axis"
rgbuffer = 0.5
finalrg = int(round(float(finalrg))) # this may result in some y-axis being relatively higher than others
newfinalrg = finalrg + rgbuffer # ...but it is necessary in case the value gets rounded down during conversion
#print newfinalrg, "is y-axis" 

xmin = 0
ymin = 0
xmax = 60000      # Max value for the x-axis in picoseconds
ymax = newfinalrg      # Max value for the y-axis in angstroms

show_legend = False
show_grid = False

folders = []
for rmsdfile2 in rmsdfiles:
	if "#" in rmsdfile2:
		x = 1 
	else:
		rmsdfp = rmsdfile2.split("rmsd")
		print rmsdfp[0]	
		folders.append(rmsdfp[0])

colors = ['#00ff00', '#0000ff', '#8a2be2', '#ff69b4', '#00ffff', '#000000', '#8b8378', '#7fffd4', '#a52a2a', '#8b7355', '#5f9ea0', '#7fff00', '#458b00', '#d2691e', '#cd5b45', '#6495ed', '#b8860b', '#006400', '#6e8b3d', '#8b4500', '#68228b', '#8fbc8f', '#483d8b', '#8b0a50', '#009acd', '#104e8b', '#b22222', '#228b22', '#eec900',  '#663399']

### LIST OF COLORS FOR REFERENCE
# green, blue, 'blueviolet', 'hotpink', cyan, black, antique white, 'aquamarine1',  'brown', 'burlywood4', 'cadetblue', 'chartreuse1',
#'chartreuse4', 'chocolate', 'coral3', 'cornflowerblue', 'darkgoldenrod', 'darkgreen', 'darkolivegreen4', 'darkorange4', 'darkorchid4',
#'darkseagreen', 'darkslateblue', 'deeppink4', 'deepskyblue3', 'dodgerblue4', 'firebrick', 'forestgreen', 'gold2', 'rebeccapurple'
###

pylab.figure(figsize=(16,8)) # 16x8 inches, 100 dpi is default so 1600x800 pixels
n = 0
for folder in folders:
    time = []
    rmsd = []
    rg = []
    print folder
    # get the rmsd data
    f = open("%s/rmsd.xvg" % folder)    
    for line in f:
        if line[0:1] != "@" and line[0:1] != "#":
            fields = line.split()
            time.append(fields[0])
            rmsd.append(fields[1])
    f.close()
    
    # get the rg data
    f = open("%s/gyrate.xvg" % folder)    
    for line in f:
        if line[0:1] != "@" and line[0:1] != "#":
            fields = line.split()
            rg.append(fields[1])
    f.close()

    pylab.plot(time, rmsd, colors[n], linestyle='--', label="rmsd %s" % folder)
    pylab.plot(time, rg, colors[n], label="rg %s" % folder)
    
    n += 1
    
pylab.axis([xmin, xmax, ymin, ymax])  
pylab.xlabel('time (ps)')
pylab.ylabel('distance (nm)')
if show_legend: pylab.legend()
pylab.title("RMSD/GYRATE %s" % Plot_title)
if show_grid: pylab.grid()
pylab.savefig("rmsd_gyrate%s.png" % FileName)
pylab.savefig("rg_gyrate%s.eps" % FileName)

