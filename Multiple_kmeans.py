#! /usr/bin/python

import commands
import sys

# data input file
# change this for each project p*_vectors_cut.txt
#cluster_data = "p1761_5vectors_nonnorm.txt"
cluster_data = sys.argv[1]

# fields to take from input file
# keep this the same
fields = 17

# starting cluster centers
# keep this the same
clusters = 100

# max iterations w/o convergence
iterations = 500

# multiplier for -mult flag
multiplier = 1

# change range to set total number of clustering trials
for n in range(1,102):
    name = "trial.%s" % (n)
    command = "./Kmeans_clustering_BChE.pl -iter %s -data %s -k %s -name %s -nume %s" % (iterations, cluster_data, clusters, name, fields)
    status, output = commands.getstatusoutput(command)
    if status != 0:
        print output
        sys.exit()
