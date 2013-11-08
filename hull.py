#!/usr/bin/python
import sys
import subprocess

hull_path = "./hull"

def get_alpha_shape(infile,radius):
    f = open(infile)
    line = f.readline()
    points = []
    while line:
        fields = line.split()
        points.append([ float( fields[0]), float( fields[1]) ])
        line = f.readline()

    command = "%s -A -aa %s -r -m10000000 -oN -oFpoints < %s" % (hull_path, str(radius), infile)
    print sys.stderr, "Running command: %s" % command
    retcode = subprocess.call(command, shell=True)
    results_file = open("points-alf")
    results_file.next()
    results_indices = [[int(i) for i in line.rstrip().split()] for line in results_file]
    results_file.close()
    return [(points[i], points[j]) for i,j in results_indices ]

def google( file ):
    infile=file

    radius=1
    alpha = get_alpha_shape(infile,radius)

    print "\nFor google:"
    idx = 0;
    for a in alpha:
        print 'new google.maps.LatLng(',a[0][1], ',', a[0][0],'),'
        print 'new google.maps.LatLng(',a[1][1], ',', a[1][0],'),'
        idx+=1

if __name__ == '__main__':
    google('google.txt')
