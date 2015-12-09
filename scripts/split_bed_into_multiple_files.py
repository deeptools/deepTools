#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script is usted to split a bed file into
smaller beds following the convention that I introduced
in which a # indicates a group in the bed file

Each bed group is saved under a different name
corresponding to the name of the group

Useful to split the results of heatmapper when the clustering
option is used.

:Authors: fidel.ramirez@gmail.com
"""

import sys

i = 0
tempArray = []

for line in sys.stdin:
    if line[0] == '#':
        clusterName = line[1:].strip()
        tempArray.append("#" + clusterName + "\n")
        open(clusterName + ".bed", 'w').write("".join(tempArray))
        tempArray = []
        continue

    tempArray.append(line)

if len(tempArray) > 0:
    clusterName = "no_name"
    open(clusterName + ".bed", 'w').write("".join(tempArray))
