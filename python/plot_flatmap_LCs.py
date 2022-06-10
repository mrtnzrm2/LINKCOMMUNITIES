from re import A
from tkinter.tix import Tree
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import sys
import json
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import copy
import pandas as pd
from os import path, makedirs

np.set_printoptions(linewidth=240)

def unique(list1):
    x = np.array(list1)
    print(np.unique(x))

def format_areas(areas):
    for i, a in enumerate(areas):
        if ("." in a and a != "pro.st."):
            areas[i] = a.replace(".", "/")
    return areas

def getAreaIndex(area, labels):
    labels = np.asarray(labels)
    index = np.where(labels == area)[0]
    if (len(index) > 0):
        return index[0]
    return -1

if len(sys.argv) != 7:
    print("> [U] : {} [merged] [model] [imputation] [distance] [resolution] [num link coms] [link type]".format(sys.argv[0]))
    print(">       [merged] : true/false. If true the merged data is used")
    print(">       [model] : Used fln: normal, original, zz_model")
    print(">       [imputation] : true/false. Original data or imputed version")
    print(">       [distance] : imputation method (choose from phys or trac)")
    print(">       [num link coms] : Number of link communities")
    print(">       [link type] : Either: h, v, or wsbm")
 

    sys.exit(0)

merged = sys.argv[1]
model = sys.argv[2]
imputation = sys.argv[3]
distance = sys.argv[4]
k = sys.argv[5]
ltype = sys.argv[6]

print("> good work needs good preparation")

if(merged == "true"):
    foldername = "merged"
else:
    foldername = "removed"
if imputation == "true":
    filename = "fln"
    filepath = "../CSV/{}/imputation/{}/{}".format(foldername, distance, model)

print("> reading data")
data = pd.read_csv("{}/{}.csv".format(filepath, filename), dtype=float).to_numpy()
numberOfLinks = np.count_nonzero(data)

memberships = pd.read_csv("../CSV/{}/tables/{}/{}/{}/{}.csv".format(foldername,distance, model, k, ltype))
memberships["AREA"] = format_areas(list(memberships["AREA"]))
nodeNames = pd.read_csv("../CSV/{}/labels/{}/{}/imputation_labels.csv".format(foldername, distance, model))["x"]
nodeNames = list(nodeNames)
nodeNames = format_areas(nodeNames)

dorsal = pd.read_csv("../CSV/{}/labels/{}/{}/dorsal.csv".format(foldername, distance, model))["x"]
dorsal = list(dorsal)
ventral = pd.read_csv("../CSV/{}/labels/{}/{}/ventral.csv".format(foldername, distance, model))["x"]
ventral = list(ventral)

print("> reading flatmap data")
index, name, sequence, view, x, y, name2 = np.loadtxt("utils/flatmap/flatmapdataframe2c.csv", skiprows=1, delimiter=',', unpack=True, dtype=str)
name = [nme.lower() for nme in name]
name2 = [nm2.lower() for nm2 in name2]

print("> reading labels and their positions")
f = open('utils/flatmap/F99-107-centres.json')
data = json.load(f)
f.close()

print("> preparing data for creating polygons")
for index, n in enumerate(name):
    name[index] = n.replace("\"", "")
sequence = sequence.astype(int)
x = x.astype(float)
y = y.astype(float)
for index, n in enumerate(name2):
    name2[index] = n.replace("\"", "")

totalAreas = len(set(name2))
start = np.zeros(totalAreas)
end = np.zeros(totalAreas)
sindex = 0
eindex = 0
for index, s in enumerate(sequence):
    if(s == 1):
        start[sindex] = index
        sindex += 1
    if(s == 3):
        end[eindex] = index
        eindex += 1

start = start.astype(int)
end = end.astype(int)

print("> defining colormap and limits")
lower = 0
upper = 1
cm = copy.copy(matplotlib.cm.get_cmap("Set3"))

print("> creating layout")
fsize = 10 
fig = plt.figure(figsize=(fsize, fsize))

left, width = 0.04, 1
bottom, height = 0.1, 0.8
smallHeightWidth = 0.02
spacing = 0.04

#adding sfpd layout
flatmap = [left, bottom, width, height]
axFlatmap = fig.add_axes(flatmap)


print("> plotting")
for index in range(totalAreas):
    px = x[start[index]:(end[index] + 1)]
    py = y[start[index]:(end[index] + 1)]

    vertices = np.column_stack((px, py))

    aindex = getAreaIndex(name[start[index]], nodeNames)
    if(aindex != -1 and memberships["COMMSHIP"][aindex] != -1):
        color = cm(memberships["COMMSHIP"][aindex])
        color = list(color)
        # color[3] = memberships["P"][aindex]
        color = tuple(color)
    else:
        color = tuple((.5, .5, .5, 0.3))

    # if(nodeNames[aindex] in dorsal):
    #     pattern = '///'
    # elif(nodeNames[aindex] in ventral):
    #     pattern = '\\\\\\'
    # else:
    #     pattern = ''
    pattern = ''

    # if nodeNames[aindex] in dorsal:
    #     axFlatmap.add_patch(plt.Polygon(vertices, closed=True, edgecolor=[1, 0, 0, 0.5], linewidth=.5, facecolor=color, hatch=pattern))
    # elif nodeNames[aindex] in ventral:
    #     axFlatmap.add_patch(plt.Polygon(vertices, closed=True, edgecolor=[0, 0, 1, 0.5], linewidth=.5, facecolor=color, hatch=pattern))
    # else:
    axFlatmap.add_patch(
        plt.Polygon(
            vertices,
            closed=True,
            edgecolor=[0, 0, 0],
            linewidth=.5,
            facecolor=color,
            hatch=pattern
        )
    )

for d in data:
    axFlatmap.text(data[d][0], data[d][1], d, fontsize=6)

axFlatmap.set_aspect("equal")
axFlatmap.autoscale()
axFlatmap.set_axis_off()

plot_main_folder = "../plots/{}/{}/{}/AVERAGE_full_l10/".format(foldername, distance, model)
subfolder = "flatmap"
if ~path.exists("{}/{}p".format(plot_main_folder, subfolder)):
    makedirs("{}/{}".format(plot_main_folder, subfolder), exist_ok=True)
plt.savefig(
    "{}/flatmap/k_{}_{}.pdf".format(
        plot_main_folder, k, ltype
    )
)
plt.close()

print("> done!")
