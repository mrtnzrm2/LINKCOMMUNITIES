from re import A
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import sys
import json
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import copy
import pandas as pd

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
    print(">       [merged] : true/false, if true the merged data is used")
    print(">       [model] : Which Latent variable was used?")
    print(">       [imputation] : Original data or imputed version")
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

memberships = pd.read_csv("../CSV/{}/tables/{}/{}/{}/{}.csv".format(foldername, distance, model, k, ltype))
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
cm = copy.copy(matplotlib.cm.get_cmap("Dark2"))

print("> creating layout")
fsize = 7.5
fig = plt.figure(
  figsize=(fsize, fsize),
  facecolor="#E7DDD7"
)

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
        color[3] = memberships["P"][aindex]
        color = tuple(color)
        ed_color = memberships["LC"][aindex]
        if ed_color == "#fc8d62":
          pattern = "/"
        elif ed_color == "#8da0cb":
          pattern = "\\"
        elif ed_color == "#66c2a5":
          pattern = "|"
        else:
          pattern = "+"
    else:
        color = tuple((1, 1, 1, 0.3))
    pattern = ''
    axFlatmap.add_patch(
      plt.Polygon(
        vertices,
        closed = True,
        edgecolor = [0, 0, 0, 0.5],
        linewidth = 0.5,
        facecolor = ed_color,
        hatch = pattern
      )
    )

for d in data:
    axFlatmap.text(data[d][0], data[d][1], d, fontsize=6)

axFlatmap.set_aspect("equal")
axFlatmap.autoscale()
axFlatmap.set_axis_off()


plt.savefig("../plots/{}/{}/{}/AVERAGE_full_l10/flatmap/{}_clc_k_{}_{}_b.pdf".format(foldername, distance, model, filename, k, ltype))
plt.close()

print("> done!")
