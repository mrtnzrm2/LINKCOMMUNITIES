import numpy as np
import pandas as pd
import matplotlib
import copy


if __name__ == "__main__":
    memb = pd.read_csv("../CSV/merged/tables/tracto2016/zz_model/4/h_1.csv")
    print(memb["COMMSHIP"] == -1)
    # print(memb["COMMSHIP"].to_numpy())
    # dorsal = pd.read_csv("../CSV/merged/labels/tracto2016/zz_model/NONULL/ventral.csv")["x"]
    # print(list(dorsal))

    # fln = pd.read_csv('../CSV/merged/imputation/tracto2016/zz_model/NONULL/fln.csv', dtype=float).to_numpy()
    # print(fln)

    # cm = copy.copy(matplotlib.cm.get_cmap("Set1"))
    # print(type(cm(3)[3]))