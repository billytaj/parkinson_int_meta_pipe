import sys
import matplotlib
# Switching to backend that doesn't display to the user
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# Turn interactive plotting off
plt.ioff()
from matplotlib import cm
import numpy as np
import pandas as pd
import os.path
import multiprocessing as mp

# Assigning script arguments
RPKM_table = sys.argv[1]
Output_folder = sys.argv[2]

# Labels and imput data
RPKM_df = pd.read_table(RPKM_table, index_col=0)
del RPKM_df["RPKM"]
del RPKM_df["Piechart"]
RPKM_df.sort_index(inplace=True)
valid_cols = []
for n in range(len(RPKM_df.columns)):
    if RPKM_df.sum()[n] != 0.0:
        valid_cols.append(n)
else:
    subset_RPKM_df = RPKM_df.iloc[:, valid_cols]
Taxon = subset_RPKM_df.columns
Abundance = subset_RPKM_df.sum()
RPKMs = [subset_RPKM_df.sum()]
for n in range(len(subset_RPKM_df)):
    RPKMs.append(subset_RPKM_df.iloc[n])

#Generation of colour map
cs=cm.nipy_spectral(np.arange(len(Taxon)) / len(Taxon))

#legend variable
legend = []

#function to produce pie chart
def pie(RPKM, Output_file):
    # Generation of figure
    for x in range(len(Taxon)):
        if max(RPKM)*RPKM[x] > 0:
            radius = 1/max(RPKM)*RPKM[x]
        else:
            radius = 0
        patches, texts, autotexts = plt.pie(Abundance, labels=Taxon, autopct='%1.1f%%\nRPKM: ' + str('%1.0f' % RPKM[x]), pctdistance=1.3, radius=radius, colors=cs)
        for n in range(len(patches)):
            if n != x:
                patches[n].set_color((0, 0, 0, 0))
            else:
                patches[n].set_label(Taxon[n] + "\n" + "Abundance: " + str(autotexts[n].get_text()) )
                legend.append(patches[n])
        for n in range(len(texts)):
            texts[n].set_color((0, 0, 0, 0))
        for n in range(len(autotexts)):
            autotexts[n].set_color((0, 0, 0, 0))

    for x in range(5):
        fifth = max(RPKM)/5
        patches, texts = plt.pie([1], labels=[str('%1.0f' % (fifth*(5-x)))], labeldistance=1.1, startangle=270, textprops=dict(fontsize="x-small"), wedgeprops=dict(edgecolor='black', linewidth=0.5), radius=0.2*(5-x))
        patches[0].set_facecolor((0, 0, 0, 0))

    patches, texts = plt.pie([1,0], startangle=90, wedgeprops=dict(edgecolor='black', linewidth=0.5))
    patches[0].set_facecolor((0, 0, 0, 0))

    #Setting final properties and building image
    plt.legend(loc="right", fontsize="xx-small", bbox_to_anchor=(0, 0.5), handles=legend)
    plt.axis('scaled')
    #plt.show()
    plt.savefig(os.path.join(Output_folder, Output_file + ".png"), format="png")

#Loop to build chart for each EC
ChartProcesses = []
for RPKM in RPKMs:
    if RPKM.name == None:
        Output_file = "All_EC"
    else:
        Output_file = RPKM.name
    process = mp.Process(target=pie, args=(RPKM, Output_file))
    process.start()
    ChartProcesses.append(process)
for item in ChartProcesses:
    item.join()
