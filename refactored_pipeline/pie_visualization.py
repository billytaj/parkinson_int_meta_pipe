import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# Labels and imput data
Taxon = ['Frogs', 'Hogs', 'Dogs', 'Logs', 'Bogs', "Pogs", "Slogs", "Cogs"]
Abundance = [5, 15, 30, 10, 20, 35, 10, 25]
RPKM = [10, 15, 45, 60, 10, 15, 45, 60]

#Generation of colour map
a=np.random.random(len(Taxon))
cs=cm.nipy_spectral(np.arange(len(Taxon)) / len(Taxon))

#legend variable
legend = []

# Generation of figure
for x in range(len(Taxon)):
    radius = 1/max(RPKM)*RPKM[x]
    patches, texts, autotexts = plt.pie(Abundance, labels=Taxon, autopct='%1.1f%%\nRPKM: ' + str(RPKM[x]), pctdistance=1.3, radius=radius, colors=cs)

    for n in range(len(patches)):
        if n != x:
            patches[n].set_color((0, 0, 0, 0))
        else:
            patches[n].set_label(Taxon[n] + "\n" + "Abundance: " + str(autotexts[n].get_text()) )#+ "\n" + "RPKM: " + str(RPKM[n]))
            legend.append(patches[n])
    for n in range(len(texts)):
        if True:#n != x:
            texts[n].set_color((0, 0, 0, 0))
    for n in range(len(autotexts)):
        if True:#n != x:
            autotexts[n].set_color((0, 0, 0, 0))
        else:
            if Abundance[n] < 0.1*max(Abundance) or RPKM[n] < 0.1*max(RPKM):
                autotexts[n].set_size("xx-small")
            elif Abundance[n] < 0.3*max(Abundance) or RPKM[n] < 0.3*max(RPKM):
                autotexts[n].set_size("x-small")
            elif Abundance[n] < 0.5*max(Abundance) or RPKM[n] < 0.5*max(RPKM):
                autotexts[n].set_size("small")

for x in range(5):
    fifth = max(RPKM)/5
    patches, texts = plt.pie([1], labels=[str('%1.0f' % (fifth*(5-x)))], labeldistance=1.1, startangle=270, textprops=dict(fontsize="x-small"), wedgeprops=dict(edgecolor='black', linewidth=0.5), radius=0.2*(5-x))
    patches[0].set_facecolor((0, 0, 0, 0))

patches, texts = plt.pie([1,0], startangle=90, wedgeprops=dict(edgecolor='black', linewidth=0.5))
patches[0].set_facecolor((0, 0, 0, 0))

#Setting final properties and building image
plt.legend(loc=6, bbox_to_anchor=(-0.35, 0.5), handles=legend)
plt.axis('scaled')
plt.show()
#plt.savefig("test.png", format="png")