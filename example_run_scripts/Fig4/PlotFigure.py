import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.offsetbox as offsetbox
import imageio
import os

SYSTEMS = ["5", "mixed", "15", "20", "40", "upper"]
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "black"]

# load images from a gif
gif_filename = "results/upper/path/figs/path34.gif"
# load it as grayscale
gif = imageio.mimread(gif_filename)
# save the images as png
for (i, frame) in enumerate(gif):
    # save image
    imageio.imsave("results/upper/path/figs/path_"+str(i)+".png", frame)


fig = plt.figure(figsize=(8, 4))
ax = plt.gca()

# put the 1st, 5th and 10th images inside the plot
for (i,x,y) in zip([0, 5, 9], [0.07, 0.555, 0.92], [20.0, -3.0, 20.0]):
    ab = offsetbox.AnnotationBbox(offsetbox.OffsetImage(mpimg.imread("results/upper/path/figs/path_"+str(i)+".png"), zoom=0.055), (x, y), frameon=False, zorder=-1)
    ax.add_artist(ab)


for (SYSTEM, COLOR) in zip(SYSTEMS, COLORS):
    # find the last energy path in the directory
    files = os.listdir("results/" + SYSTEM + "/path/figs")
    counts = []
    for file in files:
        if "energy_path" in file:
            if "interp" not in file:
                counts.append(int(file.split("energy_path")[1][:-4]))
    counter_image = np.max(counts)
    # load the data
    d = np.loadtxt("results/" + SYSTEM + "/path/figs/energy_path"+str(counter_image)+".csv", delimiter=",", skiprows=1)
    xi_arr, energy_arr = d[:, 0], d[:, 1]
    d = np.loadtxt("results/" + SYSTEM + "/path/figs/energy_path_interpolated"+str(counter_image)+".csv", delimiter=",", skiprows=1)
    xi_interp_arr, en_interp_arr = d[:, 0], d[:, 1]

    # plot the energy with the cubic interpolation
    ax.plot(xi_interp_arr, en_interp_arr, "-", color=COLOR, linewidth=2.0, alpha=1.0, zorder=0)
    
    # plot the energy as scatter
    plt.scatter(xi_arr, energy_arr, s=60, color=COLOR, linewidths=0.75, alpha=1.0, zorder=1)

ax.set_xlim(0.0, 1.0)
[i.set_linewidth(2.5) for i in ax.spines.values()]
# make axis load last
[spine.set_zorder(10) for spine in ax.spines.values()]
ax.xaxis.set_tick_params(width=1.5)
ax.yaxis.set_tick_params(width=1.5)
ax.set_xlabel(r"$\xi$", fontsize=20)
ax.set_ylabel(r"$\Delta$F [k$_\mathrm{B}$T]", fontsize=20)
ax.tick_params(labelsize=20)
#plt.legend(fontsize=13, loc="lower left", ncol=1, frameon=false, borderaxespad=0.0, handletextpad=0.0, columnspacing=0.0)
# add table as a legend
text = [["", r"$\phi$", r"$\tilde{J}_s$ [nm$^{-1}$]"], ["-", "0.05", "0.8"], ["-", "0.1", "0.4"], ["-", "0.15", "0.267"], ["-", "0.2", "0.2"], ["-", "0.4", "0.1"]]
col_widths = [0.05, 0.15, 0.2]
table = plt.table(cellText=text, loc="lower left", cellLoc="center", rowLoc="center", colLoc="center", edges="horizontal", bbox=[0.03, 0.1, 0.45, 0.45], fontsize=13, colWidths=col_widths)
table.auto_set_font_size(False)
table.set_fontsize(12)
# change first column text colors
for i in range(1, 6):
    table[(i, 0)].get_text().set_color(COLORS[i-1])
    table[(i, 0)].get_text().set_fontsize(30)
# remove most of the table borders, only keep the first row
for i in range(1, 6):
    for j in range(0, 3):
        table[(i, j)].set_linewidth(0.0)
table[(0, 0)].set_linewidth(0.0)
# add a one row, one column table below the above one
text = [["-", r"no mix, $J^u_s=0.04$ nm$^{-1}$"]]
colors_cells = [["black", "white"]]
col_widths = [0.05, 0.35]
table1 = plt.table(cellText=text, loc="lower left", cellLoc="center", rowLoc="center", colLoc="center", edges="horizontal", bbox=[0.03, 0.025, 0.45, 0.075], fontsize=15, colWidths=col_widths)
table1.auto_set_font_size(False)
table1.set_fontsize(12)
# change first column text colors
table1[(0, 0)].set_linewidth(0.0)
table1[(0, 0)].get_text().set_color("black")
table1[(0, 0)].get_text().set_fontsize(30)
plt.savefig("results/fig4.svg", bbox_inches="tight", dpi=1000)
plt.savefig("results/fig4.png", bbox_inches="tight", dpi=300)
plt.close("all")