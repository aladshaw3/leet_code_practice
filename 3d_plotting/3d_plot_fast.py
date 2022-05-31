import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.array([150.,175.,200.,225.,250.,300.,350.,400.,450.,500.,550.])
Y = np.array([0.,2.,4.,8.,16.])

X, Y = np.meshgrid(X, Y)

# Copy X into Z to create the needed space
Z_NOx_Std = X.copy()
Z_NOx_Std_data = X.copy()

# import the data

# Unaged
Z_NOx_Std[0] = np.array([33.01323194,
58.896615,
80.12748128,
89.48628891,
93.47790753,
94.56137767,
92.04206093,
91.42194361,
87.44278744,
72.28861772,
37.95195402,
])

# 2hr
Z_NOx_Std[1] = np.array([34.563798,
59.6809372,
78.65117789,
87.65817278,
91.8312902,
93.66800167,
92.74072367,
90.66722191,
87.00773441,
76.19362444,
47.46231115,
])

# 4hr
Z_NOx_Std[2] = np.array([30.4443837,
55.52648057,
76.8915174,
86.20058391,
90.47222289,
92.86866525,
92.04206093,
89.93090806,
86.35515486,
75.89763042,
47.1921945,
])

# 8hr
Z_NOx_Std[3] = np.array([30.68738285,
56.00442691,
75.28712107,
84.60677179,
88.98440184,
91.69317051,
91.13891152,
88.95529221,
84.9708952,
74.62908466,
46.42686398,
])


# 16hr
Z_NOx_Std[4] = np.array([29.29881623,
53.18576902,
73.99325306,
82.97209269,
87.55380468,
90.90950735,
90.45728933,
88.40305682,
84.37764107,
74.20623607,
46.65196119,
])

# Unaged
Z_NOx_Std_data[0] = np.array([22.79569589,
68.75456458,
88.46257898,
94.98017292,
96.56656681,
97.22583241,
96.93270014,
93.59406948,
87.24503606,
72.71146631,
62.58884367,
])

# 2hr
Z_NOx_Std_data[1] = np.array([25.59597191,
70.75948819,
88.45869738,
94.64642594,
96.80976833,
97.80574315,
97.56320067,
90.75926115,
84.37764107,
72.30341742,
57.3665884,
])

# 4hr
Z_NOx_Std_data[2] = np.array([19.68299239,
65.73492413,
87.66943789,
94.00617663,
96.39489515,
96.31678315,
96.67709182,
91.77169269,
82.00462452,
68.70920441,
53.85507192,
])

# 8hr
Z_NOx_Std_data[3] = np.array([23.62883587,
68.87221291,
87.07425861,
93.2569487,
95.70820851,
97.02207999,
95.91026686,
90.61199837,
80.02711073,
65.96068858,
49.57822492,
])


# 16hr
Z_NOx_Std_data[4] = np.array([19.23170823,
64.46040057,
85.50867832,
92.23527426,
95.13596964,
97.02207999,
95.39905021,
90.05976298,
78.84060246,
63.42359705,
46.8770584,
])


# Plot the surface.
surf = ax.plot_surface(X, Y, Z_NOx_Std, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

surf = ax.scatter(X, Y, Z_NOx_Std_data, color='black')

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

# Adjust angles of view
ax.azim = -135
ax.elev = 30

ax.set_zlim(0, 100)

# Labeling Axis
ax.set_xlabel('Temperature ($^o$C)')
ax.set_ylabel('Aging Time (hr)')
ax.set_zlabel('NOx Conversion (%)')

plt.xticks(rotation = -10)
plt.yticks(rotation = 15)

plt.savefig("NOx_Fast.png")
plt.show()

print("\nDisplaying plot. Press enter to continue...(this closes the images)")
input()
plt.close()

# Copy X into Z to create the needed space
Z_NH3_Std = Z_NOx_Std.copy()
Z_NH3_Std_data = Z_NOx_Std.copy()


# Unaged
Z_NH3_Std[0] = np.array([35.08451049,
62.98979645,
83.89134331,
92.29657473,
95.82694807,
97.83708968,
98.02329565,
97.77265059,
97.52810776,
97.86461463,
98.31177092,
])

# 2hr
Z_NH3_Std[1] = np.array([37.34093124,
64.19078981,
82.95975835,
90.68232912,
94.1488576,
96.58322862,
96.89861903,
96.52091704,
96.04497242,
96.1943627,
96.87114878,
])

# 4hr
Z_NH3_Std[2] = np.array([33.23308833,
60.5387896,
81.6270743,
89.49718677,
92.84701418,
95.83091199,
95.99546963,
95.52689334,
94.85846415,
95.13724123,
95.94825022,
])

# 8hr
Z_NH3_Std[3] = np.array([33.11737444,
61.15154131,
79.94504589,
87.73990673,
91.55947674,
94.20089262,
94.71742802,
94.10948918,
93.27645312,
93.44584687,
94.82276417,
])


# 16hr
Z_NH3_Std[4] = np.array([31.8445217,
58.21033308,
78.65117789,
86.37767415,
90.12887957,
93.26049683,
93.86540029,
93.18909686,
92.28769622,
92.3887254,
93.69727812,
])


# Unaged
Z_NH3_Std_data[0] = np.array([23.86026364,
66.85013225,
86.80254633,
92.30338589,
95.75971001,
97.83708968,
98.10849842,
98.7666743,
99.40674586,
99.36572712,
99.54980558,
])

# 2hr
Z_NH3_Std_data[1] = np.array([23.82554947,
65.41629324,
83.03739043,
90.7231961,
94.16316357,
96.58322862,
97.06902458,
97.9751369,
98.81349173,
98.73145423,
98.87451395,
])

# 4hr
Z_NH3_Std_data[2] = np.array([18.53742492,
62.86724611,
85.12051792,
89.64703235,
92.8899321,
95.87793178,
97.35871401,
98.80348999,
99.30787017,
99.47143926,
99.66235418,
])

# 8hr
Z_NH3_Std_data[3] = np.array([18.30599715,
66.91140742,
85.24990472,
91.96282775,
94.7067905,
96.0816842,
97.44391679,
98.52737229,
99.20899448,
99.57715141,
99.77490279,
])


# 16hr
Z_NH3_Std_data[4] = np.array([17.26457219,
62.00939371,
83.69726311,
91.28171146,
93.56231277,
95.92495157,
97.27351124,
98.34329383,
99.20899448,
99.57715141,
99.77490279,
])


# Plot the surface.
fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
surf2 = ax2.plot_surface(X, Y, Z_NH3_Std, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

surf2 = ax2.scatter(X, Y, Z_NH3_Std_data, color='black')

# Add a color bar which maps values to colors.
#fig2.colorbar(surf2, shrink=0.5, aspect=5)

# Adjust angles of view
ax2.azim = -135
ax2.elev = 30

# Labeling Axis
ax2.set_xlabel('Temperature ($^o$C)')
ax2.set_ylabel('Aging Time (hr)')
ax2.set_zlabel('NH$_3$ Conversion (%)')

ax2.set_zlim(0, 100)

plt.xticks(rotation = -10)
plt.yticks(rotation = 15)

plt.savefig("NH3_Fast.png")
plt.show()

print("\nDisplaying plot. Press enter to continue...(this closes the images)")
input()
plt.close()


# Copy X into Z to create the needed space
Z_N2O_Std = Z_NOx_Std.copy()
Z_N2O_Std_data = Z_NOx_Std.copy()


# Unaged
Z_N2O_Std[0] = np.array([1.041424965,
2.083355829,
2.846509615,
3.405581463,
3.67663471,
3.839949483,
3.748922047,
3.644753572,
3.599075096,
3.551928148,
3.128851221,
])

# 2hr
Z_N2O_Std[1] = np.array([0.991667995,
1.924040383,
2.730061495,
3.228491227,
3.562186937,
3.839949483,
3.817084266,
3.75520065,
3.717725923,
3.784494872,
4.029240062,
])

# 4hr
Z_N2O_Std[2] = np.array([0.853968472,
1.740214869,
2.613613374,
3.11951262,
3.462045135,
3.839949483,
3.851165376,
3.810424189,
3.796826475,
3.890207019,
4.164298388,
])

# 8hr
Z_N2O_Std[3] = np.array([0.808840057,
1.691194732,
2.471287893,
2.942422384,
3.376209305,
3.792929694,
3.868205931,
3.847239882,
3.875927026,
3.995919167,
4.299356714,
])


# 16hr
Z_N2O_Std[4] = np.array([0.74519742,
1.556389354,
2.419533173,
2.915177732,
3.304679447,
3.745909904,
3.834124821,
3.865647728,
3.895702164,
4.059346455,
4.389395598,
])


# Unaged
Z_N2O_Std_data[0] = np.array([0.578569425,
2.818657886,
3.881604021,
3.84149589,
3.705246653,
3.260038745,
2.965056528,
3.037294644,
3.104696648,
2.811943117,
2.318501265,
])

# 2hr
Z_N2O_Std_data[1] = np.array([0.630640674,
3.100523675,
3.71340118,
3.51456007,
3.347597362,
3.495137693,
3.731881493,
4.178581116,
4.31098006,
3.974776737,
3.376458152,
])

# 4hr
Z_N2O_Std_data[2] = np.array([0.356398766,
2.708362577,
3.454627578,
3.255735878,
3.133007788,
3.401098114,
3.800043712,
4.344251733,
4.449406025,
4.080488884,
3.443987315,
])

# 8hr
Z_N2O_Std_data[3] = np.array([0.458226985,
2.732872646,
3.208792657,
2.942422384,
2.947030156,
3.275712008,
3.783003157,
4.362659579,
4.449406025,
4.038204025,
3.466497037,
])


# 16hr
Z_N2O_Std_data[4] = np.array([0.446655596,
2.438751823,
2.950019056,
2.833443777,
2.80397044,
3.181672429,
3.731881493,
4.344251733,
4.429630887,
4.080488884,
3.376458152,
])


# Plot the surface.
fig3, ax3 = plt.subplots(subplot_kw={"projection": "3d"})
surf3 = ax3.plot_surface(X, Y, Z_N2O_Std, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

surf3 = ax3.scatter(X, Y, Z_N2O_Std_data, color='black')

# Add a color bar which maps values to colors.
#fig3.colorbar(surf3, shrink=0.5, aspect=5)

# Adjust angles of view
ax3.azim = -135
ax3.elev = 30

ax3.set_zlim(0, 6)

#ax3.pcolormesh(X, Y, Z_N2O_Std, vmin=0., vmax=5., cmap=cm.coolwarm)

# Labeling Axis
ax3.set_xlabel('Temperature ($^o$C)')
ax3.set_ylabel('Aging Time (hr)')
ax3.set_zlabel('NH$_3$-to-N$_2$O Conversion (%)')

plt.xticks(rotation = -10)
plt.yticks(rotation = 15)

plt.savefig("N2O_Fast.png")
plt.show()

print("\nDisplaying plot. Press enter to continue...(this closes the images)")
input()
plt.close()
