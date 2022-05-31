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
Z_NOx_Std[0] = np.array([34.39022718,
                54.90147382,
                75.02834747,
                87.09965742,
                92.46075295,
                93.73069472,
                91.99093926,
                89.50752759,
                84.37764107,
                64.88453892,
                25.15517762,
                ])

# 2hr
Z_NOx_Std[1] = np.array([32.44623391,
                52.15634614,
                72.18183785,
                84.41605923,
                90.42930498,
                92.74327914,
                91.37747929,
                88.82643728,
                84.29854052,
                71.56343239,
                35.56367262,
                ])

# 4hr
Z_NOx_Std[2] = np.array([33.92737164,
                54.04362142,
                72.18183785,
                83.51698572,
                89.09884962,
                91.5677844,
                90.62769488,
                88.03489989,
                83.46798473,
                71.03487165,
                41.24962815,
                ])

# 8hr
Z_NOx_Std[3] = np.array([34.04308552,
                52.57301731,
                70.49980944,
                81.47363684,
                87.12462553,
                90.32959661,
                89.43485605,
                86.93042911,
                82.2023759,
                69.97775018,
                40.34923931,
                ])


# 16hr
Z_NOx_Std[4] = np.array([33.23308833,
                51.22496354,
                68.94716783,
                79.83895774,
                85.69402837,
                89.34218103,
                88.75323386,
                86.19411526,
                81.60912176,
                69.55490159,
                40.34923931,
                ])

# Unaged
Z_NOx_Std_data[0] = np.array([21.75427093,
50.23230576,
75.95993243,
89.09260369,
95.00149351,
97.33554526,
96.93270014,
93.77814794,
85.16864658,
77.68205147,
66.2174107,
])

# 2hr
Z_NOx_Std_data[1] = np.array([17.9588555,
44.68077523,
72.13008313,
87.23588068,
93.44786499,
96.13340597,
94.35957637,
89.01051574,
84.29854052,
72.40912956,
56.33114123,
])

# 4hr
Z_NOx_Std_data[2] = np.array([20.04170543,
47.18080222,
72.56999825,
87.03154579,
92.81840224,
95.59581304,
93.69499474,
89.69160606,
83.05270683,
72.30341742,
59.25740496,
])

# 8hr
Z_NOx_Std_data[3] = np.array([19.92599154,
44.60724502,
70.49980944,
85.15166482,
91.41641702,
94.48301135,
92.50215591,
88.40305682,
81.60912176,
69.76632588,
55.6558496,
])


# 16hr
Z_NOx_Std_data[4] = np.array([17.9588555,
41.17583542,
68.94716783,
82.83586943,
90.12887957,
93.88742735,
92.16134481,
88.03489989,
81.01586763,
68.49778012,
53.4048775,
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

plt.savefig("NOx_Std.png")
plt.show()

print("\nDisplaying plot. Press enter to continue...(this closes the images)")
input()
plt.close()

# Copy X into Z to create the needed space
Z_NH3_Std = Z_NOx_Std.copy()
Z_NH3_Std_data = Z_NOx_Std.copy()


# Unaged
Z_NH3_Std[0] = np.array([34.82731933,
55.08490722,
75.6374936,
88.35466038,
94.28727777,
97.49227789,
97.80176844,
97.53334859,
97.29080611,
97.71661762,
98.26675148,
])

# 2hr
Z_NH3_Std[1] = np.array([33.0595175,
52.8916482,
73.08754546,
85.42547357,
91.93858498,
95.92495157,
96.50668627,
96.11594443,
95.59014425,
95.85608383,
96.71358073,
])

# 4hr
Z_NH3_Std[2] = np.array([34.15879941,
54.04362142,
72.95815866,
84.47054853,
90.35777512,
94.8905162,
95.22864467,
95.02988149,
94.26521001,
94.50296834,
95.72315301,
])

# 8hr
Z_NH3_Std[3] = np.array([34.15879941,
52.45046697,
70.88796984,
82.01852988,
88.55522269,
92.94703156,
94.03580583,
93.37317532,
92.4854476,
92.81157399,
94.59766696,
])


# 16hr
Z_NH3_Std[4] = np.array([33.23308833,
51.22496354,
69.20594144,
80.52007403,
86.8385061,
91.84990314,
93.01337255,
92.26870454,
91.29893933,
91.54302822,
93.2470837,
])


# Unaged
Z_NH3_Std_data[0] = np.array([21.23675257,
49.55958092,
75.80760131,
89.37301973,
95.0812105,
97.49227789,
97.86993065,
98.56418799,
98.81349173,
98.94287853,
99.32470837,
])

# 2hr
Z_NH3_Std_data[1] = np.array([15.8412914,
43.32046642,
72.95815866,
85.4309225,
91.94716857,
95.92495157,
96.52372683,
96.13435227,
97.01395418,
98.62574209,
98.87451395,
])

# 4hr
Z_NH3_Std_data[2] = np.array([17.1488583,
45.09744639,
72.56999825,
84.74299505,
90.38638706,
94.90618946,
95.39905021,
96.13435227,
97.82473483,
98.94287853,
99.54980558,
])

# 8hr
Z_NH3_Std_data[3] = np.array([17.9588555,
43.38174159,
70.88796984,
84.47054853,
90.84417815,
94.51435788,
95.73986131,
95.39803842,
97.2314807,
98.94287853,
99.32470837,
])


# 16hr
Z_NH3_Std_data[4] = np.array([15.87600557,
39.70523131,
65.97127142,
82.97209269,
90.27193929,
93.73069472,
93.86540029,
95.02988149,
97.42923207,
98.94287853,
99.32470837,
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

plt.savefig("NH3_Std.png")
plt.show()

print("\nDisplaying plot. Press enter to continue...(this closes the images)")
input()
plt.close()


# Copy X into Z to create the needed space
Z_N2O_Std = Z_NOx_Std.copy()
Z_N2O_Std_data = Z_NOx_Std.copy()


# Unaged
Z_N2O_Std[0] = np.array([0.69428331,
0.980402743,
1.164481206,
1.362232585,
1.570795685,
1.316554109,
1.141717169,
1.012431548,
0.984801867,
1.16283362,
1.373092982,
])

# 2hr
Z_N2O_Std[1] = np.array([0.717426087,
1.066187983,
1.319745367,
1.430344214,
1.459209107,
1.394920425,
1.226919943,
1.086062933,
1.067857446,
1.289688197,
2.048384612,
])

# 4hr
Z_N2O_Std[2] = np.array([0.771811613,
1.138492685,
1.358561407,
1.457588866,
1.48782105,
1.421564972,
1.276337552,
1.148649611,
1.117295291,
1.346772756,
2.111411831,
])

# 8hr
Z_N2O_Std[3] = np.array([0.774125891,
1.127463154,
1.371500087,
1.457588866,
1.48782105,
1.434103583,
1.28826594,
1.167057457,
1.137070429,
1.374257914,
2.151929329,
])


# 16hr
Z_N2O_Std[4] = np.array([0.728997476,
1.078443017,
1.319745367,
1.416721889,
1.430597163,
1.379247161,
1.243960498,
1.141286472,
1.107407722,
1.374257914,
2.160933218,
])


# Unaged
Z_N2O_Std_data[0] = np.array([0.347141655,
1.348053772,
2.328962412,
2.315795395,
2.131589773,
1.238187793,
0.766824964,
0.957208009,
1.364484514,
1.606824638,
1.508151308,
])

# 2hr
Z_N2O_Std_data[1] = np.array([0.356398766,
1.20099336,
2.044311451,
2.084215855,
1.888388256,
1.473286741,
1.056514395,
1.306957089,
2.017064065,
2.494806675,
2.453559591,
])

# 4hr
Z_N2O_Std_data[2] = np.array([0.354084488,
1.214473898,
1.85023125,
1.879880967,
1.759634511,
1.485825351,
1.087187394,
1.411881813,
2.195040306,
2.663946111,
2.566108196,
])

# 8hr
Z_N2O_Std_data[3] = np.array([0.370284432,
1.127463154,
1.669089729,
1.743657709,
1.673798681,
1.452911498,
1.644924751,
1.457901429,
2.195040306,
2.663946111,
2.543598475,
])


# 16hr
Z_N2O_Std_data[4] = np.array([0.335570267,
0.968147709,
1.487948208,
1.593812125,
1.559350908,
1.379247161,
1.141717169,
1.491035552,
2.254365719,
2.68508854,
2.498579033,
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

plt.savefig("N2O_Std.png")
plt.show()

print("\nDisplaying plot. Press enter to continue...(this closes the images)")
input()
plt.close()
