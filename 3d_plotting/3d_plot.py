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

# Plot the surface.
surf = ax.plot_surface(X, Y, Z_NOx_Std, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

# Adjust angles of view
ax.azim = -135
ax.elev = 30

# Labeling Axis
ax.set_xlabel('Temperature (oC)')
ax.set_ylabel('Aging Time (hr)')
ax.set_zlabel('NOx Conversion (%)')

plt.show()

print("\nDisplaying plot. Press enter to continue...(this closes the images)")
input()
plt.close()

# Copy X into Z to create the needed space
Z_NH3_Std = Z_NOx_Std.copy()


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


# Plot the surface.
fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
surf2 = ax2.plot_surface(X, Y, Z_NH3_Std, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig2.colorbar(surf2, shrink=0.5, aspect=5)

# Adjust angles of view
ax2.azim = -135
ax2.elev = 30

# Labeling Axis
ax2.set_xlabel('Temperature (oC)')
ax2.set_ylabel('Aging Time (hr)')
ax2.set_zlabel('NH3 Conversion (%)')

plt.show()

print("\nDisplaying plot. Press enter to continue...(this closes the images)")
input()
plt.close()
