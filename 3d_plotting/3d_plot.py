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
Z = X.copy()

# Fill in Z values
for rx, x_row in enumerate(X):
    for cx, val_x in enumerate(x_row):
        for ry, y_row in enumerate(Y):
            for cy, val_y in enumerate(y_row):
                #print(cy)
                pass

# Unaged
Z[0] = np.array([34.39022718,
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
Z[1] = np.array([32.44623391,
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
Z[2] = np.array([33.92737164,
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
Z[3] = np.array([34.04308552,
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
Z[4] = np.array([33.23308833,
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
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
