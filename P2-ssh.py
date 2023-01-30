import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as cdf

# load data file
df = cdf.Dataset("dt_global_allsat_phy_l4_20170515_20170914.nc", "r")

# load data from file
h = df["adt"][0,:,:]
λ = df["longitude"][:]
ϕ = df["latitude"][:]

# center of domain (Northeast Pacific)
λ0 = 145.0
ϕ0 = 32.5

# Cartesian coordinates
a = 6371e3
x = a*np.cos(np.deg2rad(ϕ0))*np.deg2rad(λ-λ0)
y = a*np.deg2rad(ϕ-ϕ0)

# cut domain
ix0 = np.argmax(x > -1e6)
ix1 = np.argmin(x < 1e6)
iy0 = np.argmax(y > -1e6)
iy1 = np.argmin(y < 1e6)
x = x[ix0:ix1]
y = y[iy0:iy1]
h = h[iy0:iy1,ix0:ix1]

# plot SSH
fig, ax = plt.subplots(1, 1)
ax.set_aspect(1)
img = ax.pcolormesh(1e-3*x, 1e-3*y, h, vmin=0, vmax=2, rasterized=True, cmap="inferno")
plt.colorbar(img, ax=ax, location="right", shrink=0.7, label="sea surface elevation (m)")
ax.set_xlabel("zonal coordinate (km)")
ax.set_ylabel("meridional coordinate (km)")
fig.tight_layout()
fig.savefig("ssh.pdf", dpi=300)
