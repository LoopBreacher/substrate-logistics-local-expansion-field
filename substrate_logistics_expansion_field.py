import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage import gaussian_filter, map_coordinates
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u

# --- Hardware Constants (Substrate Logistics) ---
H_GLOBAL = 67.42          # km/s/Mpc baseline
K_OMEGA = 1.11587e-37     # m^3 / (node * s^2)
M_P = 1.672622e-27        # kg / node
MPC_TO_METER = 3.08567758149137e22
SEC_TO_KM_S_MPC = 3.24078e-20
SOLAR_MASS_KG = 1.98847e30
M_B_SUN = 5.44            # Solar absolute B-band magnitude
MASS_TO_LIGHT_B = 0.44   # Stellar M/L_B ratio

# --- Step 1: Query & Clean VizieR Catalog ---
print("Querying Karachentsev et al. (2004) Catalog of Neighboring Galaxies (CNG)...")
viz = Vizier(row_limit=-1)
cat_list = viz.get_catalogs('J/AJ/127/2031')

if not cat_list:
    raise RuntimeError("Could not retrieve catalog from VizieR.")

table = cat_list[0]

def clean_numeric_column(col):
    out = []
    for val in col:
        if np.ma.is_masked(val) or val is None:
            out.append(np.nan)
        else:
            try:
                out.append(float(val))
            except (ValueError, TypeError):
                out.append(np.nan)
    return np.array(out, dtype=np.float64)

names = [str(val) for val in table['Name']]
dist = clean_numeric_column(table['Dist'])
bmag = clean_numeric_column(table['Bmag'])

valid_mask = np.isfinite(dist) & np.isfinite(bmag) & (dist > 0)

names_clean = np.array(names)[valid_mask]
dist_clean = dist[valid_mask]
bmag_clean = bmag[valid_mask]
ra_clean = [str(val) for val in table['RAJ2000'][valid_mask]]
dec_clean = [str(val) for val in table['DEJ2000'][valid_mask]]

# Parse Coordinates to 3D Cartesian (Mpc)
sky_coords = SkyCoord(ra=ra_clean, dec=dec_clean, unit=(u.hourangle, u.deg), distance=dist_clean*u.Mpc)
cart = sky_coords.cartesian
X, Y, Z = cart.x.value, cart.y.value, cart.z.value

# --- Step 2: Distance Modulus & True Stellar Mass Calculation ---
# Convert Apparent Magnitude (bmag) to Absolute Magnitude (M_B)
M_B = bmag_clean - 5.0 * np.log10(dist_clean * 1e6) + 5.0

# Calculate B-band Luminosity (L_B) and Compiled Stellar Mass (M_compiled)
L_B = 10.0 ** (0.4 * (M_B_SUN - M_B))
M_compiled_kg = L_B * MASS_TO_LIGHT_B * SOLAR_MASS_KG
nodes = M_compiled_kg / M_P

# --- Step 3: Spatial Grid & Field Evaluation ---
GRID_LIM = 10.0  # 10 Mpc radius
RES = 120        # 120^3 grid resolution
grid_edges = np.linspace(-GRID_LIM, GRID_LIM, RES + 1)
dx_mpc = (2 * GRID_LIM) / RES
cell_vol_m3 = (dx_mpc * MPC_TO_METER) ** 3

galaxy_points = np.vstack([X, Y, Z]).T
grid_3d, _ = np.histogramdd(galaxy_points, bins=[grid_edges, grid_edges, grid_edges], weights=nodes)

# Gaussian smoothing kernel (1.0 Mpc scale)
sigma_mpc = 1.0
sigma_cells = sigma_mpc / dx_mpc
smoothed_nodes = gaussian_filter(grid_3d, sigma=sigma_cells)

# Volumetric Node Density rho_N (nodes / m^3)
rho_N = smoothed_nodes / cell_vol_m3

# Calculate H_local Field Grid
delta_H_rad_s = np.sqrt((8.0 * np.pi * K_OMEGA * rho_N) / 3.0)
delta_H_km_s_mpc = delta_H_rad_s / SEC_TO_KM_S_MPC
H_local_grid = H_GLOBAL + delta_H_km_s_mpc

# --- Step 4: Interpolate H_local at Exact Galaxy Positions & Export CSV ---
grid_coords = (galaxy_points + GRID_LIM) / dx_mpc - 0.5
galaxy_hlocal = map_coordinates(H_local_grid, grid_coords.T, order=1, mode='nearest')

df_out = pd.DataFrame({
    'Galaxy_Name': names_clean,
    'Cartesian_X_Mpc': np.round(X, 3),
    'Cartesian_Y_Mpc': np.round(Y, 3),
    'Cartesian_Z_Mpc': np.round(Z, 3),
    'Distance_Mpc': np.round(dist_clean, 3),
    'Stellar_Mass_Msun': np.round(M_compiled_kg / SOLAR_MASS_KG, 2),
    'Evaluated_H_local': np.round(galaxy_hlocal, 3)
})

csv_filename = 'galaxy_hlocal_evaluations.csv'
df_out.to_csv(csv_filename, index=False)
print(f"File Output Saved: '{csv_filename}' ({len(df_out)} galaxies evaluated).")

# --- Step 5: Enhanced Visual Plotting ---
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(10.5, 9), dpi=120)

# Extract Z=0 equatorial plane slice
z_slice_idx = RES // 2
h_slice = H_local_grid[:, :, z_slice_idx]

# Linear colormap stretch spanning true field bounds
extent = [-GRID_LIM, GRID_LIM, -GRID_LIM, GRID_LIM]
im = ax.imshow(h_slice.T, origin='lower', extent=extent, cmap='inferno', vmin=67.0, vmax=77.0)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label(r'Local Expansion Rate $H_{local}(\vec{r})$ [km/s/Mpc]', fontsize=11)

# Contour Overlay
contour_levels = [68.0, 70.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0,]
cs = ax.contour(h_slice.T, levels=contour_levels, colors='cyan', linewidths=0.7, alpha=0.6, extent=extent)
ax.clabel(cs, inline=True, fontsize=8, fmt='%.1f')

# Overlay Galaxies (|Z| < 2 Mpc)
z_mask = np.abs(Z) < 2.0
ax.scatter(X[z_mask], Y[z_mask], c='cyan', s=16, alpha=0.8, edgecolors='none', label='CNG Galaxies (|Z| < 2 Mpc)')

# Structure Annotations
ax.plot(0, 0, 'r+', markersize=12, markeredgewidth=2, label='Local Group (Origin)')
ax.annotate('Local Group', xy=(0, 0), xytext=(0.5, -0.8), color='white', fontsize=9, fontweight='bold')

# Labels & Aesthetics
ax.set_title('Substrate Logistics Local Expansion Field ($H_{local}$)\nDensity-Driven Pointer Recycling Field (CNG Catalog)', fontsize=12, pad=12)
ax.set_xlabel('Cartesian X [Mpc]', fontsize=11)
ax.set_ylabel('Cartesian Y [Mpc]', fontsize=11)
ax.legend(loc='upper right', framealpha=0.3, fontsize=9)
ax.grid(True, linestyle='--', alpha=0.2)

# Apply tight_layout BEFORE adding the inset axes
plt.tight_layout()

# Inset Histogram showing distribution of galaxy H_local values
ax_inset = fig.add_axes([0.15, 0.15, 0.22, 0.18])
ax_inset.hist(galaxy_hlocal, bins=25, color='orangered', edgecolor='black', alpha=0.8)
ax_inset.set_title(r'Galaxy $H_{local}$ Dist.', fontsize=8)
ax_inset.set_xlabel('km/s/Mpc', fontsize=7)
ax_inset.set_ylabel('Count', fontsize=7)
ax_inset.tick_params(labelsize=7)

plt.savefig('substrate_logistics_expansion_field.png', dpi=300)
plt.show()

# Display sample of evaluated CSV data
print("\nSample Data Output (Top Galaxies by H_local):")
print(df_out.sort_values('Evaluated_H_local', ascending=False).head(10).to_string(index=False))