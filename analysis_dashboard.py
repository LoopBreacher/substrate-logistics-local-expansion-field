import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def generate_analysis_dashboard(csv_path='galaxy_hlocal_evaluations.csv', output_img='substrate_logistics_analysis_dashboard.png'):
    # Load evaluated CSV data
    df = pd.read_csv(csv_path)

    # Define spatial distance volume bins
    df['Dist_Bin'] = pd.cut(
        df['Distance_Mpc'], 
        bins=[0, 2, 4, 6, 8, 12], 
        labels=['0-2 Mpc', '2-4 Mpc', '4-6 Mpc', '6-8 Mpc', '8-12 Mpc']
    )

    # Set dark theme for publication aesthetic
    plt.style.use('dark_background')
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)

    # ---------------------------------------------------------
    # Panel 1: Distribution of H_local with Empirical Benchmarks
    # ---------------------------------------------------------
    ax1 = axes[0, 0]
    sns.histplot(df['Evaluated_H_local'], kde=True, ax=ax1, color='cyan', bins=30, stat="density", alpha=0.6)
    ax1.axvline(67.42, color='lime', linestyle='--', linewidth=2, label=r'Planck CMB Baseline ($67.42$)')
    ax1.axvline(73.04, color='magenta', linestyle='--', linewidth=2, label=r'SHOES 2022 Baseline ($73.04$)')
    ax1.axvline(df['Evaluated_H_local'].mean(), color='yellow', linestyle='-', linewidth=2, label=f"Catalog Mean ({df['Evaluated_H_local'].mean():.2f})")
    ax1.set_title(r'Distribution of Evaluated $H_{\mathrm{local}}(\vec{r})$ across 451 CNG Galaxies', fontsize=12, pad=10)
    ax1.set_xlabel(r'$H_{\mathrm{local}}(\vec{r})$ [km/s/Mpc]', fontsize=11)
    ax1.set_ylabel('Density', fontsize=11)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.grid(True, linestyle=':', alpha=0.3)

    # ---------------------------------------------------------
    # Panel 2: Radial Profile vs. Distance from Origin
    # ---------------------------------------------------------
    ax2 = axes[0, 1]
    ax2.scatter(df['Distance_Mpc'], df['Evaluated_H_local'], c=df['Evaluated_H_local'], cmap='plasma', alpha=0.7, edgecolors='none', s=25)
    
    # 30-galaxy rolling mean profile
    df_sorted_dist = df.sort_values('Distance_Mpc')
    rolling_mean = df_sorted_dist['Evaluated_H_local'].rolling(window=30, min_periods=5).mean()
    ax2.plot(df_sorted_dist['Distance_Mpc'], rolling_mean, color='cyan', linewidth=2.5, label='30-Galaxy Rolling Mean')
    ax2.axhline(67.42, color='lime', linestyle=':', label='Void Baseline (67.42)')
    ax2.axhline(73.04, color='magenta', linestyle=':', label='SHOES Value (73.04)')
    ax2.set_title(r'Radial Profile of $H_{\mathrm{local}}(\vec{r})$ vs. Distance from Origin', fontsize=12, pad=10)
    ax2.set_xlabel('Distance from Local Group [Mpc]', fontsize=11)
    ax2.set_ylabel(r'$H_{\mathrm{local}}(\vec{r})$ [km/s/Mpc]', fontsize=11)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, linestyle=':', alpha=0.3)

    # ---------------------------------------------------------
    # Panel 3: H_local vs Stellar Mass (Log Scale)
    # ---------------------------------------------------------
    ax3 = axes[1, 0]
    sc = ax3.scatter(df['Stellar_Mass_Msun'], df['Evaluated_H_local'], c=df['Distance_Mpc'], cmap='viridis', alpha=0.7, s=25)
    ax3.set_xscale('log')
    ax3.set_title(r'$H_{\mathrm{local}}(\vec{r})$ vs. Compiled Stellar Mass', fontsize=12, pad=10)
    ax3.set_xlabel(r'Stellar Mass [$M_\odot$] (Log Scale)', fontsize=11)
    ax3.set_ylabel(r'$H_{\mathrm{local}}(\vec{r})$ [km/s/Mpc]', fontsize=11)
    ax3.grid(True, linestyle=':', alpha=0.3)
    cbar3 = fig.colorbar(sc, ax=ax3)
    cbar3.set_label('Distance [Mpc]', fontsize=9)

    # ---------------------------------------------------------
    # Panel 4: Boxplot across Distance Volume Shells
    # ---------------------------------------------------------
    ax4 = axes[1, 1]
    sns.boxplot(x='Dist_Bin', y='Evaluated_H_local', data=df, ax=ax4, hue='Dist_Bin', palette='rocket', legend=False, boxprops=dict(alpha=0.7))
    ax4.axhline(67.42, color='lime', linestyle='--', linewidth=1.5, label='Planck 2018 (67.42)')
    ax4.axhline(73.04, color='magenta', linestyle='--', linewidth=1.5, label='SHOES 2022 (73.04)')
    ax4.set_title(r'$H_{\mathrm{local}}$ Variation across Volume Shells', fontsize=12, pad=10)
    ax4.set_xlabel('Spatial Distance Shell', fontsize=11)
    ax4.set_ylabel(r'$H_{\mathrm{local}}(\vec{r})$ [km/s/Mpc]', fontsize=11)
    ax4.legend(loc='upper right', fontsize=9)
    ax4.grid(True, linestyle=':', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_img, dpi=300)
    plt.close()
    print(f"Dashboard saved as '{output_img}'.")

if __name__ == '__main__':
    generate_analysis_dashboard()