import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read the CSV file
df = pd.read_csv('galaxy_hlocal_evaluations.csv')

print("--- DATAFRAME OVERVIEW ---")
print(df.info())
print("\n--- STATISTICAL SUMMARY ---")
stats = df[['Distance_Mpc', 'Stellar_Mass_Msun', 'Evaluated_H_local']].describe()
print(stats)

print("\n--- TOP 10 HIGHEST H_LOCAL GALAXIES ---")
print(df.sort_values('Evaluated_H_local', ascending=False)[['Galaxy_Name', 'Distance_Mpc', 'Stellar_Mass_Msun', 'Evaluated_H_local']].head(10))

print("\n--- TOP 10 LOWEST H_LOCAL GALAXIES ---")
print(df.sort_values('Evaluated_H_local', ascending=True)[['Galaxy_Name', 'Distance_Mpc', 'Stellar_Mass_Msun', 'Evaluated_H_local']].head(10))

# Distance-binned H_local statistics
df['Dist_Bin'] = pd.cut(df['Distance_Mpc'], bins=[0, 2, 4, 6, 8, 12], labels=['0-2 Mpc', '2-4 Mpc', '4-6 Mpc', '6-8 Mpc', '8-12 Mpc'])
print("\n--- H_LOCAL BY DISTANCE BIN ---")
bin_summary = df.groupby('Dist_Bin', observed=False)['Evaluated_H_local'].agg(['count', 'mean', 'median', 'std', 'min', 'max'])
print(bin_summary)