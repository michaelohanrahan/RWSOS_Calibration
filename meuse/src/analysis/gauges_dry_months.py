import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

# load example data
ds = xr.open_dataset(r'p:\11209265-grade2023\wflow\wflow_meuse_julia\wflow_meuse_per_catchment_N\_output\ds_obs_model_combined.nc')
obs = ds.sel(runs='Obs.')

# compute 7-day rolling mean
obs_7q = obs.rolling(time=7*24).mean()

# get the monthly min
monthly_obs = obs_7q.resample(time='M').min('time')

# # identify low flow months
# low_flow_months = monthly_obs.groupby('wflow_id').apply(
#     lambda x: x['Q'].groupby('time.month').mean().argmin(dim='month')
# )


# Extract the wflow_id names
wflow_ids = monthly_obs.wflow_id.values

# Create a DataFrame to store the results
results = []

# Loop through each wflow_id and year to find the month with the minimum discharge
for wflow_id in wflow_ids:
    wflow_data = monthly_obs.sel(wflow_id=wflow_id)
    for year in range(wflow_data['time.year'].min().item(), wflow_data['time.year'].max().item() + 1):
        year_data = wflow_data.sel(time=str(year))
        min_month = year_data['Q'].argmin().item()  # 'Q' is assumed to be the discharge variable
        min_month_date = pd.to_datetime(year_data.time.values[min_month])
        results.append([wflow_id, year, min_month_date.month])

# Convert the results to a DataFrame for easy viewing
results_df = pd.DataFrame(results, columns=['wflow_id', 'Year', 'Month'])
print(results_df)


# Count occurrences of each (Year, Month) combination
counts = results_df.groupby(['Year', 'Month']).size().reset_index(name='counts')

# Create scatter plot
plt.figure(figsize=(14, 8))
scatter = plt.scatter(counts['Year'], counts['Month'], s=counts['counts']*50, alpha=0.6, edgecolors="w", linewidth=0.5)

# Add labels to indicate the count
for i, row in counts.iterrows():
    plt.text(row['Year'], row['Month'], row['counts'], fontsize=9, ha='center', va='center')

plt.xlabel('Year')
plt.ylabel('Month')
plt.title('Frequency of Lowest Discharge Months per Year')
plt.grid(True)
plt.show()