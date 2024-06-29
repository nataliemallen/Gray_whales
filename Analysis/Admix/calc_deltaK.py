# %%
import numpy as np
import pandas as pd

# %%
# Load the log likelihood values
data = pd.read_csv('log_likelihoods.txt', sep=' ', header=None, names=['Run', 'LogLikelihood'])

# %%
# Function to extract K value and run number
def extract_k_run(run_identifier):
    parts = run_identifier.split('_')
    k_value = int(parts[1][1:])  # Extract number after 'K'
    run_number = int(parts[2][3:])  # Extract number after 'run'
    return k_value, run_number

# %%
# Apply the function to extract K and run number
data[['K', 'RunNumber']] = data['Run'].apply(lambda x: pd.Series(extract_k_run(x)))

# %%
# Calculate mean and standard deviation for each K
stats = data.groupby('K')['LogLikelihood'].agg(['mean', 'std']).reset_index()

# %%
# Calculate the difference in mean log likelihoods for successive K values
stats['DeltaL'] = stats['mean'].diff().abs()

# %%
# Shift the DeltaL column to align with the current K and calculate DeltaK
stats['DeltaK'] = stats['DeltaL'].shift(-1) / stats['std']

# %%
# Display results
print(stats)

# Save the results to a file
stats.to_csv('deltaK_results.csv', index=False)

#   K          mean            std        DeltaL        DeltaK
#0  1 -1.233738e+08       0.000001           NaN  2.183049e+12
#1  2 -1.202521e+08       1.774507  3.121663e+06  1.402293e+06
#2  3 -1.177637e+08   84482.302936  2.488378e+06  2.369075e+01
#3  4 -1.157623e+08  137382.456510  2.001449e+06  1.345976e+01
#4  5 -1.139131e+08  135268.160114  1.849135e+06           NaN

# %%
#K=2 is highest deltaK; best supported
#3311553 sites


