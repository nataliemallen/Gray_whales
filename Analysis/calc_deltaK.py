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

# %%
#K=2 is highest deltaK; best supported
#37214 sites

