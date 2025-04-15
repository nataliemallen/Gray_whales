# %%
print("hello")

# %%
import csv
import itertools
import os

# Define the groups
groups = {
    "pink_west": [
        "038001_ER-14-0163",
        "038007_ER-16-0047",
        "038016_ER-16-0058",
        "038020_ER-16-0063",
        "038023_ER-16-0067",
        "038026_ER-16-0073",
        "038057_ER-14-0165"
    ],
    "east": [
        "038029_ER-17-0175",
        "038033_ER-17-0187",
        "038037_ER-17-0199",
        "038040_ER-17-0218",
        "038044_ER-17-0225",
        "038047_ER-17-0241",
        "038045_ER-17-0229"
    ],
    "blue_west": [
        "038006_ER-16-0046",
        "038018_ER-16-0061",
        "038021_ER-16-0065",
        "037997_ER-14-0151",
        "038025_ER-16-0072",
        "038059_ER-16-0059",
        "037997_ER-14-0159"
    ]
}

# Create output directory if it doesn't exist
output_dir = "pairwise_comparisons"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Generate all possible combinations of groups (including same group comparisons)
group_combinations = list(itertools.product(groups.keys(), groups.keys()))

# Process each group combination
for group1_name, group2_name in group_combinations:
    group1 = groups[group1_name]
    group2 = groups[group2_name]
    
    # Generate pairs differently based on whether we're comparing within a group or between groups
    if group1_name == group2_name:
        # Within-group comparison: use combinations to get unique pairs
        # This ensures each pair appears only once (A,B) but not (B,A)
        pairs = list(itertools.combinations(group1, 2))
    else:
        # Between-group comparison: use product to get all possible pairs
        pairs = list(itertools.product(group1, group2))
    
    # Create filename for the comparison
    filename = f"{output_dir}/{group1_name}_vs_{group2_name}.csv"
    
    # Write pairs to CSV file
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Individual1', 'Individual2'])  # Header
        writer.writerows(pairs)
    
    print(f"Created {filename} with {len(pairs)} comparisons")

print("All pairwise comparison files have been generated!")

# %%
import numpy as np
import matplotlib.pyplot as plt

# Data organization
categories = [
    'blue west vs blue west',
    'pink west vs pink west',
    'east vs east',
    'pink west vs east',
    'blue west vs east',
    'pink west vs blue west'
]

averages = [101300024, 56205165.5, 109722160, 128890794, 123010472, 120422953]
std_devs = [17295282.5, 34032863.7, 20078813, 8537261.26, 9689838.36, 9061411.61]

# Create figure and axis
fig, ax = plt.subplots(figsize=(12, 7))

# Bar positions
x_pos = np.arange(len(categories))

# Create bars
bars = ax.bar(x_pos, averages, yerr=std_devs, align='center',
              alpha=0.8, capsize=10, error_kw={'elinewidth': 2, 'capthick': 2})

# Add some visual improvements
ax.set_ylabel('Average Value', fontsize=14)
ax.set_title('Average Values with Standard Deviation by Group', fontsize=16)
ax.set_xticks(x_pos)
ax.set_xticklabels(categories, rotation=45, ha='right', fontsize=12)
ax.yaxis.grid(True)

# Format y-axis with commas for thousands
plt.gca().yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))

# Tight layout to ensure everything fits well
plt.tight_layout()

# Show the plot
plt.show()
