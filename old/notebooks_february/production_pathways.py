import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
file_path = 'eu_industry_prod_scenarios.csv'  # Change this to the path of your CSV file
df = pd.read_csv(file_path)

# Get unique sector values
sectors = df['sector'].unique()

# Create a figure with subplots for each unique sector (excluding 'chlorine')
sectors_to_plot = [sector for sector in sectors if sector != 'chlorine']
fig, axes = plt.subplots(1, len(sectors_to_plot), figsize=(len(sectors_to_plot)*3, 6), sharex=True)

# If there's only one sector, axes will not be an array, so we ensure it's iterable
if len(sectors_to_plot) == 1:
    axes = [axes]

# Iterate over each unique sector to create the corresponding subplot
for i, sector in enumerate(sectors_to_plot):
    sector_df = df[df['sector'] == sector]
    
    # Plot the 'regain' and 'deindustrial' columns against 'year'
    axes[i].plot(sector_df['year'], sector_df['regain'], label='Regain', color='blue')
    axes[i].plot(sector_df['year'], sector_df['deindustrial'], label='Deindustrial', color='red')
    
    # Set titles and labels for each subplot
    axes[i].set_title(f'{sector.capitalize()}', fontsize=14)
    if sector == 'steel':
        axes[i].set_ylabel('Mt steel/yr', fontsize=12)
    elif sector == 'cement':
        axes[i].set_ylabel('Mt cement/yr', fontsize=12)
    elif sector == 'ammonia':
        axes[i].set_ylabel('Mt NH3/yr', fontsize=12)
    elif sector == 'hvc':
        axes[i].set_ylabel('Mt HVC/yr', fontsize=12)
    elif sector == 'methanol':
        axes[i].set_ylabel('Mt methanol/yr', fontsize=12)
        
    axes[i].set_ylim(bottom=0)
    axes[i].grid(True, linestyle='--')  # Enable grid on each subplot
    axes[0].legend()

# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig("graphs/production_projections.png", bbox_inches="tight")


# Show the plot
plt.show()
