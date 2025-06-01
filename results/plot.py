import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


os.chdir(r'C:\Phd\CUDA test\Test\test 1\EV\EV_delivery\code\code_JD_data\results')

base_file = "jd200_1" # base file name
file_path = "results_{}.csv".format(base_file)
df = pd.read_csv(file_path)





# Define the delivery targets to filter
delivery_targets = [5, 10, 15]
filtered_df = df[df['Total delivery'].isin(delivery_targets)]

# Set style and color palette
sns.set(style="whitegrid")
# palette = sns.color_palette("tab10")
palette = {
    "Optimal": "#1f77b4",  # Blue
    "CSA": "#2ca02c",      # Green 
    "NDF": "#ff7f0e",      # Orange
    "EDF": "#d62728"       # Red
}

# Plot Total Successful Delivery
plt.figure(figsize=(8, 6))
sns.barplot(data=filtered_df, x='Total delivery', y='Total successful delivery', hue='Methods', palette=palette)
plt.xlabel('Total Delivery', fontsize=22)
plt.ylabel('Total Successful Delivery', fontsize=22)
plt.legend(title='Methods', fontsize=22, title_fontsize=22, loc='lower right')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig('Successful_Delivery_small.pdf', dpi=300)
plt.close()

# Plot Avg. Cost per Successful Delivery
plt.figure(figsize=(8, 6))
sns.barplot(data=filtered_df, x='Total delivery', y='Avg. cost per successful delivery', hue='Methods', palette=palette)
plt.xlabel('Total Delivery', fontsize=22)
plt.ylabel('Avg. cost per delivery ($)', fontsize=22)
plt.legend(title='Methods', fontsize=22, title_fontsize=22, loc='lower right')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig('Avg_cost_small.pdf', dpi=300)
plt.close()

# Plot Elapsed Time (log scale)
plt.figure(figsize=(8, 6))
sns.barplot(data=filtered_df, x='Total delivery', y='Elapsed time (ms)', hue='Methods', palette=palette)
plt.yscale('log')
plt.xlabel('Total Delivery', fontsize=22)
plt.ylabel('Elapsed time (ms)', fontsize=22)
plt.legend(fontsize=22, title_fontsize=22, loc='upper left')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig('Elapsed_time_small.pdf', dpi=300)
plt.close()





## Scalability Experiments
# Filter the dataset for the new conditions
nC = 5
delivery2ev_ratio = 5

filtered_large_df = df[
    (df['Total delivery'] > 15) &
    (df['Total CP'] == nC) &
    (df['delivery2ev_ratio'] == delivery2ev_ratio)
]


# Plot Total Successful Delivery
plt.figure(figsize=(8, 6))
sns.barplot(data=filtered_large_df, x='Total delivery', y='Total successful delivery', hue='Methods', palette=palette)
plt.xlabel('Total Delivery', fontsize=22)
plt.ylabel('Total Successful Delivery', fontsize=22)
plt.legend(title='Methods', fontsize=22, title_fontsize=22, loc='upper left')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'Successful_Delivery_large_CP{nC}_R{delivery2ev_ratio}.pdf', dpi=300)
plt.close()

# Plot Avg. Cost per Successful Delivery
plt.figure(figsize=(8, 6))
sns.barplot(data=filtered_large_df, x='Total delivery', y='Avg. cost per successful delivery', hue='Methods', palette=palette)
plt.xlabel('Total Delivery', fontsize=22)
plt.ylabel('Avg. cost per delivery ($)', fontsize=22)
plt.legend(title='Methods', fontsize=22, title_fontsize=22, loc='lower right')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'Avg_cost_large_CP{nC}_R{delivery2ev_ratio}.pdf', dpi=300)
plt.close()

# # Plot Elapsed Time (log scale)
# plt.figure(figsize=(8, 6))
# sns.barplot(data=filtered_large_df, x='Total delivery', y='Elapsed time (ms)', hue='Methods', palette=palette)
# plt.yscale('log')
# plt.xlabel('Total Delivery', fontsize=22)
# plt.ylabel('Elapsed time (ms)', fontsize=22)
# plt.legend(title='Methods', fontsize=22, title_fontsize=22, loc='lower right')
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(True)
# plt.tight_layout()
# plt.savefig(f'Elapsed_time_large_CP{nC}_R{delivery2ev_ratio}.pdf', dpi=300)
# plt.close()


plt.figure(figsize=(8, 6))
sns.lineplot(data=filtered_large_df, x='Total delivery', y='Elapsed time (ms)', hue='Methods', marker='o', palette=palette)
# plt.yscale('log')
plt.xlabel('Total Delivery', fontsize=22)
plt.ylabel('Elapsed time (ms)', fontsize=22)
plt.legend(fontsize=22, title_fontsize=22, loc='lower right')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'Elapsed_time_large_CP{nC}_R{delivery2ev_ratio}.pdf', dpi=300)
plt.close()




## Scalability Experiments
# Filter the dataset for the new conditions
nC = 25
delivery2ev_ratio = 5

filtered_large_df = df[
    (df['Total delivery'] > 15) &
    (df['Total CP'] == nC) &
    (df['delivery2ev_ratio'] == delivery2ev_ratio)
]


# Plot Total Successful Delivery
plt.figure(figsize=(8, 6))
sns.barplot(data=filtered_large_df, x='Total delivery', y='Total successful delivery', hue='Methods', palette=palette)
plt.xlabel('Total Delivery', fontsize=22)
plt.ylabel('Total Successful Delivery', fontsize=22)
plt.legend(title='Methods', fontsize=22, title_fontsize=22, loc='upper left')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'Successful_Delivery_large_CP{nC}_R{delivery2ev_ratio}.pdf', dpi=300)
plt.close()

# Plot Avg. Cost per Successful Delivery
plt.figure(figsize=(8, 6))
sns.barplot(data=filtered_large_df, x='Total delivery', y='Avg. cost per successful delivery', hue='Methods', palette=palette)
plt.xlabel('Total Delivery', fontsize=22)
plt.ylabel('Avg. cost per delivery ($)', fontsize=22)
plt.legend(title='Methods', fontsize=22, title_fontsize=22, loc='lower right')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'Avg_cost_large_CP{nC}_R{delivery2ev_ratio}.pdf', dpi=300)
plt.close()

# # Plot Elapsed Time (log scale)
# plt.figure(figsize=(8, 6))
# sns.barplot(data=filtered_large_df, x='Total delivery', y='Elapsed time (ms)', hue='Methods', palette=palette)
# plt.yscale('log')
# plt.xlabel('Total Delivery', fontsize=22)
# plt.ylabel('Elapsed time (ms)', fontsize=22)
# plt.legend(title='Methods', fontsize=22, title_fontsize=22, loc='lower right')
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.grid(True)
# plt.tight_layout()
# plt.savefig(f'Elapsed_time_large_CP{nC}_R{delivery2ev_ratio}.pdf', dpi=300)
# plt.close()


plt.figure(figsize=(8, 6))
sns.lineplot(data=filtered_large_df, x='Total delivery', y='Elapsed time (ms)', hue='Methods', marker='o', palette=palette)
# plt.yscale('log')
plt.xlabel('Total Delivery', fontsize=22)
plt.ylabel('Elapsed time (ms)', fontsize=22)
plt.legend(fontsize=22, title_fontsize=22, loc='lower right')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'Elapsed_time_large_CP{nC}_R{delivery2ev_ratio}.pdf', dpi=300)
plt.close()





## Vary the number of CP along X-axis
delivery2ev_ratio = 5
nD = 200

filtered_cp_df = df[
    (df['Total delivery'] == nD) &
    (df['delivery2ev_ratio'] == delivery2ev_ratio)
]

# Plot Avg. Cost per Successful Delivery vs Total CP
plt.figure(figsize=(10, 3.5))
sns.lineplot(data=filtered_cp_df, x='Total CP', y='Avg. cost per successful delivery', hue='Methods', marker='o', palette=palette)
plt.xlabel('Total CP', fontsize=22)
plt.ylabel('Avg. cost per \ndelivery ($)', fontsize=22)
plt.legend(fontsize=22, title_fontsize=22, loc='upper left')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.show
plt.savefig(f'Avg_Cost_varying_CP_D{nD}_R{delivery2ev_ratio}.pdf', dpi=300)
plt.close()





## Vary delivery2ev_ratio along X-axis
nC = 50
nD = 200

filtered_ratio_df = df[
    (df['Total delivery'] == nD) &
    (df['Total CP'] == nC)
]

# Plot Avg. Cost per Successful Delivery vs delivery2ev_ratio
plt.figure(figsize=(10, 3.5))
sns.lineplot(data=filtered_ratio_df, x='delivery2ev_ratio', y='Avg. cost per successful delivery', hue='Methods', marker='o', palette=palette)
plt.xlabel('Delivery to EV Ratio', fontsize=22)
plt.ylabel('Avg. cost per \ndelivery ($)', fontsize=22)
plt.legend(fontsize=22, title_fontsize=22, loc='upper left')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.show
plt.savefig(f'Avg_Cost_varying_Ratio_D{nD}_CP{nC}.pdf', dpi=300)
plt.close()






mean_costs = df.groupby("Methods")["Avg. cost per successful delivery"].mean().reset_index()
print(mean_costs)



mean_succ_del = df.groupby("Methods")["Total successful delivery"].mean().reset_index()
print(mean_succ_del)