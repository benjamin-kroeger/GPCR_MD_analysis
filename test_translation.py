import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Sample data
data = pd.DataFrame(np.random.rand(10, 10),
                    index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'],
                    columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])

def create_split_heatmap(data, rows_between_sections, cols_between_sections):
    fig, ax = plt.subplots(figsize=(8, 8))

    # Create a mask to separate sections
    mask = np.zeros_like(data)
    for i in range(data.shape[0]):
        if i % (rows_between_sections + 1) == 0 and i != 0:
            mask[i, :] = True
    for j in range(data.shape[1]):
        if j % (cols_between_sections + 1) == 0 and j != 0:
            mask[:, j] = True

    # Plot the heatmap with seaborn, using the mask to create white spaces
    sns.heatmap(data, cmap='coolwarm', annot=False, cbar=True, mask=mask, ax=ax)

    # Adjust x and y tick positions
    ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
    ax.set_xticklabels(data.columns, minor=False, rotation=45, ha="right")
    ax.set_yticklabels(data.index, minor=False)

    # Remove minor ticks and labels for cleaner appearance
    ax.tick_params(axis='both', which='both', length=0)

    plt.show()

# Set the number of rows and columns between the sections
rows_between_sections = 1
cols_between_sections = 1

create_split_heatmap(data, rows_between_sections, cols_between_sections)
