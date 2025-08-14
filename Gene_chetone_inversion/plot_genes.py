import matplotlib.pyplot as plt
from adjustText import adjust_text

# Read gene data from the file
genes = []
with open("end_genes.txt") as f:
    for line in f:
        parts = line.strip().split()
        gene_name = parts[0]
        start_pos = float(parts[1])
        end_pos = float(parts[2])
        genes.append((gene_name, start_pos, end_pos))

# Create the plot
fig, ax = plt.subplots(figsize=(12, 1))

# Draw the genome line
ax.plot([min([gene[1] for gene in genes]), max([gene[2] for gene in genes])], [0, 0], color='black', lw=2)

# Create a list to store the text objects
texts = []

# Add ticks and gene numbers
for i, (gene, start, end) in enumerate(genes):
    mid_pos = (start + end) / 2
    ax.plot(mid_pos, 0, 'k|', markersize=10)  # Gene as tick at mid position
    # Create the text object (gene number)
    text = ax.text(mid_pos, 0.01, str(i), ha='center', va='bottom', fontsize=6)  # Gene number above the tick
    texts.append(text)

# Adjust layout using adjustText to avoid overlap
adjust_text(texts, ax=ax, only_move={'points': 'xy', 'text': 'xy'},
            expand_points=(1.05, 1.05), expand_text=(1.1, 1.1), lim=50, precision=0.01)

# Adjust layout
ax.set_yticks([])
ax.set_xticks([])

# Remove frame around plot
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Save the plot to a file
plt.savefig("gene_plot.png", bbox_inches='tight', dpi=300)

plt.close()
