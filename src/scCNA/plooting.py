import seaborn as sns
import matplotlib.pyplot as plt

def plot_cna_heatmap(assignments, title="CNA Assignment Heatmap"):
    plt.figure(figsize=(12, 6))
    sns.heatmap(assignments, cmap="coolwarm", cbar=True)
    plt.title(title)
    plt.xlabel("CNAs")
    plt.ylabel("Cells")
    plt.tight_layout()
    plt.show()
