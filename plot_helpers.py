from modules import *

def plot_counts_ordered(counts):
    items = sorted(counts.items(), key=lambda x: x[1], reverse=True)

    keys = [k for k, v in items]
    vals = [v for k, v in items]

    plt.figure(figsize=(18, 6))  
    plt.bar(range(len(vals)), vals)

    plt.xticks(range(len(keys)), keys, rotation='vertical', fontsize=6)

    plt.tight_layout()
    plt.show()