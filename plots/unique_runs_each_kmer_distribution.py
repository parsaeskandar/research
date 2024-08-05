import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


kmer_sizes = [10, 30, 100, 300]

for kmer_size in kmer_sizes:
    # Specify the path to your text file
    file_path = f'/Users/seeskand/Documents/pangenome-index/chr19/{kmer_size}_tags.txt'

    # Loading the data into a pandas dataframe
    df = pd.read_csv(file_path, sep='\t', header=None, names=['runs', 'uniques', 'run/unique'])
    df_filtered = df
    # Plotting the distribution
    plt.figure()
    sns.histplot(df_filtered['run/unique'], kde=True)
    # plt.xlim(0, 10)
    plt.title(f'Distribution of run/unique for kmer size {kmer_size}')
    plt.xlabel('run/unique Value')
    plt.ylabel('Frequency')
    plt.grid(True)
    # plt.show()
    plt.savefig(f'/Users/seeskand/Documents/pangenome-index/chr19/plots/{kmer_size}_tags_distribution.png', dpi=300)


