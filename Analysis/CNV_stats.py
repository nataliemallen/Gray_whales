import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

results_dir = '/scratch/gautschi/allen715/whale/new_analysis/long_updated/calls'
output_dir = '/scratch/gautschi/allen715/whale/new_analysis/cnv_analysis'
os.makedirs(output_dir, exist_ok=True)

pop1_samples = [l.strip() for l in open('/scratch/gautschi/allen715/whale/new_analysis/east_bams.txt')]
pop2_samples = [l.strip() for l in open('/scratch/gautschi/allen715/whale/new_analysis/west_bams_33.txt')]
all_samples = pop1_samples + pop2_samples


def parse_cnv_file(file_path, sample):
    df = pd.read_csv(file_path, sep='\t', header=None)
    df.columns = ['type','region','length','mean_rd','p1','p2','p3','p4','q0','dG','start_pos']
    df[['scaffold','coords']] = df['region'].str.split(':', expand=True)
    df[['start','end']] = df['coords'].str.split('-', expand=True)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['sample'] = sample
    return df[['type','scaffold','start','end','length','sample']]


def load_cnv_data(bin_size):
    dfs = []
    for s in all_samples:
        f = os.path.join(results_dir, f"{s}.{bin_size}.calls.tsv")
        if os.path.exists(f):
            dfs.append(parse_cnv_file(f, s))
    cnv = pd.concat(dfs, ignore_index=True)
    cnv['population'] = np.where(cnv['sample'].isin(pop1_samples),'ENP','WNP')
    return cnv


def create_plot(cnv, outfile):
    fig, axes = plt.subplots(2,2,figsize=(10,8))

    counts = cnv.groupby(['sample','population']).size().reset_index(name='count')

    axes[0,0].boxplot([
        counts[counts.population=='ENP']['count'],
        counts[counts.population=='WNP']['count']
    ])
    axes[0,0].set_xticklabels(['ENP','WNP'])

    axes[0,1].hist(np.log10(cnv['length']), bins=30)

    types = cnv.groupby(['population','type']).size().unstack(fill_value=0)
    types.plot(kind='bar', ax=axes[1,0])

    chrom = cnv.groupby(['scaffold','population']).size().unstack(fill_value=0)
    chrom.head(20).plot(kind='bar', ax=axes[1,1])

    plt.tight_layout()
    plt.savefig(outfile,dpi=300)
    plt.close()


def run_stats(cnv):
    counts = cnv.groupby(['sample','population']).size().reset_index(name='count')

    enp = counts[counts.population=='ENP']['count']
    wnp = counts[counts.population=='WNP']['count']

    w_stat, w_p = stats.mannwhitneyu(enp,wnp)

    type_table = cnv.groupby(['population','type']).size().unstack(fill_value=0)
    chi2, chi_p, _, _ = stats.chi2_contingency(type_table)

    regions = cnv.groupby(['type','scaffold','start','end']).size().reset_index()

    fisher = []
    for _, r in regions.iterrows():
        mask = (
            (cnv.type==r.type) &
            (cnv.scaffold==r.scaffold) &
            (cnv.start==r.start) &
            (cnv.end==r.end)
        )
        samples = set(cnv.loc[mask,'sample'])

        enp_with = len([s for s in pop1_samples if s in samples])
        wnp_with = len([s for s in pop2_samples if s in samples])

        table = [
            [enp_with, len(pop1_samples)-enp_with],
            [wnp_with, len(pop2_samples)-wnp_with]
        ]

        _, p = stats.fisher_exact(table)
        fisher.append(p)

    fisher = np.array(fisher)
    sig = np.sum(fisher < 0.05)

    return w_stat, w_p, chi2, chi_p, sig


def main():
    for bin_size in ['10kb','100kb']:
        cnv = load_cnv_data(bin_size)
        create_plot(cnv, f"{output_dir}/cnv_{bin_size}.png")
        stats = run_stats(cnv)

        with open(f"{output_dir}/summary_{bin_size}.txt","w") as f:
            f.write(str(stats))


if __name__ == "__main__":
    main()
