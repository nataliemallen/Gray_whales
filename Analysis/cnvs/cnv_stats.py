import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

results_dir = '/scratch/gautschi/allen715/2026_whales/cnvs/results'
output_dir = '/scratch/gautschi/allen715/2026_whales/cnvs/stats'
os.makedirs(output_dir, exist_ok=True)

pop1_samples = [l.strip() for l in open('/scratch/gautschi/allen715/2026_whales/lists/east.txt')]
pop2_samples = [l.strip() for l in open('/scratch/gautschi/allen715/2026_whales/lists/west_ds.txt')]
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


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR-adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]   # enforce monotonicity
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    return adj


def run_stats(cnv):
    res = {}

    # (1) overall CNV abundance per individual: Mann-Whitney U
    counts = cnv.groupby(['sample','population']).size().reset_index(name='count')
    enp = counts[counts.population=='ENP']['count']
    wnp = counts[counts.population=='WNP']['count']
    res['abundance_U'], res['abundance_p'] = stats.mannwhitneyu(enp, wnp, alternative='two-sided')
    res['abundance_median_ENP'] = enp.median()
    res['abundance_median_WNP'] = wnp.median()

    # (2) duplication vs deletion proportions: chi-square
    type_table = cnv.groupby(['population','type']).size().unstack(fill_value=0)
    res['type_table'] = type_table
    chi2_y, p_y, dof, _ = stats.chi2_contingency(type_table)                 # Yates-corrected (2x2 default)
    chi2_n, p_n, _, _   = stats.chi2_contingency(type_table, correction=False)
    res['type_chi2_yates'], res['type_p_yates'] = chi2_y, p_y
    res['type_chi2_nocorr'], res['type_p_nocorr'] = chi2_n, p_n
    res['type_dof'] = dof

    # (3) per-region differentiation: Fisher's exact + BH-FDR
    regions = cnv.groupby(['type','scaffold','start','end']).size().reset_index()
    pvals = []
    for _, r in regions.iterrows():
        mask = ((cnv.type==r.type) & (cnv.scaffold==r.scaffold) &
                (cnv.start==r.start) & (cnv.end==r.end))
        samples = set(cnv.loc[mask,'sample'])
        enp_with = len([s for s in pop1_samples if s in samples])
        wnp_with = len([s for s in pop2_samples if s in samples])
        table = [[enp_with, len(pop1_samples)-enp_with],
                 [wnp_with, len(pop2_samples)-wnp_with]]
        _, p = stats.fisher_exact(table)
        pvals.append(p)
    pvals = np.array(pvals)
    padj  = bh_fdr(pvals)
    res['fisher_n_regions']   = len(pvals)
    res['fisher_n_sig_raw']   = int(np.sum(pvals < 0.05))
    res['fisher_n_sig_fdr']   = int(np.sum(padj < 0.05))
    res['fisher_min_p']       = float(pvals.min()) if pvals.size else float('nan')
    res['fisher_min_padj']    = float(padj.min()) if padj.size else float('nan')
    return res


def main():
    for bin_size in ['10kb','100kb']:
        cnv = load_cnv_data(bin_size)
        create_plot(cnv, f"{output_dir}/cnv_{bin_size}.png")
        r = run_stats(cnv)

        lines = []
        lines.append(f"===== CNV statistics: {bin_size} =====")
        lines.append(f"total CNV calls: {len(cnv)}   (ENP {(cnv.population=='ENP').sum()}, "
                     f"WNP {(cnv.population=='WNP').sum()})")
        lines.append("")
        lines.append("(1) Overall CNV abundance per individual  -- Mann-Whitney U:")
        lines.append(f"    U = {r['abundance_U']:.1f}   p = {r['abundance_p']:.4g}   "
                     f"(median ENP {r['abundance_median_ENP']}, WNP {r['abundance_median_WNP']})")
        lines.append("")
        lines.append("(2) Duplication vs deletion proportion  -- chi-square:")
        lines.append(f"    contingency (rows=pop, cols=type):\n{r['type_table'].to_string()}")
        lines.append(f"    Yates-corrected : chi2 = {r['type_chi2_yates']:.2f}  p = {r['type_p_yates']:.4g}  df = {r['type_dof']}")
        lines.append(f"    uncorrected     : chi2 = {r['type_chi2_nocorr']:.2f}  p = {r['type_p_nocorr']:.4g}")
        lines.append("")
        lines.append("(3) Per-region differentiation  -- Fisher's exact + BH-FDR:")
        lines.append(f"    regions tested       : {r['fisher_n_regions']}")
        lines.append(f"    significant (raw)    : {r['fisher_n_sig_raw']}  (min p = {r['fisher_min_p']:.3g})")
        lines.append(f"    significant (FDR<.05): {r['fisher_n_sig_fdr']}  (min adj p = {r['fisher_min_padj']:.3g})")
        text = "\n".join(lines)

        print(text + "\n")
        with open(f"{output_dir}/summary_{bin_size}.txt","w") as f:
            f.write(text + "\n")


if __name__ == "__main__":
    main()