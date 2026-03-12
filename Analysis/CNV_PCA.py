import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os

results_dir = '/scratch/gautschi/allen715/whale/new_analysis/long_updated/calls'
bin_size = '100kb'

pop1_samples = [l.strip() for l in open('east_bams.txt')]
pop2_samples = [l.strip() for l in open('west_bams_33.txt')]
all_samples = pop1_samples + pop2_samples


def parse_file(file_path, sample):
    df = pd.read_csv(file_path, sep='\t', header=None)
    df.columns = ['type','region','length','rd','p1','p2','p3','p4','a','b','c']
    df[['scaffold','coords']] = df['region'].str.split(':', expand=True)
    df[['start','end']] = df['coords'].str.split('-', expand=True)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['sample'] = sample
    return df[['type','scaffold','start','end','sample']]


calls = []
for s in all_samples:
    f = os.path.join(results_dir,f"{s}.{bin_size}.calls.tsv")
    calls.append(parse_file(f,s))

cnv = pd.concat(calls)

regions = cnv.groupby(['type','scaffold','start','end']).size().reset_index()
regions['id'] = regions.index

sample_region = []
for _, r in regions.iterrows():
    mask = (
        (cnv.type==r.type) &
        (cnv.scaffold==r.scaffold) &
        (cnv.start==r.start) &
        (cnv.end==r.end)
    )
    for s in cnv.loc[mask,'sample'].unique():
        sample_region.append({'sample':s,'region':r.id})

sr = pd.DataFrame(sample_region)

matrix = pd.crosstab(sr['sample'], sr['region'])

pca = PCA(n_components=2)
coords = pca.fit_transform(matrix)

plt.figure(figsize=(8,6))

for samples,label in [(pop1_samples,'East'),(pop2_samples,'WNP')]:
    idx = matrix.index.isin(samples)
    plt.scatter(coords[idx,0],coords[idx,1],label=label)

plt.legend()
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig("cnv_pca.png",dpi=300)
