#!/usr/bin/env python3
"""
Aggregate all QC evidence into ONE master verdict table.
Pure stdlib (no pandas). Joins:
  - per-sample rows (qc01)           by sample name
  - ngsRelate relate.res             by index -> bam_list order -> name
  - PCAngsd inbreeding (per sample)  by bam_list order
  - cohort IBS matrix                by bam_list order (duplicate detection)

Verdict per sample:
  FAIL    -> evidence of mixture / contamination / duplicate (do not use as-is)
  REVIEW  -> technical outlier (depth/breadth/dup/error/insert) worth a look
  PASS    -> clean on every check
"""
import argparse, csv, os, statistics, sys
from math import isnan

# relatedness thresholds (kinship-coefficient scale)
KING_2ND   = 0.08838835     # >= this by KING => 2nd-degree-or-closer
THETA_UNREL= 0.04419417     # < this by model theta => effectively unrelated
RELATED_MANY_N = 5          # # of "KING-relative but theta-unrelated" partners to flag a mixture
# robust-outlier sensitivity
Z = 3.0                     # |robust z| beyond this = outlier
DUP_FRAC = 0.30             # IBS < DUP_FRAC * cohort-median => duplicate/near-identical

def med(xs):
    xs=[x for x in xs if x is not None and not (isinstance(x,float) and isnan(x))]
    return statistics.median(xs) if xs else None
def mad(xs, m=None):
    xs=[x for x in xs if x is not None and not (isinstance(x,float) and isnan(x))]
    if not xs: return None
    m = med(xs) if m is None else m
    d=[abs(x-m) for x in xs]
    return statistics.median(d) if d else None
def scale(xs, m=None):
    """robust scale = 1.4826*MAD, with stdev fallback when MAD==0."""
    xs=[x for x in xs if x is not None and not (isinstance(x,float) and isnan(x))]
    if len(xs)<2: return None
    md=mad(xs,m); s=1.4826*md if md else 0.0
    if not s:
        try: s=statistics.pstdev(xs)
        except: s=0.0
    return s or None
def rz(x, m, sc):
    if x is None or m is None or sc in (None,0): return 0.0
    return (x-m)/sc
def f(x):
    try: return float(x)
    except: return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--persample_dir',required=True)
    ap.add_argument('--relate',required=True)
    ap.add_argument('--pop')
    ap.add_argument('--inbreed')
    ap.add_argument('--ibs')
    ap.add_argument('--bamlist',required=True)
    ap.add_argument('--out',required=True)
    ap.add_argument('--flags',required=True)
    a=ap.parse_args()

    # canonical order + names from bam list (== ngsRelate index order)
    names=[os.path.basename(l.strip())[:-4] if l.strip().endswith('.bam')
           else os.path.basename(l.strip()) for l in open(a.bamlist) if l.strip()]
    idx2name={i:n for i,n in enumerate(names)}
    N=len(names)

    # per-sample rows
    rows={}
    for fn in os.listdir(a.persample_dir):
        if not fn.endswith('.qc.tsv'): continue
        with open(os.path.join(a.persample_dir,fn)) as fh:
            r=list(csv.DictReader(fh,delimiter='\t'))
            if r: rows[r[0]['sample']]=r[0]

    # relate.res: per-individual relatedness summaries
    relR0={i:[] for i in range(N)}; relKING={i:[] for i in range(N)}
    relTHETA={i:[] for i in range(N)}; n_related={i:0 for i in range(N)}
    with open(a.relate) as fh:
        # sniff delimiter
        first=fh.readline(); delim='\t' if '\t' in first else ','
        fh.seek(0)
        rd=csv.reader(fh,delimiter=delim)
        hdr=[h.strip() for h in next(rd)]
        def col(*cands):
            for c in cands:
                if c in hdr: return hdr.index(c)
            return None
        ca=col('a'); cb=col('b');
        if ca is None: ca=0
        if cb is None: cb=1
        cR0=col('R0'); cK=col('KING'); cT=col('theta')
        for parts in rd:
            if len(parts)<=max(filter(None,[ca,cb,cR0,cK,cT])): continue
            try:
                ai=int(parts[ca]); bi=int(parts[cb])
                R0=f(parts[cR0]); K=f(parts[cK]); T=f(parts[cT])
            except: continue
            for i in (ai,bi):
                if R0 is not None: relR0[i].append(R0)
                if K  is not None: relKING[i].append(K)
                if T  is not None: relTHETA[i].append(T)
            if K is not None and T is not None and K>=KING_2ND and T<THETA_UNREL:
                n_related[ai]+=1; n_related[bi]+=1

    medR0 ={i:med(relR0[i])  for i in range(N)}
    medK  ={i:med(relKING[i])for i in range(N)}
    medT  ={i:med(relTHETA[i])for i in range(N)}

    # inbreeding F (bamlist order)
    F={}
    if a.inbreed and os.path.exists(a.inbreed):
        vals=[f(x) for x in open(a.inbreed).read().split()]
        for i,v in enumerate(vals[:N]): F[i]=v

    # cohort IBS -> duplicate detection
    dup_partner={i:'' for i in range(N)}
    ref=None
    if a.ibs and os.path.exists(a.ibs):
        M=[[f(x) for x in line.split()] for line in open(a.ibs) if line.strip()]
        off=[M[i][j] for i in range(len(M)) for j in range(len(M)) if i!=j and M[i][j] is not None]
        ref=med(off)
        if ref:
            for i in range(min(N,len(M))):
                for j in range(min(N,len(M))):
                    if i<j and M[i][j] is not None and M[i][j] < DUP_FRAC*ref:
                        dup_partner[i]+=(idx2name.get(j,str(j))+';')
                        dup_partner[j]+=(idx2name.get(i,str(i))+';')

    # cohort distributions for robust outlier flags
    def colvals(key):
        return [f(rows[n][key]) for n in names if n in rows and key in rows[n]]
    dist={k:(med(colvals(k)), scale(colvals(k))) for k in
          ['mean_depth','breadth_1x','breadth_4x','pct_dup','error_rate',
           'ins_size','read_len','heterozygosity']}
    Fm=med([F[i] for i in F]); Fsc=scale([F[i] for i in F])
    R0m=med([medR0[i] for i in range(N) if medR0[i] is not None])
    R0sc=scale([medR0[i] for i in range(N) if medR0[i] is not None])

    # within-bam mixture reference: a "different individual" IBS level
    # (from cohort off-diagonal median computed above if available)
    mix_ref = ref if ref else None

    # assemble + decide
    out_cols=['sample','index','VERDICT','FLAGS',
              'mean_depth','breadth_1x','breadth_4x','pct_dup','error_rate',
              'ins_size','read_len','pct_mapped','heterozygosity','inbreed_F',
              'med_R0','med_KING','med_theta','n_related_partners',
              'n_RG','n_units','mix_ibs_max','dup_partner','distinct_SM']
    master=[]; flagged=[]
    for i,n in enumerate(names):
        r=rows.get(n,{})
        flags=[]
        # mixture / contamination (FAIL-level)
        # 1 related-to-many (the core signature)
        if n_related[i] >= RELATED_MANY_N:
            flags.append(f"RELATED_TO_MANY(n={n_related[i]})")
        # 2 low R0 outlier
        if rz(medR0[i],R0m,R0sc) < -Z:
            flags.append(f"LOW_R0(z={rz(medR0[i],R0m,R0sc):.1f})")
        # 3 excess heterozygosity
        het=f(r.get('heterozygosity')); m,sc=dist['heterozygosity']
        if rz(het,m,sc) > Z: flags.append(f"EXCESS_HET(z={rz(het,m,sc):.1f})")
        # 4 strongly negative inbreeding F
        Fi=F.get(i)
        if Fi is not None and (Fi < -0.10 or rz(Fi,Fm,Fsc) < -Z):
            flags.append(f"NEG_INBREED_F(F={Fi:.3f})")
        # 5 within-bam unit mixture
        mmax=f(r.get('mix_ibs_max'))
        if mix_ref and mmax is not None and mmax > 0.5*mix_ref:
            flags.append(f"MIXED_RG_UNITS(ibs={mmax:.3f})")
        # 6 duplicate / swap
        if dup_partner[i]:
            flags.append(f"DUPLICATE_OF({dup_partner[i].strip(';')})")
        # 7 SM inconsistency
        if r.get('n_distinct_SM','1') not in ('1','',None):
            flags.append(f"MULTI_SM({r.get('distinct_SM')})")

        contam = any(x.split('(')[0] in
                     ('RELATED_TO_MANY','LOW_R0','EXCESS_HET','NEG_INBREED_F',
                      'MIXED_RG_UNITS','DUPLICATE_OF','MULTI_SM') for x in flags)

        # technical-quality outliers (REVIEW-level)
        for key in ['mean_depth','breadth_1x','breadth_4x']:
            m,sc=dist[key]; z=rz(f(r.get(key)),m,sc)
            if z < -Z: flags.append(f"LOW_{key.upper()}(z={z:.1f})")
        for key in ['pct_dup','error_rate']:
            m,sc=dist[key]; z=rz(f(r.get(key)),m,sc)
            if z >  Z: flags.append(f"HIGH_{key.upper()}(z={z:.1f})")
        for key in ['ins_size','read_len']:
            m,sc=dist[key]; z=rz(f(r.get(key)),m,sc)
            if abs(z) > Z: flags.append(f"ODD_{key.upper()}(z={z:.1f})")

        quality = any(x.startswith(('LOW_','HIGH_','ODD_')) for x in flags)
        verdict = 'FAIL' if contam else ('REVIEW' if quality else 'PASS')

        master.append([n,i,verdict,';'.join(flags) or '-',
            r.get('mean_depth'),r.get('breadth_1x'),r.get('breadth_4x'),
            r.get('pct_dup'),r.get('error_rate'),r.get('ins_size'),r.get('read_len'),
            r.get('pct_mapped'),r.get('heterozygosity'),
            f"{Fi:.4f}" if Fi is not None else 'NA',
            f"{medR0[i]:.4f}" if medR0[i] is not None else 'NA',
            f"{medK[i]:.4f}" if medK[i] is not None else 'NA',
            f"{medT[i]:.5f}" if medT[i] is not None else 'NA',
            n_related[i],r.get('n_RG'),r.get('n_units'),r.get('mix_ibs_max'),
            dup_partner[i].strip(';'),r.get('distinct_SM')])
        if verdict!='PASS': flagged.append((n,verdict,';'.join(flags)))

    with open(a.out,'w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t'); w.writerow(out_cols); w.writerows(master)

    with open(a.flags,'w') as fh:
        npass=sum(1 for m in master if m[2]=='PASS')
        fh.write(f"QC SUMMARY: {len(master)} samples | PASS={npass} "
                 f"REVIEW={sum(1 for m in master if m[2]=='REVIEW')} "
                 f"FAIL={sum(1 for m in master if m[2]=='FAIL')}\n\n")
        if not flagged: fh.write("No flagged samples - dataset is clean.\n")
        for n,v,fl in sorted(flagged,key=lambda x:x[1]):
            fh.write(f"[{v}] {n}\n        {fl}\n")
    print(open(a.flags).read())

if __name__=='__main__':
    main()

