import os.path as op
from collections import defaultdict
import sys
import pandas as pd
from scgc import run, safe_makedir
from viruscope import run_prodigal, readfa
from scgc.utils import file_transaction


def run_hmmscan(infasta, db, outtbl, threads=1, overwrite=False):
    '''
    Args: 
    infasta: input fasta sequences as aa
    db: hmm database to compare against
    threads: how many threads to use
    overwrite: overwrite output file if it already exists
    
    Returns:
    path to the outtbl (string)
    '''
    if not overwrite and op.exists(outtbl):
        return outtbl
    with file_transaction(outtbl) as tx_out_file:
        cmd = "hmmscan --tblout {outtbl} --cpu {threads} {db} {infasta}".format(outtbl=tx_out_file, db=db,                                                                             infasta=infasta, threads=threads)
        run(cmd)
    return outtbl

def import_hmmscan_out(tbl):
    '''
    args: 
        tbl: hmmscan output
    returns:
        df: pandas dataframe of hmmscan table
    '''
    names = ['sseqid', 'accession', 'qseqid', 'accession', 'evalue', 'score', 'bias', 'd-evalue', 'd-score', 'd-bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description of target']
    df = pd.DataFrame(columns = names)
    with open(tbl) as ih:
        for i, l in enumerate(ih):
            if not l.startswith('#'):
                df.loc[i] = pd.Series(dict(zip(names, l.split())))
    return df

def hmmscan_besthits(df, max_eval=0.001):
    '''
    Keeps only best hits from a pandas dataframe with the columns 'qseqid' and 'evalue'
    Args:
        df: pandas dataframe of hmmscan output (created using the function import_hmmscan_out), it really only needs two columns with the headers 'qseqid' and 'evalue' to work
    Returns:
        pandas dataframe with only the hit with the lowest evalue per qseqid present
    '''
    df['evalue'] = [float(i) for i in df['evalue']]
    df = df[df['evalue'] < max_eval]
    return df.sort_values(by=['qseqid','evalue'], ascending=True).drop_duplicates(subset='qseqid', keep='first')


def process_hmmscan(tbl, max_eval=0.001, best_only=True):
    '''
    Processes the hmmscan output table
    Args:
        tbl (path): hmmscan output
        max_eval (numeric): maximum evalue to maintain in table
        best_only (bool): keep only the best hits if True
    Returns:
        pandas dataframe of output
    '''
    names = ['sseqid', 'accession', 'qseqid', 'accession', 'evalue', 'score', 'bias', 'd-evalue', 'd-score', 'd-bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description of target']
    df = pd.DataFrame(columns = names)
    with open(tbl) as ih:
        for i, l in enumerate(ih):
            if not l.startswith('#'):
                toks = dict(zip(names, l.split()))
                if float(toks['evalue']) <= max_eval:
                    df.loc[i] = pd.Series(dict(zip(names, l.split())))
        if best_only:
            return df.sort_values(by=['qseqid','evalue'], ascending=True).drop_duplicates(subset='qseqid', keep='first')
        else:
            return df

def hmm_vog_summary(vhmm_tops, 
                vog_annotations = "/mnt/scgc/simon/simonsproject/jb_vs_test/viral_dbs/VOG/vog.annotations.tsv", 
                vog_lca = "/mnt/scgc/simon/simonsproject/jb_vs_test/viral_dbs/VOG/vog.lca.tsv"):
    vdf = vhmm_tops[['qseqid','sseqid']].rename(columns={'sseqid':'vog','qseqid':'orf'})
    voganns = pd.read_csv(vog_annotations, sep="\t")
    voganns.rename(columns={'#GroupName':'vog','ConsensusFunctionalDescription':'func_desc'}, inplace=True)
    voglca = pd.read_csv(vog_lca, skiprows=1, sep="\t", names=['vog','genomes_in_group_and_lca','genomes_total_in_lca','lca'])
    summary = pd.merge(vdf, 
                       voganns, 
                       on='vog', 
                       how='left').merge(voglca, 
                                         on='vog', 
                                         how='left')
    summary['contig'] = ["_".join(i.split("_")[:-1]) for i in summary['orf']]
    summary['orf_num'] = [int(i.split("_")[-1]) for i in summary['orf']]
    summary = summary.sort_values(by=['contig','orf_num']).reset_index()
    return summary[['orf','contig', 'vog','func_desc','lca']]

def hmm_hits_per_contig(summary, hits_category='vogs'):

    out = summary.groupby('contig')['orf'].count().reset_index()
    out = out.rename(columns={'orf':'orfs_matching_{}'.format(hits_category)})
    return out


def hmm_vog_lca_per_contig(summary):
    best_guesses = []

    clca_count = summary.groupby(['contig','lca'])['orf'].count().reset_index()
    for c in clca_count['contig'].unique():
        subdf = clca_count[clca_count['contig'] == c]
        maxcount = 0
        for i, l in subdf.iterrows():
            if l['orf'] > maxcount:
                maxcount = l['orf']
                best_guess = l
            elif l['orf'] == maxcount:
                if len(l['lca'].split(";")) > len(best_guess['lca'].split(";")):
                    best_guess = l
                else:
                    continue
        best_guesses.append(best_guess)
    return pd.DataFrame(best_guesses)[['contig','lca','orf']].rename(columns={'orf':'hit_count', 'lca':'vog_lca'})

def total_orfs(prot_fasta):
    contig_orf_count = defaultdict(lambda:0)
    for name, seq in readfa(open(prot_fasta)):
        contig_orf_count["_".join(name.split()[0].split("_")[:-1])] += 1

    orf_ct_df = pd.DataFrame(list(contig_orf_count.items()), columns=['contig','total_orfs'])
    return orf_ct_df


def vog_hmm(infasta, outdir, prot=False, threads=1, max_eval = 0.001, 
            db='/mnt/scgc/simon/simonsproject/jb_vs_test/viral_dbs/VOG/hmm/vog_hmm',
            vog_annotations = "/mnt/scgc/simon/simonsproject/jb_vs_test/viral_dbs/VOG/vog.annotations.tsv",
            vog_lca = "/mnt/scgc/simon/simonsproject/jb_vs_test/viral_dbs/VOG/vog.lca.tsv", overwrite=False):
    
    ''' uses VOGs and tables from: http://vogdb.org/download '''
    print(vog_annotations)
    print(db)
    
    outdir = safe_makedir(outdir)
    name = op.basename(infasta).split(".")[0]
    
    if prot is False:
        print("running prodigal", file=sys.stderr)
        proddir = safe_makedir(op.join(outdir, "prodigal"))
        prots = run_prodigal(infasta, proddir)
    else:
        prots = infasta
    
    out_hmm = op.join(outdir, '{}_vs_vog_hmmer.tbl'.format(name))
    orf_summary = op.join(outdir, '{}_vs_vog_orf_summary.csv'.format(name))
    contig_summary = op.join(outdir, '{}_vs_vog_contig_summary.csv'.format(name))
    
    out_hmm = run_hmmscan(prots, db, out_hmm, threads=threads, overwrite=overwrite)
    
    df = hmmscan_besthits(import_hmmscan_out(out_hmm), max_eval=max_eval)
    vsum = hmm_vog_summary(df, vog_annotations=vog_annotations, vog_lca=vog_lca)
    vsum.to_csv(orf_summary, index=False)
    
    chits = hmm_hits_per_contig(vsum)
    lcahits = hmm_vog_lca_per_contig(vsum)
    total_orfs_per_contig = total_orfs(prots)
    
    csum = pd.merge(chits, lcahits, on='contig', how='outer').merge(total_orfs_per_contig, on='contig', how='outer')
    csum.to_csv(contig_summary)
    return contig_summary, orf_summary

# plotting

## plotting
from matplotlib import gridspec

def plot_contig_results(cfile, title='',outdir=None):
    outdir = safe_makedir(outdir)
    
    def _grab_fam(i, position, debug=False):
        try:
            if len(i.split(";")) >= position + 1:
                return i.split(';')[position]
            else:

                return "_{}".format(i.split(';')[-1])
        except Exception as inst:
            if debug:
                print(i, inst)
            return 'Unknown'
    
    cdf = pd.read_csv(cfile)
    cdf['pct_vog_match'] = cdf['orfs_matching_vogs'] / cdf['total_orfs'] * 100
    
    print(round(len(cdf.dropna()) / len(cdf) * 100), "% had hits to the VOG db", sep="", file=sys.stdout)
    print(round(len([i for i in cdf['vog_lca'].dropna() if 'Caudovirales' in i]) / len(cdf) * 100), 
          "% classified as Caudovirales", sep="", file=sys.stdout)
    famcounts = pd.DataFrame(pd.Series(Counter([_grab_fam(i, 2) for i in cdf['vog_lca']]), 
                                       name='vir_fam_count')).reset_index().sort_values(by='vir_fam_count', ascending=False)
    sns.set_style("whitegrid")
    
    gridkw = dict(width_ratios=[3, 1])
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=gridkw, figsize=(15, 5))
    
    a = sns.barplot(x="index", y="vir_fam_count", data=famcounts, ax=ax1)
    a.set_xticklabels(a.get_xticklabels(), rotation=90)
    a.set_xlabel('Virus Group')
    a.set_ylabel('Contig Count')
    a.set_title('{} Contig VOG LCA'.format(title))
    
    caudo_cats = pd.DataFrame(pd.Series(Counter([_grab_fam(i, 3) for i in list(cdf.dropna()[cdf['vog_lca'].dropna().str.contains('Caudovirales')]['vog_lca'])]), 
                                        name='group_count')).reset_index()
    
    a1 = sns.barplot(x="index", y="group_count", data=caudo_cats, ax=ax2)
    a1.set_xticklabels(a1.get_xticklabels(), rotation=90)
    a1.set_xlabel('Virus Group')
    a1.set_ylabel('')
    a1.set_title('{} Contig VOG LCA\nCaudovirales Categories'.format(title))
    if outdir:
        fig.save_fig(op.join(outdir, 'vog_lca.png'))
        caudo_cats.to_csv(op.join(outdir, 'caudo_cats.csv'))
        famcounts.to_csv(op.join(outdir, 'family_table.csv', index=False))
    return fig, famcounts, caudo_cats
