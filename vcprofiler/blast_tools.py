import os
import os.path as op
import pandas as pd
import yaml
from scgc.utils import file_transaction, run, safe_makedir, file_exists
from viruscope import run_blast, readfa


def blastp(fasta, out_file, blast_db, threads=1, numalign=10, evalue=0.001, override=False):
    """
    align sequences using blastp, using -m 6 tab delimited output format.
    fasta : file path as string
    out_file : result file path with tsv extension
    options : blastp options except for db, num_threads, query, out, and outfmt
    returns : outfile path
    """
    fields = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
              'qend', 'sstart', 'send', 'evalue', 'bitscore', 'score', 'nident',
              'positive', 'gaps', 'ppos', 'qframe', 'sframe', 'qlen', 'slen']

    if file_exists(out_file) and not override:
        return out_file

    with file_transaction(out_file) as tx_out_file:
        cmd = run_blast(fasta, out_file, db=blast_db, num_alignments=numalign, evalue=evalue, threads=threads, fields=fields)

        run(cmd)

    return out_file

def import_blastp(blastpath, fields=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
      'qend', 'sstart', 'send', 'evalue', 'bitscore', 'score', 'nident', 'positive', 'gaps', 'ppos', 'qframe', 'sframe', 'qlen', 'slen']):
    return pd.read_csv(blastpath, names=fields, sep="\t")


def get_best_hits(blastdf, min_align_pct = 0, min_pct_id = 0, min_bitscore = 0):
    blastdf = blastdf[(blastdf['length'] / blastdf['qlen'] >= min_align_pct / 100) & (blastdf['pident'] >= min_pct_id) &  (blastdf['bitscore'] >= min_bitscore)]
    return blastdf.sort_values(by=['qseqid','bitscore'], ascending=False).drop_duplicates(subset=['qseqid'], keep='first')

def summarise_vog_blast(bdf):
    bdf['vog'] = [i.split("__")[-1].split(":")[0] for i in bdf['sseqid']]
    bdf['vog_desc'] = [i.split("__")[-1].split(":")[1].replace("_"," ") for i in bdf['sseqid']]
    bdf['contig'] = ["_".join(i.split("_")[:-1]) for i in bdf['qseqid']]
    outdf = bdf[['contig','qseqid','vog','vog_desc']].rename(columns={'qseqid':'orf'})
    outdf['orf_num'] = [int(i.split("_")[-1]) for i in outdf['orf']]
    outdf = outdf.sort_values(by=['contig','orf_num']).reset_index()
    return outdf[['contig','orf','vog','vog_desc']]

def hits_per_contig(bdf, hits_category):
    bdf['contig'] = ["_".join(i.split("_")[:-1]) for i in bdf['qseqid']]
    out = bdf.groupby('contig')['qseqid'].count().reset_index()
    out = out.rename(columns={'qseqid':'orfs_matching_{}'.format(hits_category)})
    return out

def vog_lca_per_contig(bdf_top, vog_lca_file):
    vog_lca = pd.read_csv(vog_lca_file, sep="\t")
    vog_lca = vog_lca.rename(columns=dict(zip(vog_lca.columns, ['vog','GenomesInGroupAndLCA', 'GenomesInLCA', 'LCA'])))

    best_guesses = []
    vbdf = summarise_vog_blast(bdf_top).merge(vog_lca[['vog','LCA']], on='vog', how='left')

    clca_count = vbdf.groupby(['contig','LCA'])['orf'].count().reset_index()
    for c in clca_count['contig'].unique():
        subdf = clca_count[clca_count['contig'] == c]
        maxcount = 0
        for i, l in subdf.iterrows():
            if l['orf'] > maxcount:
                maxcount = l['orf']
                best_guess = l
            elif l['orf'] == maxcount:
                if len(l['LCA'].split(";")) > len(best_guess['LCA'].split(";")):
                    best_guess = l
                else:
                    continue
        best_guesses.append(best_guess)
    return pd.DataFrame(best_guesses)[['contig', 'LCA','orf']].rename(columns={'orf':'hit_count', 'LCA':'vog_lca'})

'''
def alt_databases(fasta, prot=False, yaml_file=None, outdir="./", threads=1,
                min_align_pct=95, min_pct_id=35,min_bitscore=50):
    outdir = safe_makedir(outdir)

    if not prot:
        proddir = safe_makedir(op.join(outdir, "prodigal"))
        prots = run_prodigal(fasta, proddir)
    else:
        print("proteins provided.  Cannot run comparison against ssDNA viruses.", file=sys.stderr)
        prots = fasta
    
    name = op.basename(fasta).split(".")[0].split("_")[0]
    djrout = op.join(outdir, "{}_contigs_w_djr_proteins.csv".format(name))
    
    if not prot:
        print("comparing contigs to ssDNA virus database using blastn", file=sys.stderr)
        blastout = op.join(outdir, '{}_vs_ssdna_virs_blastn.out'.format(name))
        blastdb = config['databases']['blast']['ssdna_virs']
        blastout = blastn(fasta, blastout, blastdb, threads=threads, numalign=10, evalue=0.001)
        ssbdf = get_best_hits(import_blastp(blastout), min_align_pct=20, min_pct_id=95, min_bitscore=0)
        print("{} contigs hit ssDNA viruses".format(len(ssbdf), file=sys.stdout)
    
    print("comparing orfs to djr proteins", file=sys.stderr)
    blastdb = config['databases']['blast']['djr']
    out = op.join(outdir, "{}_vs_djr_proteins_blastp.out".format(name))
    blastout = blastp(prots, out, blastdb, threads=threads)
    djrbdf = get_best_hits(import_blsatp(blastout), min_align_pct=75, min_pct_id=35, min_bitscore=50)
    djr_chit = hits_per_contig(djrbdf)
    djrbdf.to_csv(djrout, index=False)
    print("{} contigs contain proteins resembling djr proteins.".format(len(djr_chit))
'''