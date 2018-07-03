import pandas as pd
import os.path as op
import networkx as nx
import random
import community
import sys

from scgc.utils import run

# running
def create_contig_msh(infasta, out_prefix=None, threads=1, k=21):
    if out_prefix is None:
        out_prefix = op.join(op.dirname, op.basename(infasta).split(".")[0])

    if op.exists('{}.msh'.format(out_prefix)):
        print("msh already exists.", file=sys.stderr)
        return '{}.msh'.format(out_prefix)

    cmd = 'mash sketch -p {threads} -i -k {k} -o {outprefix} {infasta}'.format(threads=threads, k=k, outprefix=out_prefix, infasta=infasta)
    print(cmd)
    run(cmd)
    return '{}.msh'.format(out_prefix)

def run_mash(fa1, fa2, threads=1, k=21, dist=0.1, outtbl='mashout.tsv'):
    if op.exists(outtbl):
        print("output table {} exists, aborting mash".format(outtbl), file=sys.stderr)
        return outtbl

    if not fa1.endswith('msh'):
        msh1 = create_contig_msh(fa1, threads=threads, k=k)
    else:
        msh1 = fa1

    if not fa2.endswith('msh'):
        msh2 = create_contig_msh(fa2, threads=threads, k=k)
    else:
        msh2 = fa2

    cmd = 'mash dist -d {dist} -p {threads} {msh1} {msh2} > {outtbl}'.format(dist=dist,
                                                                            threads=threads,
                                                                            msh1=msh1, msh2=msh2,
                                                                            outtbl=outtbl)
    run(cmd)
    return outtbl

# munging
def load_mash(mashin, keep_self=False, shorten_names=True):
    names = ['q1','q2','distance','pvalue','shared_hashes']
    df = pd.read_csv(mashin, names=names, sep="\t")
    if not keep_self:
        df = df[df['q1'] != df['q2']]
    if shorten_names:
        df['q1'] = [op.basename(i).split()[0] for i in df['q1']]
        df['q2'] = [op.basename(i).split()[0] for i in df['q2']]
    return df


def cluster_mash(mashin, distance = 0.05, keep_self=False, shorten_names=True,
                output=None):
    df = load_mash(mashin, keep_self=keep_self, shorten_names=shorten_names)
    subdf = df[df['distance'] <= distance]
    nodes = set(list(subdf['q1'].unique()) + list(subdf['q2'].unique()))
    g = nx.Graph()
    g.add_nodes_from(nodes)
    edges = [(l.q1, l.q2, {"distance":l.distance,"pvalue":l.pvalue}) for i, l in subdf.iterrows()]
    g.add_edges_from(edges)
    partition = community.best_partition(g)
    values = [partition.get(node) for node in g.nodes()]
    if output:
        with open(output, "w") as oh:
            for p in partition:
                print(p, partition[p], sep=',', file=oh)
    return pd.DataFrame.from_dict(partition, orient='index').reset_index().rename(columns={"index":'sample', 0:'d{}_group'.format(distance)})

def cluster_pairs(df, c1, c2):
    nodes = set(list(df[c1].unique()) + list(df[c2].unique()))
    g = nx.Graph()
    g.add_nodes_from(nodes)
    edges = [(l[c1], l[c2]) for i, l in df.iterrows()]
    g.add_edges_from(edges)
    partition = community.best_partition(g)
    values = [partition.get(node) for node in g.nodes()]
    return pd.DataFrame.from_dict(partition, orient='index').reset_index().rename(columns={"index":'sample', 0:'group'})
