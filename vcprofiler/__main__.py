import os.path as op
from scgc import safe_makedir
from viruscope import run_prodigal, readfa
import yaml
import pkg_resources
from collections import defaultdict
import click
import pandas as pd
import sys

from .blast_tools import (blastp, get_best_hits, import_blastp, summarise_vog_blast,
hits_per_contig, vog_lca_per_contig)
from .hmmer_tools import vog_hmm


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.0.0')
@click.pass_context
def cli(obj):
    """Virus contig profiling"""
    pass

@cli.command('og-blast', short_help='examine phage orfs for membership in IMG/VR and VOG OGs')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--prot', is_flag=True)
@click.option('--yaml-file', "/mnt/scgc/scgc_nfs/lab/julia/notebooks/VCProfiler/vcprofiler.yaml", show_default=True)
@click.option('--outdir', type=click.Path(), default='./')
@click.option('--threads', default=1)
@click.option('--min-align-pct', default=95)
@click.option('--min-pct-id', default=35)
@click.option('--min-bitscore', default=50)
def vir_og_comparison(fasta, prot, yaml_file=None, outdir="./", threads=1,
                min_align_pct=95, min_pct_id=35,min_bitscore=50):
    outdir = safe_makedir(outdir)

    if prot is False:
        print("running prodigal", file=sys.stderr)
        proddir = safe_makedir(op.join(outdir, "prodigal"))
        prots = run_prodigal(fasta, proddir)
    else:
        prots = fasta

    if yaml_file is None:
        resource_package = __name__
        resource_path = 'vcprofiler.yaml'  
        yaml_file = pkg_resources.resource_string(resource_package, resource_path)

    config = yaml.load(open(yaml_file))
    name = op.basename(fasta).split(".")[0]

    contigsummary = op.join(outdir, "{}_contig_summary.csv".format(name))

    # first blast against VOGs from the University of Vienna
    print("running blast comparison against VOGs", file=sys.stderr)
    outblast = op.join(outdir, "{}_vs_vog_blastp.out".format(name))
    outsummary = op.join(outdir, "{}_vog_orf_summary.csv".format(name))


    blastdb = config['databases']['vog']['blast']
    lca = config['databases']['vog']['taxonomy']

    vogblast = blastp(prots, outblast, blastdb, threads=threads)
    vbdf = get_best_hits(import_blastp(vogblast), min_align_pct=min_align_pct,
                        min_pct_id=min_pct_id, min_bitscore=min_bitscore)
    summarise_vog_blast(vbdf).to_csv(outsummary, index=False)

    vog_contig_summary = pd.merge(hits_per_contig(vbdf, "orfs_matching_vogs"),
            vog_lca_per_contig(vbdf, lca), on='contig', how='outer')

    #next blast against img/vr
    print("comparing orfs to img/vr database", file=sys.stderr)
    outblast = op.join(outdir, "{}_vs_imgvr_blastp.out".format(name))
    blastdb = config['databases']['imgvr']['blast']['prot']
    out_contigsummary = op.join(outdir, '{}_vs_imgvr_contg_summary.csv')

    iblast = blastp(prots, outblast, blastdb, threads=threads)
    ibdf = get_best_hits(import_blastp(iblast), min_align_pct=min_align_pct,
                        min_pct_id=min_pct_id, min_bitscore=min_bitscore)
    ichit = hits_per_contig(ibdf, "imgvr_orf_hits")

    # summarise orfs per contig
    print("summarising findings", file=sys.stderr)
    contig_orf_count = defaultdict(lambda:0)
    for name, seq in readfa(open(prot)):
        contig_orf_count["_".join(name.split()[0].split("_")[:-1])] += 1

    orf_ct_df = pd.DataFrame(list(contig_orf_count.items()), columns=['contig','total_orfs'])

    df1 = pd.merge(vog_contig_summary, ichit, on='contig', how='outer')
    pd.merge(df1, orf_ct_df, on='contig', how='outer').to_csv(contigsummary)
    

@cli.command('vog-hmm', short_help='examine phage orfs for membership in IMG/VR and VOG OGs')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outdir', default="./vog_hmm", show_default=True)
@click.option('--yaml-file', default="/mnt/scgc/scgc_nfs/lab/julia/notebooks/VCProfiler/vcprofiler.yaml", show_default=True)
@click.option('--prot', is_flag=True, help="use if your input file protein sequence")
@click.option('--threads', default=1, show_default=True)
@click.option('--max-eval', default=0.001, show_default=True)
@click.option('--overwrite', is_flag=True, help="overwrite existing hmm search")
def run_vog_hmm(fasta, outdir, yaml_file, prot, threads=1, max_eval = 0.001, overwrite=False):
    
    config = yaml.load(open(yaml_file))
    db = config['databases']['vog']['hmm']
    vog_annotations = config['databases']['vog']['function']
    vog_lca = config['databases']['vog']['taxonomy']

    vog_hmm(fasta, outdir, prot=prot, threads=1, max_eval = 0.001, 
            db=db, 
            vog_annotations = vog_annotations,
            vog_lca = vog_lca, overwrite=overwrite)
              
              
              

if __name__=='__main__':
    cli()
