import pandas as pd

def import_gff(gff):
    names = ['contig', 'method','type','start','stop','dot','strand','value','desc']
    bdf = []
    
    with open(gff) as ih:
        for l in ih:
            if l.startswith("#"):
                continue
            else:                
                bdf.append(dict(zip(names, l.strip().split("\t"))))
                    
    return pd.DataFrame(bdf)