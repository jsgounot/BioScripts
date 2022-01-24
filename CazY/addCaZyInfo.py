# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-12-16 12:23:50
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-12-16 13:31:50

import pandas as pd

diamond = snakemake.input[0]
names = ['Gene ID', 'CAZy ID', '%Identical', 'Length', 'Mismatches', 'Gap Open', 'Gene Start', 'Gene End', 'CAZy Start', 'CAZy End', 'E Value', 'Bit Score']
diamond = pd.read_csv(diamond, sep='\t', names=names)

diamond['CAZy GB'] = diamond['CAZy ID'].apply(lambda idv: idv.split('|')[0])
diamond['CAZy ID'] = diamond['CAZy ID'].apply(lambda idv: idv.split('|')[1])

def setcazygroup(cid):
    cid = ''.join(i for i in cid if i.isalpha())
    mapper = {
        'GH': 'Glycoside Hydrolase',
        'PL': 'Polysaccharide Lyase',
        'GT': 'Glycosyltransferases',
        'CBM': 'Carbohydrate-binding modules',
        'CE': 'Carbohydrate esterases',
        'AA': 'Auxiliary Activities'
    }
    return mapper[cid]

diamond['cazy_group'] = diamond['CAZy ID'].apply(setcazygroup)

cmeta = snakemake.params['cmeta']
names = ['CAZy ID', 'Domain', 'Species', 'CAZy GB']
cmeta = pd.read_csv(cmeta, sep='\t', names=names)

diamond = diamond.merge(cmeta, on='CAZy GB', how='left', suffixes=('_dia', '_meta'))

outfile = snakemake.output[0]
diamond.to_csv(outfile, sep='\t', index=False)