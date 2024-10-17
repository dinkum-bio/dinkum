#! /usr/bin/env python
import sys
import csv

import dinkum
from dinkum.vfg import Gene, Receptor, Ligand
from dinkum.vfn import Tissue
from dinkum import Timecourse

dinkum.reset()

# gataE, gcm, foxA

# set it all up!                                                            
pmar1 = Gene(name='pmar1')
hesC = Gene(name='hesC')
alx1 = Gene(name='alx1')
delta = Ligand(name='delta')
notch = Receptor(name='su(h)', ligand=delta)
tbr = Gene(name='tbr')
ets1 = Gene(name='ets1')
gataE = Gene(name='gataE')
gcm = Gene(name='gcm')
foxA = Gene(name='foxA')

embryo = Tissue(name='rest of embryo')
micromere = Tissue(name='micromeres')
embryo.add_neighbor(neighbor=micromere)

# maternal genes
early_ubiq = Gene(name='ub1')
late_ubiq = Gene(name='ub2')

## set up maternal gene expression

# early ubiq
early_ubiq.is_present(where=micromere, start=1)
early_ubiq.is_present(where=embryo, start=1)

# late ubiq
late_ubiq.is_present(where=micromere, start=4)
late_ubiq.is_present(where=embryo, start=4)

# pmar1: maternal in micromeres only
pmar1.is_present(where=micromere, start=1)

# hesC: present in both at the beginning
hesC.is_present(where=micromere, start=1, duration=1)
hesC.is_present(where=embryo, start=1, duration=1)

# notch: present everywhere (?)
notch.is_present(where=micromere, start=1)
notch.is_present(where=embryo, start=1)

## set up all downstream genes

# hesC: early, if not for pmar1
hesC.and_not(activator=early_ubiq, repressor=pmar1)

# alx &c: late, if not for hesC
alx1.and_not(activator=late_ubiq, repressor=hesC)
delta.and_not(activator=late_ubiq, repressor=hesC)
tbr.and_not(activator=late_ubiq, repressor=hesC)
ets1.and_not(activator=late_ubiq, repressor=hesC)

# gataE, gcm, foxA: all activated by notch
gataE.activated_by_and(sources=[late_ubiq, notch])
gcm.activated_by_and(sources=[late_ubiq, notch])
foxA.activated_by_and(sources=[late_ubiq, notch])

###

from dinkum import vfg, vfn

print(vfn._tissues)
print(vfg._rules)

def fill_ten(ls):
    return list(ls) + ['']*(10 - len(ls))

fp = open('out.btp.csv', 'w', newline='')
w = csv.writer(fp)

w.writerow(fill_ten(['# Model Commands']))
w.writerow(fill_ten(['# Command Type', 'Model Name', 'Parent Model']))
w.writerow(fill_ten(["model", "root"]))

base_name = "double neg"
w.writerow(fill_ten(["model", base_name, "root"]))

for tissue in vfn._tissues:
    w.writerow(fill_ten(["model", tissue.name, base_name]))
           
w.writerow(fill_ten(['# Region Commands']))
w.writerow(fill_ten(['# Command Type', 'Model Name', 'Region Name', 'Region Abbreviation']))

for n, tissue in enumerate(vfn._tissues):
    w.writerow(fill_ten(["region", tissue.name, f"a{n}", f"a{n}"]))

w.writerow(fill_ten(['# Standard Interactions']))
w.writerow(fill_ten(["# Command Type","Model Name","Source Type","Source Name","Target Type","Target Name","Sign","Source Region Abbrev","Target Region Abbrev"]))

for ix in vfg._rules:
    model_name = base_name
    source_type = 'gene'
    source_region = 'a0'
    target_region = 'a0'
    for (target, source, sign) in ix.btp_links():
        source_name = source.name
        target_name = target.name
        assert sign in ["positive", "negative"]

        w.writerow(["general",model_name, source_type, source_name, 'gene',
                    target_name, sign, source_region, target_region])

w.writerow(fill_ten(["# Signals"]))
w.writerow(fill_ten(["# Standalone nodes"]))
w.writerow(fill_ten(["# Interactions for Submodels"]))
w.writerow(["# Command Type","Model Name","Source Type","Source Name","Target Type","Target Name","Sign","Source Region Abbrev","Target Region Abbrev",])
fp.close()
