'''
(c) 2018 satot & Yoshitaka Moriwaki, BILAB

License: BSD-2-Clause
'''

from __future__ import absolute_import
from __future__ import with_statement

from pymol import cmd, stored, CmdException

def comp(thre,val): return thre - 1e-12 < val

def contactmap(selection='', chain='', align='hmmrenumbered.map', pcut=0.999, offset=0):
    '''
DESCRIPTION

    visualize the results from map_align algorithm in PyMOL.
    https://github.com/sokrypton/map_align

ARGUMENTS

    selection = string: atom selection {default: all}

    chain = string: chain ID you want to map {default: ''}

    align = string: name of map_align result txt file {default: 'hmmrenumbered.map'}

    pcut = float: probability cutoff to display {default: 0.999}

    offset = int: slide the sequence number to fit its original number manually {default: 0}

EXAMPLE

    contactmap 5dir, A, /path/to/hmmrenumbered.map, 0.998, offset=-7

    '''
    if selection == '':
        print("ERROR: choose the model to display the contact map onto")
        raise CmdException
    model = selection

    stored.resdict = {}
    cmd.iterate(model, "stored.resdict[ int(resi) ] = resn")
    resnlist = stored.resdict.keys()

    with open(align) as f:
        l = f.readlines()[0:]
        for s in l:
            sec = s.split()

            a = int(sec[1]) + int(offset)
            b = int(sec[2]) + int(offset)
            probability = float(sec[4])

            # not display if probability < cutoff
            if not comp(float(pcut), probability) : continue
            if a * b == 0: continue
            if not a in resnlist or not b in resnlist: continue

            cara = carb = "CB"
            if stored.resdict[a] == "GLY":
                cara = "CA"
            if stored.resdict[b] == "GLY":
                carb = "CA"

            #cmd.distance(f"{a}_{b}_{model}",f"(/{model}//{chain}/{a}/{cara})",(f"(/{model}//{chain}/{b}/{carb})"))
            cmd.distance('{0}_{1}_{2}'.format(a, b, model), '/{0}//{1}/{2}/{3}'.format(model, chain, a, cara), '/{0}//{1}/{2}/{3}'.format(model, chain, b, carb))

    cmd.set('dash_color', 'cyan')
    cmd.set('dash_radius', '0.1')
    cmd.set('dash_gap', '0.05')
    print("===finished===")

cmd.extend('contactmap', contactmap)

# tab-completion of arguments
cmd.auto_arg[0].update({
    'contactmap'       : cmd.auto_arg[0]['zoom'],
})
