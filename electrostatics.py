'''
Electrostatics (simple alternative to the APBS Tools Plugin)

(c) 2012 Thomas Holder

License: BSD-2-Clause
'''

from __future__ import print_function

from pymol import cmd, CmdException

template_apbs_in = '''
read
    mol pqr "{pqrfile}"
end
elec
    mg-auto
    mol 1

    fgcent {fgcent} # fine grid center
    cgcent mol 1    # coarse grid center
    fglen {fglen}
    cglen {cglen}
    dime {dime}
    lpbe          # l=linear, n=non-linear Poisson-Boltzmann equation
    bcfl sdh      # "Single Debye-Hueckel" boundary condition
    pdie 2.0      # protein dielectric
    sdie 78.0     # solvent dielectric
    chgm spl2     # Cubic B-spline discretization of point charges on grid
    srfm smol     # smoothed surface for dielectric and ion-accessibility coefficients
    swin 0.3
    temp 310.0    # temperature
    sdens 10.0
    calcenergy no
    calcforce no
    srad {srad}   # solvent radius

    ion charge +1 conc 0.15 radius 2.0
    ion charge -1 conc 0.15 radius 1.8

    write pot dx "{mapfile}"
end
quit
'''

def validate_apbs_exe(exe):
    '''Get and validate apbs executable.
    Raise CmdException if not found or broken.'''
    import os, subprocess

    if exe:
        exe = cmd.exp_path(exe)
    else:
        try:
            import freemol.apbs
            exe = freemol.apbs.get_exe_path()
        except:
            pass
        if not exe:
            exe = cmd.exp_path('$SCHRODINGER/utilities/apbs')
            if not os.path.exists(exe):
                exe = "apbs"

    try:
        r = subprocess.call([exe, "--version"],
                stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        if r < 0:
            raise CmdException("Broken executable: " + exe)
    except OSError:
        raise CmdException("Cannot execute: " + exe)

    return exe

def map_new_apbs(name, selection='all', grid=0.5, buffer=10.0, state=1,
        preserve=0, exe='', assign=-1, focus='', quiet=1, _template=''):
    '''
DESCRIPTION

    Create electrostatic potential map with APBS.

    For more control over parameters and a graphical user interface I
    recommend to use the APBS Tools Plugin instead.

    In case of missing atoms or residues I recommend to remodel the input
    structure with modeller before calculating the electrostatic potential.

    If selection has no charges and radii, they will be automatically assigned
    with PyMOL (not with pdb2pqr).

SEE ALSO

    apbs_surface, map_new (coulomb), APBS Tools Plugin
    '''
    import tempfile, os, shutil, glob, subprocess
    from pymol.util import protein_assign_charges_and_radii
    from .modelling import add_missing_atoms

    selection = '(%s) and not solvent' % (selection)
    grid, buffer, state = float(grid), float(buffer), int(state)
    preserve, assign, quiet = int(preserve), int(assign), int(quiet)

    # path to apbs executable
    exe = validate_apbs_exe(exe)

    # temporary directory
    tempdir = tempfile.mkdtemp()
    if not quiet:
        print(' Tempdir:', tempdir)

    # filenames
    pqrfile = os.path.join(tempdir, 'mol.pqr')
    infile = os.path.join(tempdir, 'apbs.in')
    stem = os.path.join(tempdir, 'map')

    # temporary object
    tmpname = cmd.get_unused_name('mol' if preserve else '_')
    cmd.create(tmpname, selection, state, 1)

    # partial charges
    assign = [assign]
    if assign[0] == -1:
        # auto detect if selection has charges and radii
        cmd.iterate(tmpname,
                'assign[0] *= (elec_radius * partial_charge) == 0.0',
                space=locals())
    if assign[0]:
        cmd.remove('hydro and model ' + tmpname)
        add_missing_atoms(tmpname, quiet=quiet)
        protein_assign_charges_and_radii(tmpname)
    elif not quiet:
        print(' Notice: using exsiting charges and radii')

    cmd.save(pqrfile, tmpname, 1, format='pqr', quiet=quiet)

    # grid dimensions
    extent = cmd.get_extent(tmpname)
    extentfocus = cmd.get_extent(focus) if focus else extent
    fglen = [(emax-emin + 2*buffer) for (emin, emax) in zip(*extentfocus)]
    cglen = [(emax-emin + 4*buffer) for (emin, emax) in zip(*extent)]

    if not preserve:
        cmd.delete(tmpname)

    apbs_in = {
        'pqrfile': pqrfile,
        'fgcent': 'mol 1',
        'fglen': '%f %f %f' % tuple(fglen),
        'cglen': '%f %f %f' % tuple(cglen),
        'srad': cmd.get('solvent_radius'),
        'mapfile': stem,
    }

    if focus:
        apbs_in['fgcent'] = '%f %f %f' % tuple((emax + emin) / 2.0
                for (emin, emax) in zip(*extentfocus))

    try:
        # apbs will fail if grid does not fit into memory
        # -> on failure repeat with larger grid spacing
        for _ in range(3):
            dime = [1 + max(64, n / grid) for n in fglen]
            apbs_in['dime'] = '%d %d %d' % tuple(dime)

            # apbs input file
            with open(infile, 'w') as f:
                f.write((_template or template_apbs_in).format(**apbs_in))

            # run apbs
            r = subprocess.call([exe, infile], cwd=tempdir)
            if r == 0:
                break

            if r in (-6, -9):
                grid *= 2.0
                if not quiet:
                    print(' Warning: retry with grid =', grid)
                continue

            raise CmdException('apbs failed with code ' + str(r))

        dx_list = glob.glob(stem + '*.dx')
        if len(dx_list) != 1:
            raise CmdException('dx file missing')

        # load map
        cmd.load(dx_list[0], name, quiet=quiet)
    except OSError:
        raise CmdException('Cannot execute "%s"' % (exe))
    finally:
        if not preserve:
            shutil.rmtree(tempdir)
        elif not quiet:
            print(' Notice: not deleting %s' % (tempdir))

def apbs_surface(selection='all', maximum=None, minimum=None, map_name=None,
        ramp_name=None, grid=0.5, quiet=1):
    '''
DESCRIPTION

    Show electrostatic potential on surface (calculated with APBS).

    Important: surface_color is a object property, so when calculating
    surface potential for different selections and visualize them all
    together, you should first split them into separate objects.

USAGE

    apbs_surface [ selection [, maximum [, minimum ]]]

EXAMPLE

    fetch 2x19, bsync=0
    split_chains
    apbs_surface 2x19_A, 10
    apbs_surface 2x19_B, 10

SEE ALSO

    map_new_apbs, APBS Tools Plugin, isosurface, gradient,
    util.protein_vacuum_esp
    '''
    quiet = int(quiet)

    if ramp_name is None:
        ramp_name = cmd.get_unused_name('ramp')
    if map_name is None:
        map_name = cmd.get_unused_name('map')

    map_new_apbs(map_name, selection, float(grid), quiet=quiet)

    if maximum is not None:
        maximum = float(maximum)
        minimum = -maximum if minimum is None else float(minimum)
        kwargs = {'range': [minimum, (minimum+maximum)*0.5, maximum]}
    else:
        kwargs = {'selection': selection}

    cmd.ramp_new(ramp_name, map_name, **kwargs)

    object_names = cmd.get_object_list('(' + selection + ')')
    for name in object_names:
        cmd.set('surface_color', ramp_name, name)

    cmd.show('surface', selection)
    cmd.set('surface_solvent', 0)
    cmd.set('surface_ramp_above_mode', 1)

def volume_esp(name, map, stops=[0.1, 1.0], neg='red', pos='blue',
        opacity=0.2, quiet=1):
    '''
DESCRIPTION

    Create a volume object from a map object with default coloring
    for electrostatic potential (similar to positive and negative
    isosurface).

ARGUMENTS

    name = string: name for the new volume object

    map = string: name of the map object to use

    stops = list of floats: 2 or 3 values in standard deviations for creating
    the volume ramp {default: [0.1, 1.0]}

    neg = string: color for negative volume {default: red}

    pos = string: color for positive volume {default: blue}

    opacity = float: maximum opacity in volume ramp {default: 0.2}

SEE ALSO

    volume
    '''
    from .setting import set_temporary

    opacity, quiet = float(opacity), int(quiet)

    if isinstance(stops, str):
        stops = cmd.safe_list_eval(stops)

    try:
        from pymol.colorramping import ColorRamp
    except ImportError:
        print(' Warning: volume_esp is deprecated')
        stdevD = cmd.get_volume_histogram(map, 0)[3]
        stops = [s * stdevD for s in stops]
        ramp = [
            -stops[1], neg, opacity,
            -stops[0], neg, 0.0,
            stops[0], pos, 0.0,
            stops[1], pos, opacity,
        ]
        if len(stops) == 3:
            ramp = [-stops[2], neg, opacity] + ramp + [stops[2], pos, opacity]
        cmd.volume(name, map, ramp, quiet=quiet)
        return

    c_neg = cmd.get_color_tuple(neg)
    c_pos = cmd.get_color_tuple(pos)

    c_pos_0 = c_pos + (0.0,)
    c_pos_1 = c_pos + (opacity,)
    c_neg_0 = c_neg + (0.0,)
    c_neg_1 = c_neg + (opacity,)

    if len(stops) == 2:
        cstops = [(c_neg_1, -999), (c_neg_1, -stops[1]), (c_neg_0, -stops[0]),
                (c_pos_0, stops[0]), (c_pos_1, stops[1]), (c_pos_1, 999)]
    elif len(stops) == 3:
        cstops = [(c_neg_0, -999), (c_neg_0, -stops[2]), (c_neg_1, -stops[1]),
                (c_neg_0, -stops[0]), (c_pos_0, stops[0]), (c_pos_1, stops[1]),
                (c_pos_0, stops[2]), (c_pos_0, 999)]
    else:
        print(' Error: need 2 or 3 stops')
        raise CmdException

    cmd.volume(name, map, quiet=quiet)

    # get_volume_histogram returns zeros without refresh
    with set_temporary(suspend_updates='off'):
        cmd.refresh()

    minD, maxD, meanD, stdevD = cmd.get_volume_histogram(name)[:4]

    v_ramp = []
    c_ramp = ColorRamp(360)
    for c, s in cstops:
        i = int(360 * ((s * stdevD) - minD) / (maxD - minD))
        i = min(max(i, 0), 359)
        v_ramp.append(i)
        v_ramp.extend(c)
        c_ramp.addColor(i, c)

    cmd.set_volume_ramp(name, v_ramp)
    cmd.volume_color(name, c_ramp.getRamp())
    cmd.recolor(name)

def volume_fofc(name, map, stops=[2.5, 3.0, 4.0], neg='red', pos='green',
        opacity=0.4, quiet=1):
    '''
DESCRIPTION

    Create a difference density volume object.

SEE ALSO

    volume_esp
    '''
    return volume_esp(**locals())

cmd.extend('map_new_apbs', map_new_apbs)
cmd.extend('apbs_surface', apbs_surface)
cmd.extend('volume_esp', volume_esp)
cmd.extend('volume_fofc', volume_fofc)

cmd.auto_arg[0].update([
    ('apbs_surface', cmd.auto_arg[0]['zoom']),
])
cmd.auto_arg[1].update([
    ('map_new_apbs', cmd.auto_arg[0]['zoom']),
    ('volume_esp', cmd.auto_arg[1]['volume']),
    ('volume_fofc', cmd.auto_arg[1]['volume']),
])

# vi:expandtab:smarttab
