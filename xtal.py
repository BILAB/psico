'''
Crystallographic symmetry related commands.

(c) 2010-2012 Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def cellbasis(angles, edges):
    '''
DESCRIPTION

    API only. For the unit cell with given angles and edge lengths calculate
    the basis transformation (vectors) as a 4x4 numpy.array
    '''
    from math import cos, sin, radians, sqrt
    import numpy

    rad = [radians(i) for i in angles]
    basis = numpy.identity(4)
    basis[0][1] = cos(rad[2])
    basis[1][1] = sin(rad[2])
    basis[0][2] = cos(rad[1])
    basis[1][2] = (cos(rad[0]) - basis[0][1]*basis[0][2])/basis[1][1]
    basis[2][2] = sqrt(1 - basis[0][2]**2 - basis[1][2]**2)
    edges.append(1.0)
    return basis * edges # numpy.array multiplication!

def supercell(a=1, b=1, c=1, object=None, color='green', name='supercell', withmates=1):
    '''
DESCRIPTION

    Draw a supercell, as requested by Nicolas Bock on the pymol-users
    mailing list (Subject: [PyMOL] feature request: supercell construction
    Date: 04/12/2010 10:12:17 PM (Mon, 12 Apr 2010 14:12:17 -0600))

USAGE

    supercell a, b, c [, object [, color [, name [, withmates]]]]

ARGUMENTS

    a, b, c = integer: repeat cell in x,y,z direction a,b,c times
    {default: 1,1,1}

    object = string: name of object to take cell definition from

    color = string: color of cell {default: blue}

    name = string: name of the cgo object to create {default: supercell}

    withmates = bool: also create symmetry mates in displayed cells
    {default: 1}

SEE ALSO

    show cell

    '''
    import numpy
    from pymol import cgo

    if object is None:
        object = cmd.get_object_list()[0]
    withmates = int(withmates)

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]

    basis = cellbasis(cell_angles, cell_edges)
    assert isinstance(basis, numpy.ndarray)

    ts = list()
    for i in range(int(a)):
        for j in range(int(b)):
            for k in range(int(c)):
                ts.append([i,j,k])

    obj = [
        cgo.BEGIN,
        cgo.LINES,
    ]

    for t in ts:
        shift = basis[0:3,0:3] * t
        shift = shift[:,0] + shift[:,1] + shift[:,2]

        for i in range(3):
            vi = basis[0:3,i]
            vj = [
                numpy.array([0.,0.,0.]),
                basis[0:3,(i+1)%3],
                basis[0:3,(i+2)%3],
                basis[0:3,(i+1)%3] + basis[0:3,(i+2)%3]
            ]
            for j in range(4):
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j]).tolist())
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j] + vi).tolist())

        if withmates:
            groupname = 'm%d%d%d' % tuple(t)
            symexpcell(groupname + '_', object, *t)
            cmd.group(groupname, groupname + '_*')

    obj.append(cgo.END)

    cmd.delete(name)
    cmd.load_cgo(obj, name)
    cmd.color(color, name)

def symexpcell(prefix='mate', object=None, a=0, b=0, c=0):
    '''
DESCRIPTION

    Creates all symmetry-related objects for the specified object that
    occur with their bounding box center within the unit cell.

USAGE

    symexpcell prefix, object, [a, b, c]

ARGUMENTS

    prefix = string: prefix of new objects

    object = string: object for which to create symmetry mates

    a, b, c = integer: create neighboring cell {default: 0,0,0}

SEE ALSO

    symexp
    '''
    import numpy
    from pymol import xray

    if object is None:
        object = cmd.get_object_list()[0]

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]
    spacegroup = sym[6]

    basis = cellbasis(cell_angles, cell_edges)
    basis = numpy.matrix(basis)

    extent = cmd.get_extent(object)
    center = sum(numpy.array(extent)) * 0.5
    center = numpy.matrix(center.tolist() + [1.0]).T
    center_cell = basis.I * center

    extra_shift = [[float(i)] for i in (a,b,c)]

    spacegroup = xray.space_group_map.get(spacegroup, spacegroup)

    i = 0
    matrices = xray.sg_sym_to_mat_list(spacegroup)
    for mat in matrices:
        i += 1

        mat = numpy.matrix(mat)
        shift = numpy.floor(mat * center_cell)
        mat[0:3,3] -= shift[0:3,0]
        mat[0:3,3] += extra_shift

        mat = basis * mat * basis.I
        mat_list = list(mat.flat)

        name = '%s%d' % (prefix, i)
        cmd.create(name, object)
        cmd.transform_object(name, mat_list, 0)
        cmd.color(i+1, name)

def pdbremarks(filename):
    '''
DESCRIPTION

    API only. Read REMARK lines from PDB file. Return dictionary with
    remarkNum as key and list of lines as value.
    '''
    remarks = dict()
    if not cmd.is_string(filename):
        f = filename
    elif filename[-3:] == '.gz':
        import gzip
        f = gzip.open(filename)
    else:
        f = open(filename)
    for line in f:
        recname = line[0:6]
        if recname == 'REMARK':
            num = int(line[7:10])
            lstring = line[11:]
            remarks.setdefault(num, []).append(lstring)
    return remarks

def biomolecule(name=None, filename=None, prefix=None, number=1, suffix=None,
        quiet=0):
    '''
DESCRIPTION

    Create biological unit (quaternary structure) as annotated by the REMARK
    350 BIOMOLECULE record.

USAGE

    biomolecule name [, filename [, prefix [, number ]]]

ARGUMENTS

    name = string: name of object and basename of PDB file, if
    filename is not given {default: first loaded object}

    filename = string: file to read remarks from {default: <name>.pdb}

    prefix = string: prefix for new objects {default: <name>}

EXAMPLE

    fetch 1rmv, bsync=0
    biomolecule 1rmv
    '''
    import os
    from .importing import local_mirror_pdb

    try:
        import numpy
    except ImportError:
        numpy = None

    number, quiet = int(number), int(quiet)

    if name is None:
        name = cmd.get_object_list()[0]
    if prefix is None:
        prefix = name
    if suffix is None:
        suffix = str(number)
    if filename is None:
        candidates = [
            '%s.pdb' % (name),
            '%s/%s.pdb' % (cmd.get('fetch_path'), name),
            local_mirror_pdb(name),
        ]
        for filename in candidates:
            if os.path.exists(filename):
                break
        else:
            print('please provide filename')
            raise CmdException
        if not quiet:
            print('loading from %s' % (filename))

    remarks = pdbremarks(filename)
    if 350 not in remarks:
        print('There is no REMARK 350 in ' + filename)
        raise CmdException

    current = 1
    biomt = {current: {}}
    chains = tuple()

    for line in remarks[350]:
        if line.startswith('BIOMOLECULE:'):
            current = int(line[12:])
            biomt[current] = {}
        elif line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
            chains = tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('                   AND CHAINS:'):
            chains += tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('  BIOMT'):
            row = int(line[7])
            num = int(line[8:12])
            vec = line[12:].split()
            vec = list(map(float, vec))
            biomt[current].setdefault(chains, dict()).setdefault(num, []).extend(vec)

    if number not in biomt or len(biomt[number]) == 0:
        print(' Error: no BIOMOLECULE number %d' % (number))
        raise CmdException

    if numpy is not None:
        mat_source = numpy.reshape(cmd.get_object_matrix(name), (4,4))
        mat_source = numpy.matrix(mat_source)

    for chains, matrices in biomt[number].items():
        for num in matrices:
            mat = matrices[num][0:12]
            mat.extend([0,0,0,1])
            copy = '%s_%s_%d' % (prefix, suffix, num)
            if not quiet:
                print('creating %s' % (copy))
            cmd.create(copy, 'model %s and chain %s' % (name, '+'.join(chains)))
            cmd.alter(copy, 'segi="%d"' % (num))

            if numpy is not None:
                mat = mat_source * numpy.reshape(mat, (4,4)) * mat_source.I
                mat = list(mat.flat)

            cmd.transform_object(copy, mat)

    cmd.disable(name)
    cmd.group('%s_%s' % (prefix, suffix), '%s_%s_*' % (prefix, suffix))

cmd.extend('supercell', supercell)
cmd.extend('symexpcell', symexpcell)
cmd.extend('biomolecule', biomolecule)

# tab-completion of arguments
cmd.auto_arg[0]['biomolecule'] = cmd.auto_arg[0]['pseudoatom']
cmd.auto_arg[3]['supercell'] = cmd.auto_arg[0]['pseudoatom']

# vi: ts=4:sw=4:smarttab:expandtab
