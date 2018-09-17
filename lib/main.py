"""
Main library for moldraw
"""

import svgwrite
import numpy as np
from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from lib.math import *

atom_dict = {
    1:  ['H',  1.00, 'rgb(182,225,255)', 'rgb(176,248,251)', 'black',  .01],
    5:  ['B',  1.30, 'rgb(255,158,158)', 'rgb(255,179,179)', 'black',  .02],
    6:  ['C',  1.23, 'rgb(128,128,128)', 'rgb(179,179,179)', 'black',  .00],
    7:  ['N',  1.54, 'rgb(128,000,128)', 'rgb(255,000,255)', 'white',  .00],
    8:  ['O',  1.77, 'rgb(212,000,000)', 'rgb(255,085,085)', 'black',  .00],
    9:  ['F',  1.80, 'rgb(098,146,000)', 'rgb(143,160,085)', 'black',  .05],
    14: ['Si', 1.85, 'rgb(215,254,220)', 'rgb(200,200,200)', 'black', -.05],
    15: ['P',  1.95, 'rgb(010,154,000)', 'rgb(046,189,070)', 'black',  .05],
    16: ['S',  2.15, 'rgb(238,206,000)', 'rgb(253,196,000)', 'black',  .00],
    17: ['Cl', 2.25, 'rgb(140,205,000)', 'rgb(143,200,085)', 'black', -.10],
    26: ['Fe', 2.25, 'rgb(140,140,140)', 'rgb(200,200,209)', 'black', -.10],
    35: ['Br', 2.35, 'rgb(230,085,000)', 'rgb(243,160,085)', 'black', -.10],
    53: ['I',  2.55, 'rgb(175,000,205)', 'rgb(215,000,165)', 'black',  .10],
    }

partial_charge = {
    -3: 'rgb(000,000,255)',
    -2: 'rgb(050,000,255)',
    -1: 'rgb(100,000,255)',
     1: 'rgb(255,000,100)',
     2: 'rgb(255,000,050)',
     3: 'rgb(255,002,000)',
     }

bond_type = {
    'SINGLE':1,
    'DOUBLE':2,
    'TRIPLE':3
    }

bond_arrange = {
        1: [0],
        2: [-1, 1],
        3: [-1.5, 0, 1.5]
        }


class atom:
    """ Represents an atom """

    def __init__(self, element, radius, pos, scale_factor=1.0):
        self.element = element
        self.radius = radius
        self.pos = pos
        self.scale_factor = scale_factor
        self.font_size = scale_factor * radius / 2.5

    def rotate(self, anchor, angle):
        relative_pos = self.pos - anchor
        self.pos = np.dot(rotate_matrix(angle), relative_pos) + anchor

    def scale_pos(self, anchor, N):
        self.new_pos = np.dot(scale_matrix(N), self.pos-anchor)
        self.pos = self.new_pos + anchor

    def mirror(self, angle):
        self.pos = np.dot(mirror_matrix(angle), self.pos)

    def draw(self, canvas):
        constants = atom_dict[self.element]
        [element_symbol, rad, main_color,
         grad_color, font_color, font_xoffset] = constants[:6]

        rad *= self.radius * self.scale_factor / 3.0
        gradients = [
                     {
                      'grad': canvas.defs.add(canvas.linearGradient((0,0), (0,1))),
                      'offset': [1.0, 1.0, 1.0],
                      'color': 3*[main_color],
                      'opacity': [1.0, 1.0, 1.0],
                      'position': (self.pos[0], self.pos[1]),
                      'radius': (rad, rad)
                     },
                     {
                      'grad': canvas.defs.add(canvas.linearGradient((0,0), (0,1))),
                      'offset': [.0, .85, 1.0],
                      'color': 3*['white'],
                      'opacity': [1.0, .0, .0],
                      'position': (self.pos[0], self.pos[1]-rad/2),
                      'radius': (rad/1.5, rad/3)
                     },
                     {
                      'grad': canvas.defs.add(canvas.linearGradient((0,0), (0,1))),
                      'offset': [.0, .6, 1.0],
                      'color': 2*['white'] + [grad_color],
                      'opacity': [.0, .0, 1.0],
                      'position': (self.pos[0], self.pos[1]+rad/2),
                      'radius': (rad/2, rad/2.5)
                     }
                    ]
        for g in gradients:
            for i in range(3):
                g['grad'].add_stop_color(offset=g['offset'][i], color=g['color'][i], opacity=g['opacity'][i])
            canvas.add(canvas.ellipse(center = g['position'],
                                      r = g['radius'],
                                      fill = g['grad'].get_paint_server()
                                      ))

        radius_vec = np.array([-self.radius/7.0, self.radius/7.0])
        label_pos = self.pos + radius_vec
        canvas.add(canvas.text(element_symbol,
                               insert=label_pos,
                               font_family = 'Arial',
                               font_size = self.font_size,
                               font_weight = 'bold',
                               fill = font_color))


class bond:
    def __init__(self, atom1, atom2, type, scale_factor):
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type
        self.scale_factor = scale_factor / 10.0

    def draw(self, canvas):
        dir_vector = self.atom2.pos - self.atom1.pos
        perp_vector = perp_2d_vec(dir_vector) * self.scale_factor
        lines = bond_arrange[self.type]

        for separation_factor in lines:
            start = self.atom1.pos + perp_vector * separation_factor
            end   = self.atom2.pos + perp_vector * separation_factor
            canvas.add(canvas.line(start = start,
                                   end   = end,
                                   stroke = 'black',
                                   stroke_width = 5))

class molecule:
    """ Represents a molecule """

    def __init__(self,
                 atoms=[], atom_radii=100.0,
                 bonds=[], scale_factor=1.0):
        self.atoms = atoms
        self.atom_radii = atom_radii
        self.bonds = bonds
        self.scale_factor = scale_factor
        self.center = None

    def create_from_smiles(self, smiles_code):
        # Creating base rdkit molecule object
        m = Chem.MolFromSmiles(smiles_code)
        m = AllChem.AddHs(m, False, False)
        confID = AllChem.Compute2DCoords(m, False, True)
        conf = m.GetConformer(confID)
        AllChem.WedgeMolBonds(m, conf)
        num_atoms= m.GetNumAtoms()

        # Getting all atomic coordinates and transforming
        for i in range(num_atoms):
            element = m.GetAtomWithIdx(i).GetAtomicNum()
            pos = np.array([m.GetConformer().GetAtomPosition(i).x,
                            m.GetConformer().GetAtomPosition(i).y])
            new_atom = atom(element, self.atom_radii, pos)
            self.atoms.append(new_atom)
        self.calculate_center()
        self.scale(self.scale_factor)

        # Getting all bonds with types
        for mbond in m.GetBonds():
            atom1 = self.atoms[mbond.GetBeginAtomIdx()]
            atom2 = self.atoms[mbond.GetEndAtomIdx()]
            type = mbond.GetBondTypeAsDouble()
            new_bond = bond(atom1, atom2, type, self.scale_factor)
            self.bonds.append(new_bond)

    def rotate(self, angle):
        for atom in self.atoms:
            atom.rotate(self.center, angle)

    def scale(self, N):
        for atom in self.atoms:
            atom.scale_pos(self.center, N)

    def mirror(self, angle):
        for atom in self.atoms:
            atom.mirror(angle)

    def draw(self, canvas):
        for bond in self.bonds:
            bond.draw(canvas)
        for atom in self.atoms:
            atom.draw(canvas)

    def find_minmax_points(self):
        min_x, max_x = 1E10, -1
        min_y, max_y = 1E10, -1
        for i, atom in enumerate(self.atoms):
            if atom.pos[0] < min_x:
                min_x = atom.pos[0]
            if atom.pos[0] > max_x:
                max_x = atom.pos[0]
            if atom.pos[1] < min_y:
                min_y = atom.pos[1]
            if atom.pos[1] > max_y:
                max_y = atom.pos[1]
        return [min_x, max_x, min_y, max_y]

    def calculate_center(self):
        minmax = self.find_minmax_points()
        self.c = np.array([[minmax[1]-minmax[0]],
                           [minmax[3]-minmax[2]]])
        self.center = self.c.transpose().flatten()

    def set_new_center(self, new_center):
        trans = new_center - self.center
        for atom in self.atoms:
            atom.pos = atom.pos + trans
        self.calculate_center()


def create_canvas(filename, size, debug=False):
    """ Returns an svg canvas object """

    canvas = svgwrite.Drawing(filename=filename+'.svg',
                              size=size,
                              debug=debug)
    return canvas
