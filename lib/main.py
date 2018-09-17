"""
Main library for moldraw
"""

import svgwrite
import numpy as np
from lib.math import *

atom_dict = {
    'H':  [1,  1.00, 'rgb(182,225,255)', 'rgb(176,248,251)', 'black',  .01],
    'B':  [5,  1.30, 'rgb(255,158,158)', 'rgb(255,179,179)', 'black',  .02],
    'C':  [6,  1.23, 'rgb(128,128,128)', 'rgb(179,179,179)', 'black',  .00],
    'N':  [7,  1.54, 'rgb(128,000,128)', 'rgb(255,000,255)', 'white',  .00],
    'O':  [8,  1.77, 'rgb(212,000,000)', 'rgb(255,085,085)', 'black',  .00],
    'F':  [9,  1.80, 'rgb(098,146,000)', 'rgb(143,160,085)', 'black',  .05],
    'Si': [14, 1.85, 'rgb(215,254,220)', 'rgb(200,200,200)', 'black', -.05],
    'P':  [15, 1.95, 'rgb(010,154,000)', 'rgb(046,189,070)', 'black',  .05],
    'S':  [16, 2.15, 'rgb(238,206,000)', 'rgb(253,196,000)', 'black',  .00],
    'Cl': [17, 2.25, 'rgb(140,205,000)', 'rgb(143,200,085)', 'black', -.10],
    'Fe': [26, 2.25, 'rgb(140,140,140)', 'rgb(200,200,209)', 'black', -.10],
    'Br': [35, 2.35, 'rgb(230,085,000)', 'rgb(243,160,085)', 'black', -.10],
    'I':  [53, 2.55, 'rgb(175,000,205)', 'rgb(215,000,165)', 'black',  .10],
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


class atom:
    """ Represents an atom """

    def __init__(self, id, element, radius, pos, scale_factor=1.0):
        self.id = id
        self.element = element
        self.radius = radius
        self.pos = pos
        self.scale_factor = scale_factor
        self.font_size = scale_factor * radius / 2.5

    def rotate(self, anchor, angle):
        relative_pos = self.pos - anchor
        self.pos = np.dot(rotate_matrix(angle), relative_pos) + anchor 

    def scale(self, N):
        self.radius = self.radius * N
        self.font_size = self.radius / 2.5
        self.pos = np.dot(scale_matrix(N), self.pos)

    def mirror(self, angle):
        self.pos = np.dot(mirror_matrix(angle), self.pos)

    def draw(self, canvas):
        constants = atom_dict[self.element]
        [atomic_num, rad, main_color,
         grad_color, font_color, font_xoffset] = constants[:6]

        rad = self.radius * self.scale_factor / 3.0
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
        canvas.add(canvas.text(self.element,
                               insert=label_pos,
                               font_family = 'Arial',
                               font_size = self.font_size,
                               font_weight = 'bold',
                               fill = font_color))

class molecule:
    """ Represents a molecule """

    def __init__(self, atoms):
        self.atoms = atoms

    def rotate(self, anchor, angle):
        for atom in self.atoms:
            atom.rotate(anchor, angle)

    def scale(self, N):
        for atom in self.atoms:
            atom.scale(N)

    def mirror(self, angle):
        for atom in self.atoms:
            atom.mirror(angle)

    def draw(self, canvas):
        for atom in self.atoms:
            atom.draw(canvas)

def create_canvas(filename, size, debug=False):
    """ Returns an svg canvas object """

    canvas = svgwrite.Drawing(filename=filename+'.svg',
                              size=size,
                              debug=debug)
    return canvas
