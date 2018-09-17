#!/usr/bin/env python3

"""
Creates svg files of molecules from smiles code.

Further description will come when script is fully ready.
Written by Peleg Bar Sapir
"""

from lib.main import *
from lib.math import *
import rdkit
import svgwrite
import numpy as np

if __name__ == '__main__':
    print('hi')
    canvas = create_canvas('test', (600, 600))
    H1 = atom(0, 'H', 150.0, np.array([100, 500]))
    H2 = atom(1, 'H', 150.0, np.array([450, 350]))
    O =  atom(2, 'O', 150.0, np.array([150, 50]))
    N =  atom(3, 'N', 150.0, np.array([250, 400]))
    mol = molecule([H1, H2, O, N])
    mol.rotate(np.array([300, 300]), np.pi)
    mol.draw(canvas)
    canvas.save()
