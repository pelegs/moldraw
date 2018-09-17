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
    mol = molecule(scale_factor=105.0)
    mol.create_from_smiles('')
    mol.rotate(np.pi/2)
    canvas = create_canvas('test', size=[800,800])
    mol.draw(canvas)
    canvas.save()
