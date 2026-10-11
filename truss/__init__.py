"""2D truss finite element solver (direct stiffness method)."""
from .model import Truss
from .solver import Result, assemble_stiffness, solve
