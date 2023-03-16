from .compile import generate_code

try:
    from . import problem
except:
    f = open("problem.py", 'w+')