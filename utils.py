import os
import platform
from string import Template
from subprocess import Popen, PIPE
from compiler import generate_code

import sys


def machine_name():
    # Set architecture from platform
    try:
        with open("/sys/devices/cpu/caps/pmu_name", "r") as pmu:
            machine = pmu.readlines()[0].strip()
    except:
        machine = platform.processor()
    return machine


def create_ouput(problem, header):
    path = "output/"
    out_file = path + str(problem) + ".txt"

    if not os.path.exists(out_file):
        if not os.path.isdir(path):
            os.mkdir(path)
        with open(out_file, "a") as f:
            f.write(header + "\n")
    return out_file


def compile_form(problem):
    try:
        import ffcx
        import ffcx.codegeneration
    except ImportError:
        print("ffcx is not available")

    with open("forms/" + problem.problem + ".ufl", 'r') as f:
        src = Template(f.read())
        d = {'degree': str(problem.degree),
             'vdegree': str(problem.degree + 1),
             'cell': problem.cell_type}
        result = src.substitute(d)
        print("==========================")
        print(result)
        print("==========================")

        with open("compiler/problem.py", "w") as f2:
            f2.writelines(result)
    action = problem.form_rank == 1
    generate_code(action, problem.scalar_type,
                  problem.num_dofs, problem.batch_size)


_build_cmd = "cd compiler && rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS='{flag}' -DCMAKE_CXX_FLAGS='{flag}' .. && make"


def compile_cpp(flag):
    build_cmd = _build_cmd.format(form_compiler="ffcx", flag=flag)
    if os.system(build_cmd) != 0:
        raise RuntimeError("build failed")


def run_code(mpi_size):
    run = f"mpirun -n {mpi_size} ./compiler/build/benchmark"
    with Popen(run.split(), stdout=PIPE) as p:
        out = p.stdout.read().decode("ascii").strip()
    out = out.split(", ")
    out = [float(o) for o in out]
    return out
