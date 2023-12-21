import os
import platform
import yaml
from string import Template
from subprocess import Popen, PIPE

import sys

_build_cmd = "cd {form_compiler} && rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"


def set_compiler(compiler):
    os.environ["CXX"] = compiler["cpp"][0]
    os.environ["CC"] = compiler["cc"][0]

    try:
        with Popen([compiler["cpp"][0], "-dumpversion"], stdout=PIPE) as p:
            compiler_version = p.stdout.read().decode("ascii").strip()
    except:
        compiler_version = compiler["version"][0]
    return compiler_version


def parse_compiler_configuration(file):
    # Read Compiler configuration file
    with open(file, "r") as stream:
        try:
            compilers = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return compilers


def machine_name():
    # Set architecture from platform
    try:
        with open("/sys/devices/cpu/caps/pmu_name", "r") as pmu:
            machine = pmu.readlines()[0].strip()
    except:
        machine = platform.processor()
    return machine


def create_output(problem, out_file):
    header = "machine,problem,compiler,version,flags,degree,fcomp,scalar,batch_size,rank,cell_type,ncells,time"

    if not os.path.exists(out_file):
        with open(out_file, "a") as f:
            f.write(header)
    return out_file


def run(problem: str, degree: int, nrepeats: int, flag: str, action: bool,
        scalar_type: str, global_size: int, batch_size: int, mpi_size: int,
        cell_type: str):

    try:
        import ffcx
        import ffcx.codegeneration
    except ImportError:
        print("ffcx is not available")

    with open("forms/" + problem + ".ufl", 'r') as f:
        src = Template(f.read())
        d = {'degree': str(degree), 'vdegree': str(
            degree + 1), "cell": cell_type}
        result = src.substitute(d)

        with open("ffcx/problem.py", "w") as f2:
            f2.writelines(result)

    sys.path.insert(1, 'ffcx/')
    from compile import generate_code
    generate_code(action, scalar_type, global_size, batch_size)

    run = f"mpirun -n {mpi_size} ./ffcx/build/benchmark"
    build = _build_cmd.format(form_compiler="ffcx", flag=flag)

    if os.system(build) != 0:
        raise RuntimeError("build failed")
    result = [Popen(run.split(), stdout=PIPE).stdout.read().decode("ascii").strip() for i in range(nrepeats)]
    print(result)

    return result
