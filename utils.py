import os
import platform
import yaml
from string import Template
from subprocess import Popen, PIPE

import sys


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


def create_ouput(problem):
    header = "machine,problem,compiler,version,flags,degree,fcomp,rank,ncells,time"
    path = "output/"
    out_file = path + str(problem) + ".txt"

    if not os.path.exists(out_file):
        if not os.path.isdir(path):
            os.mkdir(path)
        with open(out_file, "a") as f:
            f.write(header)
    return out_file


def run_ffcx(problem: str, degree: int, nrepeats: int,
             flag: list, matrix_free: bool):
    try:
        import ffcx
        import ffcx.codegeneration
    except ImportError:
        print("ffcx is no available")

    with open("forms/" + problem + ".ufl", 'r') as f:
        src = Template(f.read())
        d = {'degree': str(degree), 'vdegree': str(degree + 1)}
        result = src.substitute(d)

        with open("ffcx/problem.py", "w") as f2:
            f2.writelines(result)

    sys.path.insert(1, 'ffcx/')
    from compile import generate_code
    generate_code(matrix_free)

    run = "./ffcx/build/benchmark"
    build = f"cd ffcx && rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"

    if os.system(build) != 0:
        raise RuntimeError("build failed")
    result = []
    for i in range(nrepeats):
        with Popen([run], stdout=PIPE) as p:
            out = p.stdout.read().decode("ascii").strip()
        result.append(out)
    print(result)

    return result


def run_tsfc(problem: str, degree: int, nrepeats: int,
             flag: list, matrix_free: bool):
    try:
        import tsfc
        import ffcx
    except ImportError:
        print("tsfc is no available")

    with open("forms/" + problem + ".ufl", 'r') as f:
        src = Template(f.read())
        d = {'degree': str(degree), 'vdegree': str(degree + 1)}
        result = src.substitute(d)
        with open("tsfc/problem.py", "w") as f2:
            f2.writelines(result)

    sys.path.insert(1, 'tsfc/')
    from generate_tsfc import generate_code
    generate_code(matrix_free)

    run = "./tsfc/build/benchmark"
    build = f"cd tsfc && rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"

    if os.system(build) != 0:
        raise RuntimeError("build failed")
    result = []
    for i in range(nrepeats):
        with Popen([run], stdout=PIPE) as p:
            out = p.stdout.read().decode("ascii").strip()
        result.append(out)
    print(result)

    return result


def run_ffc(problem: str, degree: int, nrepeats: int,
            flag: list, matrix_free: bool):
    try:
        import ffc
    except ImportError:
        print("ffc is no available")

    with open("forms/" + problem + ".ufl", 'r') as f:
        src = Template(f.read())
        d = {'degree': str(degree), 'vdegree': str(degree + 1)}
        result = src.substitute(d)

        with open("ffc/problem.py", "w") as f2:
            f2.writelines(result)
        with open("ffc/problem.ufl", "w") as f2:
            f2.writelines(result)

    sys.path.insert(1, 'ffc/')
    from compile_ffc import generate_code
    generate_code(matrix_free)

    run = "./ffc/build/benchmark"
    build = f"cd ffc && rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"

    if os.system(build) != 0:
        raise RuntimeError("build failed")
    result = []
    for i in range(nrepeats):
        with Popen([run], stdout=PIPE) as p:
            out = p.stdout.read().decode("ascii").strip()
        result.append(out)
    print(result)

    return result


def run_cross(problem: str, degree: int, nrepeats: int,
              flag: list, matrix_free: bool, batch_size: int):
    try:
        import ffcx
        import ffcx.codegeneration
    except ImportError:
        print("ffcx is no available")

    with open("forms/" + problem + ".ufl", 'r') as f:
        src = Template(f.read())
        d = {'degree': str(degree), 'vdegree': str(degree + 1)}
        result = src.substitute(d)

        with open("cross/problem.py", "w") as f2:
            f2.writelines(result)

    sys.path.insert(1, 'cross/')
    from compile import generate_code
    generate_code(matrix_free, batch_size)

    run = "./cross/build/benchmark"
    build = f"cd cross && rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"

    if os.system(build) != 0:
        raise RuntimeError("build failed")
    result = []
    for i in range(nrepeats):
        with Popen([run], stdout=PIPE) as p:
            out = p.stdout.read().decode("ascii").strip()
        result.append(out)
    print(result)

    return result
