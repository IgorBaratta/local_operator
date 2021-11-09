import os
from string import Template
from subprocess import Popen, PIPE


def create_ouput(problem):
    header = "machine,problem,compiler,version,flags,degree,method,rank,ncells,time"
    path = "output/"
    out_file = path + str(problem) + ".txt"

    if not os.path.exists(out_file):
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

    os.environ["UFC_INCLUDE_DIR"] = ffcx.codegeneration.get_include_path()

    with open("forms/" + problem + ".ufl", 'r') as f:
        src = Template(f.read())
        d = {'degree': str(degree), 'vdegree': str(degree + 1)}
        result = src.substitute(d)
        with open("problem.ufl", "w") as f2:
            f2.writelines(result)

    build = f"cd ffcx && rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"

    if os.system(f"ffcx problem.ufl -o ffcx/") != 0:
        raise RuntimeError("ffcx failed")
    if os.system(build) != 0:
        raise RuntimeError("build failed")
    result = []
    for i in range(nrepeats):
        with Popen(["./ffcx/build/benchmark-mf"], stdout=PIPE) as p:
            out = p.stdout.read().decode("ascii").strip()
        result.append(out)

    return result
