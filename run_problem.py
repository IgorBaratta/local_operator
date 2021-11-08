import os
import ffcx.codegeneration
import sys
import platform
from subprocess import Popen, PIPE
from string import Template

if len(sys.argv) > 1:
    problem = sys.argv[1]
else:
    problem = "Lagrange.ufl"


#########################
# COMPILERS AND FLAGS
#########################
compilers = [["g++-11", "gcc-11"], ["clang++-12", "clang-12"]]
opt_flags = ["\"-Ofast -march=native -mprefer-vector-width=256\"",
             "\"-Ofast -march=native -fno-tree-vectorize\"",
             "\"-O3 -march=native\""]


# Set architecture from platform
try:
    with open("/sys/devices/cpu/caps/pmu_name", "r") as pmu:
        machine = pmu.readlines()[0].strip()
except:
    machine = platform.processor()

os.environ["UFC_INCLUDE_DIR"] = ffcx.codegeneration.get_include_path()

family = problem.split(".")[0]
nrepeats = 3
degrees = [1, 2, 3, 4]

ffc_opts = {"ffcx": ""}

title = "machine,problem,compiler,version,flags,degree,method,rank,ncells,time"
out_file = "output/" + str(family) + ".txt"
if not os.path.exists(out_file):
    with open(out_file, "a") as f:
        f.write(title)

for flag in opt_flags:
    for compiler in compilers:
        for degree in degrees:
            try:
                with Popen([compiler[0], "-dumpversion"], stdout=PIPE) as p:
                    compiler_version = p.stdout.read().decode("ascii").strip()
            except:
                compiler_version = "unknown"
            os.environ["CXX"] = compiler[0]
            os.environ["CC"] = compiler[1]
            # Uses PE_ENV if on Cray
            compiler_name = os.environ.get("PE_ENV", compiler[0])

            d = {'degree': str(degree), 'vdegree': str(degree + 1)}
            with open("forms/" + problem, 'r') as f:
                src = Template(f.read())
                result = src.substitute(d)
                print(result)
                with open("problem.ufl", "w") as f2:
                    f2.writelines(result)
                with open("problem.py", "w") as f2:
                    f2.writelines(result)

            build = f"cd ffcx && rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"
            text = f"\n{machine}, {family}, {compiler_name}, {compiler_version}, {flag}, {degree}, "
            for opt in ffc_opts:
                print(f"ffcx {ffc_opts[opt]} problem.ufl")
                if os.system(f"ffcx {ffc_opts[opt]} problem.ufl -o ffcx/") != 0:
                    raise RuntimeError("ffcx failed")
                if os.system(build) != 0:
                    raise RuntimeError("build failed")

                for i in range(nrepeats):
                    # Rank 0
                    with Popen(["./ffcx/build/benchmark-mf"], stdout=PIPE) as p:
                        result = p.stdout.read().decode("ascii").strip()
                    text1 = text + f"\"{opt}\", {1}, {result}"

                    print(i, text1)
                    with open(out_file, "a") as file:
                        file.write(text1)

                    # Rank 1
                    with Popen(["./ffcx/build/benchmark"], stdout=PIPE) as p:
                        result = p.stdout.read().decode("ascii").strip()
                    text1 = text + f"\"{opt}\", {2}, {result}"

                    print(i, text1)
                    with open(out_file, "a") as file:
                        file.write(text1)
