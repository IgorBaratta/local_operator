import os
import ffcx.codegeneration
import sys
import platform
from string import Template

if len(sys.argv) > 1:
    problem = sys.argv[1]
else:
    problem = "Lagrange.ufl"


#########################
# COMPILERS AND FLAGS
#########################
compilers = [["g++", "gcc"], ["clang++", "clang"]]
opt_flags = ["\"-Ofast -march=native\""]
machine = None

if not machine:
    machine = platform.processor()

os.environ["UFC_INCLUDE_DIR"] = ffcx.codegeneration.get_include_path()

family = problem.split(".")[0]
nrepeats = 10
degrees = [1, 2, 3]
# if family == "Lagrange":
#    degrees = [1, 2, 3, 4, 5]


title = "machine,problem,compiler,flags,degree,method,ncells,time"
out_file = str(family) + ".txt"
if not os.path.exists(out_file):
    with open(out_file, "a") as f:
        f.write(title)

for flag in opt_flags:
    for compiler in compilers:
        for degree in degrees:
            os.environ["CXX"] = compiler[0]
            os.environ["CC"] = compiler[1]
            compiler_name = os.environ.get("PE_ENV", compiler[0]) # Uses PE_ENV if on Cray

            d = {'degree': str(degree), 'vdegree': str(degree + 1)}
            with open(problem, 'r') as f:
                src = Template(f.read())
                result = src.substitute(d)
                print(result)
                with open("problem.ufl", "w") as f2:
                    f2.writelines(result)

            build = f"rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"
            ffc_opts = {"ffcx": "",
                        "fused": "--fuse_loops",
                        "fused + full_tables": "--fuse_loops --full_tables",
                        "fused + hoist": "--fuse_loops --full_tables --code_hoisting problem.ufl"}

            text = f"\n{machine}, {family}, {compiler_name}, {flag}, {degree}, "
            for opt in ffc_opts:
                print(f"ffcx {ffc_opts[opt]} problem.ufl")
                if os.system(f"ffcx {ffc_opts[opt]} problem.ufl") !=0 : raise RuntimeError("ffcx failed")
                if os.system(build) !=0 : raise RuntimeError("build failed")

                for i in range(nrepeats):
                    text1 = text + f"\"{opt}\", "
                    print(i, text1)
                    with open(out_file, "a") as file:
                        file.write(text1)
                    if os.system(f"./build/benchmark >>{out_file}") != 0 : raise RuntimeError("benchmark failed")
