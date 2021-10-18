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
compilers = [["g++", "gcc"]]
opt_flags = ["\"-Ofast -march=native\""]

# Set machine name, or leave as None to get architecture from platform
machine = None
if not machine:
    machine = platform.processor()

os.environ["UFC_INCLUDE_DIR"] = ffcx.codegeneration.get_include_path()

family = problem.split(".")[0]
nrepeats = 3
degrees = [1, 2, 3, 4]

ffc_opts = {"ffcx": ""}

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
            # Uses PE_ENV if on Cray
            compiler_name = os.environ.get("PE_ENV", compiler[0])

            d = {'degree': str(degree), 'vdegree': str(degree + 1)}
            with open(problem, 'r') as f:
                src = Template(f.read())
                result = src.substitute(d)
                print(result)
                with open("problem.ufl", "w") as f2:
                    f2.writelines(result)
                with open("problem.py", "w") as f2:
                    f2.writelines(result)

            # Generate TSFC kernels
            os.system('python3 tsfc_kernel.py')

            build = f"rm -rf build && mkdir build && cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"
            text = f"\n{machine}, {family}, {compiler_name}, {flag}, {degree}, "
            for opt in ffc_opts:
                print(f"ffcx {ffc_opts[opt]} problem.ufl")
                if os.system(f"ffcx {ffc_opts[opt]} problem.ufl") != 0:
                    raise RuntimeError("ffcx failed")
                if os.system(build) != 0:
                    raise RuntimeError("build failed")

                for i in range(nrepeats):
                    text1 = text + f"\"{opt}\", "
                    print(i, text1)
                    with open(out_file, "a") as file:
                        file.write(text1)
                    if os.system(f"./build/benchmark >>{out_file}") != 0:
                        raise RuntimeError("benchmark failed")
                    text2 = text + f"\"tsfc\", "
                    print(i, text2)
                    with open(out_file, "a") as file:
                        file.write(text2)
                    if os.system(f"./build/benchmark >>{out_file} 1") != 0:
                        raise RuntimeError("benchmark failed")
