import os
import ffcx.codegeneration
import sys

if len(sys.argv) > 1:
    problem = int(sys.argv[1])
else:
    problem = 0

#########################
# COMPILERS AND FLAGS
#########################
compilers = [["g++-10", "gcc-10"], ["clang++-12", "clang-12"]]
opt_flags = ["\"-Ofast -march=native\""]


os.environ["UFC_INCLUDE_DIR"] = ffcx.codegeneration.get_include_path()

if problem == 0:
    forms = ["a = inner(grad(u), grad(v))*dx"]
    family = "Lagrange"
    degrees = [1, 2, 3, 4, 5]
    nrepeats = 10

elif problem == 1:
    forms = ["a = inner(curl(u), curl(v))*dx"]
    family = "N1curl"
    degrees = [1, 2, 3, 4]
    nrepeats = 10

elif problem == 2:
    forms = ["a = inner(k * grad(u), grad(v)) * dx"]
    family = "DG"
    degrees = [1, 2, 3]
    nrepeats = 10

else:
    raise RuntimeError("Problems form 0 to 2")

title = "method,compiler,flags,degree,ncells,time"
for d in degrees:
    out_file = str(family) + str(d)
    with open(out_file, "a") as file:
        file.write(title)

for flag in opt_flags:
    for compiler in compilers:
        for form in forms:
            for degree in degrees:
                out_file = str(family) + str(degree)
                os.environ["CXX"] = compiler[0]
                os.environ["CC"] = compiler[1]

                with open("problem.ufl", "r") as file:
                    lines = file.readlines()

                lines[0] = "degree = " + str(degree) + "\n"
                lines[1] = "family = \"" + family + "\"\n"
                lines[9] = form

                with open("problem.ufl", "w") as file:
                    file.writelines(lines)

                text = f"\n{compiler[0]}, {flag}, {degree}, "
                print(text)
                build = f"cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"
                os.system("ffcx problem.ufl")
                os.system("rm -rf build")
                os.system("mkdir build")

                os.system(build)
                for i in range(nrepeats):
                    text = f"\nffcx, {compiler[0]}, {flag}, {degree}, "
                    with open(out_file, "a") as file:
                        file.write(text)
                    os.system(f"./build/benchmark >>{out_file}")

                # run fused loops case
                os.system("ffcx --fuse_loops problem.ufl")
                os.system("rm -rf build")
                os.system("mkdir build")
                os.system(build)

                for i in range(nrepeats):
                    text = f"\nfused, {compiler[0]}, {flag}, {degree}, "
                    with open(out_file, "a") as file:
                        file.write(text)
                    os.system(f"./build/benchmark >>{out_file}")

                # run fused loops case
                os.system("ffcx --fuse_loops --full_tables problem.ufl")
                os.system("rm -rf build")
                os.system("mkdir build")
                os.system(build)

                for i in range(nrepeats):
                    text = f"\nfused + full_tables, {compiler[0]}, {flag}, {degree}, "
                    with open(out_file, "a") as file:
                        file.write(text)
                    os.system(f"./build/benchmark >>{out_file}")

                # run fused loops case
                os.system("ffcx --fuse_loops --full_tables --code_hoisting problem.ufl")
                os.system("rm -rf build")
                os.system("mkdir build")
                os.system(build)

                for i in range(nrepeats):
                    text = f"\nfused + hoist, {compiler[0]}, {flag}, {degree}, "
                    with open(out_file, "a") as file:
                        file.write(text)
                    os.system(f"./build/benchmark >>{out_file}")
