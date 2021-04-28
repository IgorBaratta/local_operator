import os
import ffcx.codegeneration


forms = ["a = inner(grad(u), grad(v)) * dx"]
degrees = [1, 2, 3, 4]
nrepeats = 3
out_file = "results"

compilers = [["g++", "gcc"], ["clang++", "clang"]]
opt_flags = ["\"-O3 \"", "\"-O3 -march=native \"",
             "\"-Ofast -march=native \"", "\"-Ofast \"", 
             "\"-O2 -march=native \"", "\"-O2 \""]

os.environ["UFC_INCLUDE_DIR"] = ffcx.codegeneration.get_include_path()

for flag in opt_flags:
    for compiler in compilers:
        for form in forms:
            for degree in degrees:
                os.environ["CXX"] = compiler[0]
                os.environ["CC"] = compiler[1]
                print(compiler[1])

                with open("lagrange.ufl", "r") as file:
                    lines = file.readlines()

                lines[0] = "degree = " + str(degree) + "\n"
                lines[8] = form

                with open("lagrange.ufl", "w") as file:
                    file.writelines(lines)

                text = f"\n{compiler[0]}, {flag}, {degree}, "
                print(text)

                os.system("ffcx lagrange.ufl")
                os.system("rm -rf build")
                os.system("mkdir build")
                build = f"cd build && cmake -DCMAKE_C_FLAGS={flag} -DCMAKE_CXX_FLAGS={flag} .. && make"
                os.system(build)
                for i in range(nrepeats):
                    text = f"\nffcx, {compiler[0]}, {flag}, {degree}, "
                    with open(out_file, "a") as file:
                        file.write(text)
                    os.system(f"./build/benchmark >>{out_file}")

                # run fused loops case
                os.system(f"cp lagrange{degree}.c lagrange.c")
                os.system("rm -rf build")
                os.system("mkdir build")
                os.system(build)

                for i in range(nrepeats):
                    text = f"\nfused, {compiler[0]}, {flag}, {degree}, "
                    with open(out_file, "a") as file:
                        file.write(text)
                    os.system(f"./build/benchmark >>{out_file}")
