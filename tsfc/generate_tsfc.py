import problem
import tsfc
import loopy

form, = tsfc.compile_form(problem.L, coffee=False)

headers = "#include <cmath>"
headers += "\n\nusing namespace std;\n"

kernel = loopy.generate_code_v2(form.ast).device_code()
kernel = kernel.replace("void", "static void")

with open("tsfc_kernel.cpp", "w") as file:
    file.write(headers)
    file.write(kernel)
