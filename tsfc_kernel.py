import problem
import tsfc

form = tsfc.compile_form(problem.a)[0]
form_int = form.ast.gencode()

d = problem.degree
if problem.family == "Lagrange":
    dofs = int((d+1)*(d+2)*(d+3)/6)
elif problem.family == "N1curl":
    dofs = int(d*(d+2)*(d+3)/2)

print(f"#define DOFS {dofs}")
print("#include <cmath>")
print("#define restrict __restrict")
print("using namespace std;")
print("\n")
print(form_int)
