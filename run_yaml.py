import yaml


file = "problem.yaml"
with open(file, "r") as stream:
    try:
        problems = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


for problem in problems:
    assert isinstance(problem, dict), 'Argument of wrong type!'
    fields = problem.keys()

    name = problem.get("name", "Laplacian")
    assert name in ['Laplacian', 'Mass', 'Elasticity', 'N1curl', 'Stokes']

    degrees = problem.get("degree", [1, 2, 3])
    assert isinstance(degrees, list)

    num_repeats = problem.get("num_repeats", 3)
    assert num_repeats > 0

    num_dofs = [float(e) for e in problem.get("num_dofs", 1)]
    assert isinstance(num_dofs, list)

    cell_type = problem.get("cell_type", "tetrahedron")
    assert cell_type in ['tetrahedron', 'hexahedron']

    mpi_processes = problem.get("mpi_processes", [1])
    assert isinstance(mpi_processes, list)

    compilers = problem.get("compiler")


    scalar_type_list = problem.get("scalar_type")
    for scalar_type in scalar_type_list:
        for degree in degrees:
            utils.compile_form()
            for compiler in compilers:
                print(compiler)






