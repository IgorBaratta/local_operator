import yaml
import utils
import utils
import argparse
from dataclasses import dataclass


@dataclass
class Results:
    machine: str | None = None
    problem: str | None = None
    compiler: str | None = None
    flags: str | None = None
    scalar_type: str | None = None
    batch_size: int | None = None
    form_rank: int | None = None
    degree: int | None = None
    cell_type: str | None = None
    num_dofs: int | None = None
    mpi_processes: int | None = None
    time: float | None = None
    num_cells: int | None = None

    def flush(self):
        out = ", ".join(str(value) for _, value in self.__dict__.items())
        return out

    def header(self):
        out = ",".join(str(key) for key, _ in self.__dict__.items())
        return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Run local assembly benchmark.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--file_name', dest='file_name', type=str, default="problem.yaml",
                        help="Configuration file describing the problem.")

    args = parser.parse_args()
    file_name = args.file_name

    with open(file_name, "r") as stream:
        try:
            problems = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for problem in problems:

        assert isinstance(problem, dict), 'Argument of wrong type!'
        fields = problem.keys()

        benchmark = Results()

        name = problem.get("name", "Laplacian")
        out_file = utils.create_ouput(name, benchmark.header())
        benchmark.problem = name
        assert name in ['Laplacian', 'Mass', 'Elasticity', 'Helmholtz']

        degrees = problem.get("degree", [1, 2, 3])
        assert isinstance(degrees, list)

        num_repeats = problem.get("num_repeats", 3)
        assert num_repeats > 0

        num_dofs = [float(e) for e in problem.get("num_dofs", 1)]
        assert isinstance(num_dofs, list)

        cell_type = problem.get("cell_type", "tetrahedron")
        benchmark.cell_type = cell_type
        assert cell_type in ['tetrahedron', 'hexahedron']

        mpi_processes = problem.get("mpi_processes", [1])
        assert isinstance(mpi_processes, list)

        action = problem.get("action", True)
        assert isinstance(action, bool)
        benchmark.form_rank = 1 if action else 2

        compilers = problem.get("compiler")

        scalar_type_list = problem.get("scalar_type")

        benchmark.machine = utils.machine_name()

        for type_t in scalar_type_list:
            benchmark.batch_size = type_t["batch_size"]
            benchmark.scalar_type = type_t["type"]
            for degree in degrees:
                benchmark.degree = degree
                for n in num_dofs:
                    benchmark.num_dofs = n
                    utils.compile_form(benchmark)
                    for compiler in compilers:
                        flags = compiler["flags"]
                        benchmark.compiler = utils.set_compiler(compiler)
                        for flag in flags:
                            utils.compile_cpp(flag)
                            benchmark.flags = flag
                            for i in range(num_repeats):
                                for num_procs in mpi_processes:
                                    benchmark.mpi_processes = num_procs
                                    benchmark.num_cells, benchmark.time = utils.run_code(
                                        num_procs)
                                    row = benchmark.flush()
                                    with open(out_file, "a") as file:
                                        file.write(row + "\n")
