import utils
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Run local assembly benchmark.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--form_compiler', dest='form_compiler', type=str,
                        default="ffcx", choices=['ffcx', 'ffc', 'tsfc'],
                        help="Form Compiler to use")

    parser.add_argument('--scalar_type', dest='scalar_type', type=str,
                        default="double", choices=['double', 'float', '_Float16', 'double _Complex', 'float _Complex'],
                        help="Scalar type to use")

    parser.add_argument('--problem', dest='problem', type=str,
                        default="Laplacian", choices=['Laplacian', 'Mass', 'Elasticity', 'N1curl', 'Stokes'],
                        help="Problem to run")

    parser.add_argument('--conf', dest='conf', type=str, default="compilers.yaml",
                        help="Configuration file describing the compilers and flags.")
    
    parser.add_argument('--output_file', dest='output_file', type=str, default="output/output.csv",
                        help="Configuration file describing the compilers and flags.")

    parser.add_argument('--degree', dest='degree', default=range(1, 4), nargs='+',
                        help='Polynomial degree to evaluate the operators.')

    parser.add_argument('--nrepeats', dest='nrepeats', type=int, default=3,
                        help='Number of times to run each experiment.')

    parser.add_argument('--batch_size', dest='batch_size', type=int, default=None, choices=[None, 1, 2, 4, 8, 16, 32, 64],
                        help='')

    parser.add_argument('--global_size', dest='global_size', type=int, default=1e6,
                        help='Global number of dofs (assuming shared are dofs are duplicated).')

    parser.add_argument('--action', dest='action', action='store_true',
                        help='Specify whether to run the problems with matrix free approach.')

    parser.add_argument('--mpi_size', dest='mpi_size', type=int, default=1,
                        help='The number of mpi processes to use.')

    parser.add_argument('--cell_type', dest='cell_type', type=str,
                        default="tetrahedron", choices=['tetrahedron', 'hexahedron'],
                        help="Cell type to use")

    args = parser.parse_args()
    form_compiler = args.form_compiler
    problem = args.problem
    conf_file = args.conf
    degrees = [int(d) for d in args.degree]
    nrepeats = args.nrepeats
    action = args.action
    batch_size = args.batch_size
    global_size = args.global_size
    scalar_type = args.scalar_type
    mpi_size = args.mpi_size
    cell_type = args.cell_type
    output_file = args.output_file

    machine = utils.machine_name()
    out_file = utils.create_output(problem, output_file)
    compilers = utils.parse_compiler_configuration(conf_file)

    # Set rank to 1 for matrix free, 2 otherwise
    rank = 1 if action else 2

    for c_name in compilers:
        compiler = compilers[c_name]
        compiler_version = utils.set_compiler(compiler)
        flags = compiler["flags"]
        for flag in flags:
            flag = "\"" + ''.join(map(str, flag)) + "\""
            for degree in degrees:
                text = f"\n{machine}, {problem}, {c_name}, {compiler_version}, {flag}, {degree}, {form_compiler}, {scalar_type}, {batch_size}, {cell_type}, "
                results = utils.run(problem, degree, nrepeats, flag, action,
                                    scalar_type, global_size, batch_size,
                                    mpi_size, cell_type)
                for result in results:
                    row = text + f"{rank}, {result}"
                    with open(out_file, "a") as file:
                        file.write(row)
