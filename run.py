import utils
import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Run local assembly benchmark.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--form_compiler', dest='form_compiler', type=str,
                        default="ffcx", choices=['ffcx', 'ffc', 'tsfc'],
                        help="Form Compiler to use")

    parser.add_argument('--problem', dest='problem', type=str,
                        default="Lagrange", choices=['Lagrange', 'Elasticity', 'N1curl', 'Stokes'],
                        help="Problem to run")

    parser.add_argument('--conf', dest='conf', type=str, default="compilers.yaml",
                        help="Configuration file describing the compilers and flags.")

    parser.add_argument('--degree', dest='degree', default=range(1, 4), nargs='+',
                        help='Polynomial degree to evaluate the operators.')

    parser.add_argument('--nrepeats', dest='nrepeats', type=int, default=3, choices=range(1, 11),
                        help='Polynomial degree to evaluate the operators')

    parser.add_argument('--matrix_free', dest='mf', action='store_true',
                        help='Specify whether to run the problems with matrix freee approach.')

    args = parser.parse_args()
    form_compiler = args.form_compiler
    problem = args.problem
    conf_file = args.conf
    degrees = [int(d) for d in args.degree]
    nrepeats = args.nrepeats
    mf = args.mf

    machine = utils.machine_name()
    out_file = utils.create_ouput(problem)
    compilers = utils.parse_compiler_configuration(conf_file)

    # Set rank to 1 for matrix free, 2 otherwise
    rank = 1 if mf else 2

    # TODO: Add ffcx options
    opt = ""

    for c_name in compilers:
        compiler = compilers[c_name]
        compiler_version = utils.set_compiler(compiler)
        # Get set of flags to use
        flags = compiler["flags"]
        for flag in flags:
            flag = "\"" + ''.join(map(str, flag)) + "\""
            for degree in degrees:
                text = f"\n{machine}, {problem}, {c_name}, {compiler_version}, {flag}, {degree}, {form_compiler}, "
                if form_compiler == "ffcx":
                    results = utils.run_ffcx(
                        problem, degree, nrepeats, flag, mf)
                elif form_compiler == "tsfc":
                    results = utils.run_tsfc(
                        problem, degree, nrepeats, flag, mf)
                else:
                    results = utils.run_ffc(
                        problem, degree, nrepeats, flag, mf)

                for result in results:
                    row = text + f"{rank}, {result}"
                    with open(out_file, "a") as file:
                        file.write(row)
