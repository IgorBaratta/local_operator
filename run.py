from utils import run_ffcx, create_ouput

import argparse
import yaml
import os
from subprocess import Popen, PIPE
import platform


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Run local assembly benchmark.')

    parser.add_argument('--form_compiler', dest='form_compiler', type=str,
                        default="ffcx", choices=['ffcx', 'ffc', 'tsfc'],
                        help="Form Compiler to use")

    parser.add_argument('--problem', dest='problem', type=str,
                        default="Lagrange", choices=['Lagrange', 'Elasticity', 'N1Curl', 'Stokes'],
                        help="Problem to run")

    parser.add_argument('--conf', dest='conf', type=str, default="compilers.yaml",
                        help="Configuration file describing the compilers and flags.")

    parser.add_argument('--degree', dest='degree', default=range(4), nargs='+',
                        help='Polynomial degree to evaluate the operators.')

    parser.add_argument('--nrepeats', dest='nrepeats', default=3, choices=range(1, 11),
                        help='Polynomial degree to evaluate the operators')

    args = parser.parse_args()
    form_compiler = args.form_compiler
    problem = args.problem
    conf_file = args.conf
    degrees = [int(d) for d in args.degree]
    nrepeats = args.nrepeats

    # Set architecture from platform
    try:
        with open("/sys/devices/cpu/caps/pmu_name", "r") as pmu:
            machine = pmu.readlines()[0].strip()
    except:
        machine = platform.processor()

    # Read Compiler configuration file
    with open(conf_file, "r") as stream:
        try:
            compilers = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    opt = ""

    out_file = create_ouput(problem)

    for c_name in compilers:
        compiler = compilers[c_name]
        flags = compiler["flags"]
        try:
            with Popen([compiler["cpp"][0], "-dumpversion"], stdout=PIPE) as p:
                compiler_version = p.stdout.read().decode("ascii").strip()
        except:
            compiler_version = compiler["version"][0]

        os.environ["CXX"] = compiler["cpp"][0]
        os.environ["CC"] = compiler["cc"][0]

        for flag in flags:
            flag = "\"" + ''.join(map(str, flag)) + "\""
            for degree in degrees:
                text = f"\n{machine}, {problem}, {c_name}, {compiler_version}, {flag}, {degree}, "
                results = run_ffcx(problem, degree, nrepeats, flag)
                for result in results:
                    row = text + f"\"{opt}\", {1}, {result}"
                    with open(out_file, "a") as file:
                        file.write(row)
