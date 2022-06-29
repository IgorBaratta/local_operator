Python >=3.7



# Install dependencies ffcx
python3 -m venv env/ffcx
source env/ffcx/bin/activate

python3 -m pip install fenics-ffcx
python3 -m pip install pyyaml

# Install dependencies for running ffc experiments
```bash
python3 -m venv env/ffc
source env/ffc/bin/activate
```

```bash
python3 -m pip install fenics-ffc
python3 -m pip install pyyaml
```


# Install dependencies tsfc
```bash
git clone git@github.com:firedrakeproject/tsfc.git
python3 -m venv env/tsfc

source env/tsfc/bin/activate
python3 -m pip install -r tsfc/requirements.txt
python3 -m pip install tsfc

python3 -m pip install git+https://github.com/FEniCS/basix.git
python3 -m pip install git+https://github.com/FEniCS/ffcx.git

python3 -m pip install pyyaml
```

#

advisor --collect=survey --project-dir=./advi --search-dir src:r=. -- ./build/benchmark
advisor -collect tripcounts -flop -stacks --project-dir=./advi --search-dir src:r=. -- ./build/benchmark
advisor --report=roofline --with-stack --project-dir=./advi --report-output=./out/roofline.html

