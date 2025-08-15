micromamba create -n openff_ptm_prototype -c conda-forge openff-toolkit-examples python=3.11 pyxdg rustworkx
micromamba run -n openff_ptm_prototype pip install git+https://github.com/openforcefield/openff-pablo.git@main
micromamba run -n openff_ptm_prototype jupyter lab 2025_08_15.ipynb
