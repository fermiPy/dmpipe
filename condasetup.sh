# Setup script for conda-based ST installations

unset LD_LIBRARY_PATH
unset DYLD_LIBRARY_PATH
unset PYTHONPATH

export PATH=$HOME/miniconda/bin:$PATH

if [[ -n $CONDA_ENV_PATH ]]; then
   export PATH=$CONDA_ENV_PATH/bin:$PATH
fi

