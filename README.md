# FastPMRunner
 
 Python code for running a suite of FastPM simulations. Prepare for running FastPM simulations in Latin hypercube for training emulators.

## Installation for FastPM

[FastPM GitHub website](https://github.com/fastpm/fastpm) has explained several installation methods. Feel free to pick one you like.

For testing purpose, below uses the binary installation via Anaconda:
```python
conda create -n cfastpm
conda activate cfastpm

conda install -c bccp cfastpm nbodykit
```

In the `cfastpm` environment, there should be a `fastpm` binary to run the simulation.

To run `FastPMRunner`, we also need matplotlib for plotting the output power spectrum (to see if the code works!)
```
conda install matplotlib
```

An example code for running 10 simulations with hubble parameter sampling in a range [0.65, 0.75] is provided in `example/one_parameter.py`:
```bash
python -c "from examples.one_parameter import *; one_parameter_hubble()"
```

After it's finished, you will find 10 simulation outputs will be saved in `simulation_files/hubble_xxxx/` folder.

An example plot of how hubble parameter (h) changes the matter power spectrum is generated in `simulation_files/pk_one_param_hubble.pdf`:

![Screenshot from 2021-04-11 04-11-20](https://user-images.githubusercontent.com/23435784/114302045-28704280-9a7c-11eb-9e53-35bbee8439d6.png)
