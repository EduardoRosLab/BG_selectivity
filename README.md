
This repository contains the modules and experiments used for the paper "A basal ganglia computational model to explain the paradoxical sensorial improvement in the presence of Huntington's disease" (not yet published) by Álvaro González-Redondo, Francisco Naveros, Eduardo Ros and Jesús A. Garrido.

# Description of the model

The model includes five neuronal populations and nine neuronal types (all of them implemented as Izhikevich neuron models, but with different parameters in order to capture their particular cell dynamics). The total number of simulated neurons is 5,494 divided as follows: the MSN layer contains 2,292 neurons, with half of them (1,146) expressing D1 receptor and the other half D2 receptor. The STN layer contains 47 neurons, the GPe 155 neurons, and the SNr 3,000 neurons.

The neuron populations in our BG model have been connected following a channel structure. As a general norm, the neurons in every channel are only allowed to synapse neurons in the same channel.

More details can be found in the paper. 

# Running the model

## Compiling and installing the NEST module

From the repository folder:

```
export NEST_INSTALL_DIR=/path/to/nest/
cd nest_modules
mkdir mb
cd mb
cmake -Dwith-nest=${NEST_INSTALL_DIR}/bin/nest-config ../IzhikevichCond
make
make install
```

If everything went right you can load the NEST module in PyNEST by doing `nest.InstallModule('izhikevich_cond_module')`.


## Notebooks

There are three notebooks. The first one make a single run and show the results. The second notebook run all the experiments with every value of DA and HD (it can take more than a day), and then save the data. The third notebook loads this data and runs ICA on it.
