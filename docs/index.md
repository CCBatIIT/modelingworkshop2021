The workshop will be a one-week event consisting of 10 modules (~3 hours each) that include both theory and computer laboratory. Topics are as follows:

| Day | Module 1 | Module 2 |
| --- | -------- | -------- |
| 1 | **Structural analysis and visualization**. Protein structure and representation. Structural alignment. Homology modeling. Electrostatics calculations. | **Biomolecular potential energy functions** [Willow]. Molecular mechanics force fields. Mixed quantum mechanics/molecular mechanics. |
| 2 | **Molecular docking**. Scoring functions.  Common optimization algorithms. | **Molecular simulation** [Spiridon]. Markov chain Monte Carlo. The Metropolis-Hastings Algorithm. Hamiltonian Monte Carlo and molecular dynamics. Constrained dynamics including torsion and rigid-body dynamics. Integrators, thermostats, and barostats. |
| 3 | **Analysis of molecular simulations**. Equilibration versus production. Visualization of trajectories. Time series, averages, and histograms of properties including RMSD, potential energy, and distances. Principal components analysis. | **Markov state models**. Conformational clustering. Markov chains. Microstates and empirical transition matrices. Computing equilibrium populations and kinetics. |
| 4 | **Simulating thermodynamic processes**. Thermodynamic cycles. Umbrella sampling. Replica exchange. | **Analysis of thermodynamic process simulations**. Potential of mean force. Statistical estimators. |
| 5 | **Binding free energy calculations**. Alchemical methods. | **Analysis of binding free energy calculations**. Quality metrics. Pose prediction. |

Lectures will be recorded and computer laboratory exercises will be posted online for participants to complete asynchronously and as a free online resource.

## Preparation

Participants should have their own laptop computers capable of running virtual machines. A virtual machine that includes all software necessary for the workshop will be provided by the instructors. To use the virtual machine, first install [Virtual Box](https://www.virtualbox.org/).

The virtual machine includes:
- [Visual Molecular Dynamics (VMD)](https://www.ks.uiuc.edu/Research/vmd/) for visualizing biological macromolecules and molecular dynamics trajectories.
- [AutoDock Tools & AutoDock Vina](http://autodock.scripps.edu/) for performing molecular docking.
- [Robosample](https://github.com/spirilaurentiu/Robosample) for performing constrained molecular dynamics.
- [Anaconda python 3](https://www.anaconda.com/). The Python programming language is widely used in the scientific community. It includes many packages for data analysis and visualization. Anaconda provides a high-quality management of virtual environments. The virtual machine includes the following packages in two conda environments, _msmbuilder_ and _analysis_ (which contains everything else):
  - [jupyter notebook](https://jupyter.org/) to create and view documents that contain live code, equations, visualizations, and narrative text.
  - [matplotlib](https://matplotlib.org/) to create data visualizations.
  - [OpenMM](http://openmm.org/) to perform molecular simulations.
  - [pySCF](https://sunqm.github.io/pyscf/) to perform quantum chemistry calculations.
  - [ProDy](http://prody.csb.pitt.edu/) to analyze protein structural dynamics.
  - [pymbar](https://pymbar.readthedocs.io/en/master/), implementing the multistate Bennett acceptance ratio method for estimating expectation values and free energy differences.
  - [yank](http://getyank.org/latest/) for performing alchemical free energy calculations.
  - [msmbuilder](http://msmbuilder.org/3.8.0/) to create statistical models for biomolecular dynamics.

We will also use these web servers:
- Homology modelling
  - [I-TASSER](https://zhanglab.ccmb.med.umich.edu/I-TASSER/), one of most accurate in competitions
  - [Swiss Model](https://swissmodel.expasy.org/), one of the fastest
- [Advanced Poisson-Boltzmann Solver (APBS)](https://server.poissonboltzmann.org/) server

# Who?

This workshop will introduce advanced undergraduate students, graduate students, and research scientists in chemistry, biology, and related fields to computational methods for modeling biological macromolecules. Workshop participants will mostly be from Vietnam, the United States of America, and Romania. The workshop will be primarily taught by David Minh with guest lectures from Laurentiu Spiridon and Soohaeng Yoo Willow.

Up to 25 participants will join the workshop in person, including up to 10 students as part of a class at the Illinois Institute of Technology. An additional 25 students may participate online and receive support from the instructors. If the workshop is held entirely online due to COVID-19, there will be up to 50 total participants.

# When?

The workshop will be held from June 28 to July 2, 2021. Modules will start at the following times, which are in Indochina Time (ICT), Central Daylight Time (CDT), and Eastern European Standard Time (EEST):

| Module | ICT | CDT | EEST |
| ------ | --- | --- | ---- |
| 1 | 8 am | 8 pm (-1 day) | 4 am EEST |
| 2 | 2 pm | 2 am CDT | 10 am EEST |

Workshop instructors will hold online office hours at the following times:

| Instructor | ICT | CDT | EEST |
| ---------- | --- | --- | ---- |
| David, starting day 2 | 7 am | 7 pm (-1 day) | 3 am |
| Laurentiu | 1 pm | 1 am | 9 am |
| Soohaeng | 1 am | 1 pm (-1 day) | 9 pm |

All participants, in-person and remote, are welcome to work through the modules synchronously. However, it will probably make sense for remote participants, based on their time zone, to complete modules at a more convenient time. Participants in CDT are encouraged to join module 1 live and complete module 2 the following morning. Participants in EEST are encouraged to complete module 1 immediately before participating live in module 2. All participants can visit all workshop instructors for online office hours.

A Chicago-based practice workshop will be held from January 11 to 15, 2021. Modules will start at the following times, which are in Indochina Time (ICT), Central Standard Time (CST), and Eastern European Time (EET):

| Module | ICT | CST | EEST |
| ------ | --- | --- | ---- |
| 1 | 2 am | 1 pm (-1 day) | 9 pm |
| 2 | 9 pm | 8 am | 4 pm |

# Where?

If international travel returns to normal, the workshop will be held at the University of Nha Trang in Nha Trang, Vietnam. The workshop will be held in an internet-connected classroom with a projection system and capacity for about 25 students. Room furniture will be set up in a way to facilitate instructors roaming the room and students gathering in groups.

If travel to Vietnam continues to be restricted due to COVID-19, the workshop will be held online.

# Why?

The purpose of this workshop will be
- to further the education of workshop participants in the field of computational biochemistry;
- to promote international scientific collaboration, especially between the United States, Vietnam, and Romania; and
- to provide a rich cultural exchange experience that will foster this collaboration.
