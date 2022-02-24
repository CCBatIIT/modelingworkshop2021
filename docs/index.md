---
layout: default
title: About
---

This page describes a workshop in June 2021. For information about more recent or upcoming workshops, go to the [home page](https://ccbatiit.github.io/modelingworkshop/index.html).

# What?

The workshop will be a one-week event comprising 10 periods (~3 hours each) that include both theory and computer laboratory. If you would like to be a formal participant in the workshop, please submit an [application](https://forms.gle/Levf3mRaPi5efjhg8).

The workshop will be primarily taught by David Minh with guest lectures from Andrew Howard, Laurentiu Spiridon, and Soohaeng Yoo Willow. Topics are as follows:

| Day | Period 1 | Period 2 |
| --- | -------- | -------- |
{% for day in (1..5) %} | {{ day }} | {% for period in (1..2) %} {% for module in site.data.modules2021 %} {% if module.day == day %} {% if module.period == period %}
{%- if module.bold %}<b>{% endif %}{{ module.title }} {% if module.bold %}</b>{% endif -%}
{% if module.teacher %}({{ module.teacher }}){% endif -%}
{% if module.description %}. {{ module.description }} {% else %} {% endif -%}
[{% if module.slides == "ppt" %}[ppt](https://github.com/CCBatIIT/modelingworkshop/raw/main/slides/{{ module.basename }}.ppt)/{% elsif module.slides == "pdf"%}{% else %}[key](https://github.com/CCBatIIT/modelingworkshop/raw/main/slides/{{ module.basename }}.key)/{% endif -%}[pdf](https://github.com/CCBatIIT/modelingworkshop/raw/main/slides/{{ module.basename }}.pdf)].
{%- endif %} {% endif %} {% endfor %} | {% endfor %}
{% endfor %}

Lectures may be recorded and computer laboratory exercises may be posted online for participants to complete asynchronously and as a free online resource.

During the workshop, we will be using the following software:
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

For your convenience, David's PhD student Jim Tufts has prepared [a virtual machine that contains all of this software](https://jimtufts.com/data/pdb/LinuxMint.tar.gz) and a [compiled version of Robosampler](https://jimtufts.com/data/pdb/Robosample.tar.gz). To use this virtual machine, you will first need to install [Virtual Box](https://www.virtualbox.org/). The user name is _chemuser_ and initial password is _vm_.

We will also use these web servers:
- Homology modelling
  - [I-TASSER](https://zhanglab.ccmb.med.umich.edu/I-TASSER/), one of most accurate in competitions
  - [Swiss Model](https://swissmodel.expasy.org/), one of the fastest
- [Advanced Poisson-Boltzmann Solver (APBS)](https://server.poissonboltzmann.org/) server

# Who?

This workshop will introduce advanced undergraduate students, graduate students, and research scientists in chemistry, biology, and related fields to computational methods for modeling biological macromolecules. Workshop participants will mostly be from the United States of America, Romania, and Vietnam.

<!--Up to 25 participants will join the workshop in person, including up to 10 students as part of a class at the Illinois Institute of Technology.-->
Up to 15 students may participate as part of a class at the Illinois Institute of Technology: Chem 456, Computational Biochemistry and Drug Design. An additional 25 students may participate online and receive support from the instructors.
<!--If the workshop is held entirely online due to COVID-19, there will be up to 50 total participants.-->

Participants should have a laptop computer capable of running the virtual machine.

If you would like to be a formal participant in the workshop, please submit an [application](https://forms.gle/Levf3mRaPi5efjhg8). Formal participants will be able to attend the workshop live and participate in live discussions. Priority for formal participation will be given to Illinois Tech students and individuals working with a workshop instructor or collaborator.

After the workshop, lecture recordings and slides may posted online.

# When?

The workshop will be held from June 28 to July 2, 2021. Modules will start at the following times, which are in Central Daylight Time (CDT), Eastern European Standard Time (EEST), and Indochina Time (ICT):

| Module | CST  | EEST | ICT  |
| ------ | ---- | ---- | ---- |
| 1      | 9 am | 5 pm | 9 pm |
| 2      | 1 pm | 9 pm | 1 am (+1 d) |

All participants are welcome to work through the modules synchronously. Participants in CDT are especially encouraged to participate live. However, it will probably make sense for remote participants, based on their time zone, to complete modules at a more convenient time. All participants can reach out to all workshop instructors for online office hours. Online office hours will be held at the following times:

Workshop instructors will hold online office hours at the following times:

| Instructor | CDT  | EEST | ICT  |
| ---------- | ---- | ---- | ---- |
| David      | 7 pm | 3 am (+1 d) | 7 am (+1 d) |
| Laurentiu  | 1 am | 9 am | 1 pm |
| Soohaeng   | 8 am | 4 pm | 8 pm |

<!--
These were originally planned times for the workshop in Nha Trang:
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

These were times for the practice workshop in 1/2021:
| Module | ICT | CST | EEST |
| ------ | --- | --- | ---- |
| 1 | 2 am | 1 pm (-1 day) | 9 pm |
| 2 | 9 pm | 8 am | 4 pm |
-->

# Where?

The workshop will be held on the Illinois Tech campus in Wishnick Hall 210.

# Why?

The purpose of this workshop will be
- to further the education of workshop participants in the field of computational biochemistry;
- to promote international scientific collaboration, especially between the United States, Vietnam, and Romania
<!-- - to provide a rich cultural exchange experience that will foster this collaboration. -->
