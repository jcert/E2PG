# E2PG

This repository complements the paper titled "Epidemic Population Games for Policy Design: Two Populations with Viral Reservoir Case Study", and accommodates Julia files for simulating an epidemic compartmental model coupled to a population game.

## Concept

We extend to two populations a recently proposed system theoretic framework for studying an epidemic influenced by the strategic behavior of a single population's agents. Our framework couples the well-known susceptible-infected-susceptible (SIS) epidemic model with a population game that captures the strategic interactions among the agents of two large populations. This framework can also be employed to study a situation where a population of _ nonstrategic_ agents _such as animals_  serves as a disease reservoir. Equipped with the framework, we investigate the problem of designing a suitable control policy that assigns dynamic payoffs to incentivize the agents to adopt costlier and more effective strategies subject to a long-term budget constraints. We formulate a nonconvex constrained optimization program for minimizing the disease transmission rate at an endemic equilibrium, and explain how to obtain an approximate solution efficiently. Then, we propose a dynamic payoff mechanism and prove the convergence to an (approximately) optimal social state of the population game which minimizes the basic reproduction number, using a Lyapunov function.

[preprint](https://arxiv.org/abs/)


## Requirements
- Julia 1.6.7
- (*optional*) jupyter-notebook 6.2

## How to use
Following the [guide on environments](https://pkgdocs.julialang.org/v1.2/environments/), you can open Julia in a terminal, press `]` to access the package manager, type `activate .` and then `instantiate`. 
After installing all the required software you can press backspace to exit the package manager, now you should have all the required libraries to run the code. To run the code either use Jupyter notebook for the interactive plot or open Julia and then type `include("two_pop_epg.jl")` to run the main simulation (that will generate Figures 1 and 2).


## two_pop_epg.jl
File to generate the plots used for Figures 1 and 2


## two\_pop\_epg.ipybn
Jupyter Notebook file with an interactive plot, in it you can change the parameters of the game and epidemic disease.

This simulation follows the structure of example 1, using the the dynamic payoff to reach an equilibrium. The third and fourth figures exemplify Remark 2 of the paper, showing that we can perturb the optimal social state found through the optimization problem with minimal effect on the spectral radius and find another target social state for which the theorem holds.

Here are a few screenshots of the interactive plot being used 


![Screenshot from 2023-03-28 02-36-54](https://user-images.githubusercontent.com/13306869/228150023-35b32b11-b17c-4f66-90e7-0137b4f77f95.png)
![Screenshot from 2023-03-28 02-38-20](https://user-images.githubusercontent.com/13306869/228150026-7a5b593f-4e82-4544-8c0a-91223c70d72c.png)
![Screenshot from 2023-03-28 02-39-57](https://user-images.githubusercontent.com/13306869/228150034-42e19bc6-00cf-42b4-87dc-b1a95648f4b0.png)
![Screenshot from 2023-03-28 02-40-02](https://user-images.githubusercontent.com/13306869/228150035-255d000c-30b1-4bf7-b108-c09ccddf4314.png)






