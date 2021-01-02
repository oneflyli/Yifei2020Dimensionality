Code_Figure5.m: 1.Produce the averaged density along y=50. 
                2.Produce the evolution of total population desity with the 2-dimensional initial condition.


Code_Figure6.m: 1.Produce the averaged column density. 
                2.Produce the evolution of total population desity with the 1-dimensional initial condition.

Code_Figure6.m: 1. Produce snapshots of one realisation. Although the results are generated with the 0-dimensional initial condition, you can produce snapshots with other initial conditions by changing parameters 'shape' and 'distribution'. 
                2. Produce the evolution of total population desity with the 0-dimensional initial condition.

PhaseDiagram_Discrete: Produce the phase diagram of P/M and C_0 in the discrete model. As varying P/M and C(0) with small steps is time-consuming, this code only provides few values of these two parameters.

PhaseDiagram_Continuous: Produce the phase diagram of P/M and C_0 in the discrete model. Again, as varying P/M and C(0) with small steps is time-consuming, this code only provides few values of these two parameters.

Please run the 5 scripts above to obtain results of discrete/continuum models.

There are three external functions:

iteration.m:Algorithms of movement mechanism, growth mechanism, growth crowding function and periodic boundary conditions for discrete simulations.

LineApproach_odefun_1D.m: Method of lines for solving a 1-dimensional RDE

LineApproach_odefun.m: Method of lines for solving a 2-dimensional RDE