## (Under Construction!)

# Random trees of polytopes for approximate optimal control of piecewise affine systems
## Sadra Sadraddini and Russ Tedrake
### MIT CSAIL, July 2018

### Abstract
Piecewise affine (PWA) systems are widely used to model highly nonlinear behaviors such as contact dynamics in robot locomotion and manipulation. Existing control techniques for PWA systems have computational drawbacks, both in offline design and online implementation. 
In this paper, we introduce a method to obtain feedback control policies and a corresponding  set of admissible initial conditions for discrete-time PWA systems such that all the closed-loop trajectories reach a goal polytope, while a cost function is optimized. 
The idea is conceptually similar to LQR-trees \cite{tedrake2010lqr}, which consists of 3 steps: (1) open-loop trajectory optimization, (2) feedback control for computation of "funnels" of states around trajectories, and (3) repeating (1) and (2) in a way that the funnels are grown backward from the goal in a tree fashion and fill the state-space as much as possible. We show PWA dynamics can be exploited to combine step (1) and (2) into a single step that is tackled using mixed-integer convex programming, which makes the method more suitable for dealing with hard constraints. Illustrative examples on contact-based dynamics are presented. 

### Paper
The full version (corrections made) is available [here](https://github.com/sadraddini/PWA-Control/blob/master/paper.pdf)

### Dependencies:
* Python 3.6 or later
* Gurobi Version 7.0 or later

### Folder descriptions
* Main: source code of tools
* Examples: 
    * The Bouncing ball 
    * Inverted pendulum with Wall, Example 1 in this [paper](http://groups.csail.mit.edu/robotics-center/public_papers/Marcucci17.pdf)
    * 4-Mode System from Example 2 in in this [paper](https://www.researchgate.net/profile/Michal_Kvasnica/publication/4143171_Computation_of_invariant_sets_for_piecewise_affine_discrete_time_systems_subject_to_bounded_disturbances/links/54d0b5930cf298d65668244c/Computation-of-invariant-sets-for-piecewise-affine-discrete-time-systems-subject-to-bounded-disturbances.pdf)

### How to use examples:
(under development)

### How to use guide:
(under development)