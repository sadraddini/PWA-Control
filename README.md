# Random Trees of Polytopes for Approximate Optimal Control of Piecewise Affine Systems
## Sadra Sadraddini and Russ Tedrake
### [Robot Locomotion Group](http://groups.csail.mit.edu/locomotion/), MIT CSAIL

![Bouncing Ball](https://github.com/sadraddini/PWA-Control/raw/master/Examples/Bouncing_ball/figures/ball_iterations.gif)
![Inverted Pendulum With Wall](https://github.com/sadraddini/PWA-Control/raw/master/Examples/Inv_pendulum_wall/figures/inv_pendulum_wall_iterations.gif)

### Abstract
Piecewise affine (PWA) systems are widely used to model highly nonlinear behaviors such as contact dynamics in robot locomotion and manipulation. Existing control techniques for PWA systems have computational drawbacks, both in offline design and online implementation. 
In this paper, we introduce a method to obtain feedback control policies and a corresponding  set of admissible initial conditions for discrete-time PWA systems such that all the closed-loop trajectories reach a goal polytope, while a cost function is optimized. 
The idea is conceptually similar to LQR-trees [Tedrake et. al. (2010)](https://groups.csail.mit.edu/robotics-center/public_papers/Tedrake10.pdf), which consists of 3 steps: (1) open-loop trajectory optimization, (2) feedback control for computation of "funnels" of states around trajectories, and (3) repeating (1) and (2) in a way that the funnels are grown backward from the goal in a tree fashion and fill the state-space as much as possible. We show PWA dynamics can be exploited to combine step (1) and (2) into a single step that is tackled using mixed-integer convex programming, which makes the method more suitable for dealing with hard constraints. Illustrative examples on contact-based dynamics are presented. 

### Paper
The full version (corrections made) is available [here](https://github.com/sadraddini/PWA-Control/blob/master/paper.pdf)

### Dependencies:
* [Gurobi](http://www.gurobi.com/) Version 7.0 or later. Free Academic use.

### Folder Descriptions
* Main: source code of methods
* Convex_Hull (*new*): methods related to disjuctive programming for faster tree extention 
* Examples: 
    * Bouncing ball 
    * Inverted pendulum with wall, Example 1 in [paper](http://groups.csail.mit.edu/robotics-center/public_papers/Marcucci17.pdf)
    * Inverted pendulum with two walls - one on each side
    * Stablizing 4-cell PWA system: Example 2 in [paper](https://www.researchgate.net/profile/Michal_Kvasnica/publication/4143171_Computation_of_invariant_sets_for_piecewise_affine_discrete_time_systems_subject_to_bounded_disturbances/links/54d0b5930cf298d65668244c/Computation-of-invariant-sets-for-piecewise-affine-discrete-time-systems-subject-to-bounded-disturbances.pdf)
    * Hybrid stabilization of Planar Pushing (*new*) ![Pushing](https://raw.githubusercontent.com/sadraddini/PWA-Control/master/Examples/pushing_box/figures/pushing.gif)

### How to use guide:
The user may use the following to formulate a PWA control problem and obtain a controller. The following guide is the general picture, and it does not include minor details. The reader is encouraged to check the examples. 

* Step I: define the PWA system. Specify the modes, the valuations for affine dynamics in each mode, and the half-space representation of the polytope corresponding to it. 
So far, we have only considered cases where PWA cells are constructed in state space, not joint state-control space. Also, the user is asked to provide the bounding box of each PWA cell, which is used for sampling.

* Step II: Define a goal polytope. The problem is designing a control policy to get into the goal. The user may also define a cost function. The default is time optimality. 

* Step III: run the polytope tree algorithm. The tree grows incrementally, and it may take a long time such that polytopes cover a large of portion of the state space.

* Step IV: two controllers are obtained:
    * The first controller is simple. It is based on matrix multiplications and is essentially an affine feedback in each polytope. It does not work properly for states out of the tree. 
    * The second controller solves a small convex program to keep the system within the tree, or close to the tree, while decreasing the value function. This controller can handle states outside of the tree, but does not provide any guarantee that the state gets into the goal- unless the state falls into the tree.

### Visualization
Current version only supports visualization for 2D problems. For higher dimensions, projections to 2D are performed. See examples. 

### Contact us
If you have any questions regarding this work, or you are willing to contribute, please contact [Sadra Sadraddini](sadra@mit.edu) 

Last updated on Sep 17, 2018. 