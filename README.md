# Trajectory Mesh Learning (tML)
In Jiang et al 2021 preprint (https://www.researchsquare.com/article/rs-146357/private/preview) we describe a reinforcement learning agent that learns to scale foraging trajectories. The learning rules the govern such a model are described in the text and are related in part to a previously described model (http://www.nature.com/nature/journal/v533/n7603/full/nature17639.html). Example Matlab scripts that implement these rules in the context of the tasks described in the paper. These data correspond to part of the data in Figure 2 of the preprint. Additional code for analysis of larger scale model runs, grid search optimizations, more detailed scripts for generating figures, and behavioral data will be placed here upon publication.

These scripts were modified to have minimal dependencies and best illustrate how the learning process described below works. They were tested on Matlab 2018b. A few remaining conveience functions are included in the Dependencies folder and will need to be added to user path to run scripts.

## Model description from the preprint

The tML agent model is based around the idea that an agent can learn to scale the parameters of a structured representation of foraging trajectories. For a central place forager we might demand a trajectory that forms a closed out and back loop which begins at a ‘home’ location, transits through an extrema and returns home. The goal of the learning agent is to update the heading and amplitude of this trajectory so that it reliably intercepts a target location according to the specific rules of the environment. For example, interception may need to occur at the trajectory extrema or perhaps anywhere along the trajectory or perhaps for some fixed duration. In the specific cases for this study we consider interception either at the extrema (STF task) or for a fixed duration (NTF task) that correspond to the practical requirements of our real-time behavior analysis used in the experimental task designs. We note that similar results to those reported have been obtained with a range of different simulated environments.

Returning to the notion of a structured, closed-loop trajectory, we consider the problem as a control signal that determines behavior at each time step. First consider a locomoting animal. At each time step we assume it is governed by a heading angle and an instantaneous speed. Under such a model a closed-loop trajectory will be produced by a smooth rotation in heading angle (a linear function from -pi to pi). For a fixed speed this would produce a rotation about a circle. However, to produce the observed, roughly elliptical paths speed is inhomogeneous and reaches maxima along specific heading angles - outward runs (pi/2) and return runs (-pi/2). Given the expected bell-shaped distribution of speeds that minimize jerk along a trajectory this can be modeled as a sequence of gaussian speed profiles. These dynamics for heading and speed can be generated by an artificial neural network, but for simplicity we have used simple generative functions perturbed by noise. A schematic of the model architecture can be found in Fig. 2 and Fig. 7.

Θ(t) = L(-π,π) + ω[i] + ε                                                                                                        (1)

S(t) =G(tau,σ) × a[i] + ε                                                                                                        (2)  

Where, L is a linear mapping across the range {-π,π} spanning time t, ω[i] is the heading offset sampled from a distribution of mean Ω[i] and constant variance for each trial i. G is either a single (NTF) or double (STF) peaked gaussian function with offsets tau={tau1} or {tau1, tau2} and width 𝝈, scaled by a gain a[i] sampled from a distribution with a mean A[i] and constant variance for each trial i. ε is a normally distributed and smoothed noise term matched to observed variability in observed behavioral trajectories. 

Given a model of structured trajectories defined by continuous speed (S(t)) and heading (𝛳(t)), the learning problem for an agent is to learn to scale the amplitude (A(i)) and orientation (Ω(i)) offsets of trajectories trial by trial in order to reliably intercept target locations. Behavioral data indicated bidirectional and rapid learning for changes in the scaling of movement trajectories, thus we used a modified version of a learning rule (MeSH) previously described to account for rapid, bidirectional movement scaling 32,35.

A[i+1] = A[i] +𝛼 (a[i]-A[i]) 𝜐[i] - β(a[i]-A[0]))                                                               (3)

Ω[i+1] = Ω[i] +𝛼 (ω[i]-Ω[i]) 𝜐[i] - β(ω[i]-Ω[0]))                                                               (4)

where , iis the index of the ith trial, 𝜐[i]is a smoothed estimate of the local reward rate, a[i]is the magnitude of the speed on the current trial as sampled from a normal distribution centered on A[i] with rate parameters 𝛼 and β. Learning rate parameters and the standard deviation of the distribution 𝒩(A,σ) were explored using grid search optimization. The equivalent learning rule is also expressed for Ω(i) in equation 4.
