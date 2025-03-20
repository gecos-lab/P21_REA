Fracture areal intensity (P21) is defined as the ratio between the total sum of fracture trace lengths and the sampling area. It serves as the stopping criterion in 2D Discrete Fracture Networks (DFNs), meaning that the stochastic generation of fractures ceases once the selected intensity is saturated.

Due to the heterogeneous distribution of fractures in natural outcrops, the characterization of P21 is closely tied to the concept of Representative Elementary Area (REA).

Using the digitized fracture traces (polyline shapefile) and the interpretation boundary (polygon shapefile), this code calculates the P21 parameter and determines the REA through the following steps:

- Several hexagonal grids are distributed across the outcrop surface, with side lengths ranging from 1 m to 26 m (depend on the outcrop size), increasing in 1-meter steps. Only whole hexagons are considered.

- Tukey's boxplot method (1977) is used to analyze the distribution of P21 values at different scan area sizes.

- REA can be statistically defined by checking if there are not significant differences between the mean and standard deviation of a P21 sample collected at a certain scan area size, the previous and successive P21 sample. This is defined with a qualitative approach based on on the difference between the interquartile range (IQR) of two subsequential P21 samples.

- To account for “far out” data, that are not included in the IQR, we added the range between the whiskers calculated as the difference between the upper whisker length (Q3 +1.5IQR) and the lower whisker length (Q1 - 1.5IQR).

- REA correspond to a plateau in the parameter value variation.

Dependencies: Fracability (https://github.com/gecos-lab/FracAbility)
