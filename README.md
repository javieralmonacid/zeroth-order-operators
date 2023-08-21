# Characterization of singular flows of zeroth-order pseudo-differential operators via elliptic eigenfunctions: a numerical study

## Article abstract

The propagation of internal gravity waves in stratified media, such as those found in ocean basins and lakes, leads to the development of geometrical patterns called "attractors". These structures accumulate much of the wave energy and make the fluid flow highly singular. In more analytical terms, the cause of this phenomenon has been attributed to the presence of a continuous spectrum in some nonlocal zeroth-order pseudo-differential operators. In this work, we analyze the generation of these attractors from a numerical analysis perspective. First, we propose a high-order pseudo-spectral method to solve the evolution problem (whose long-term behaviour is known to be not square-integrable). Then, we use similar tools to discretize the corresponding eigenvalue problem. Since the eigenvalues are embedded in a continuous spectrum, we compute them using viscous approximations. Finally, we explore the effect that the embedded eigenmodes have on the long-term evolution of the system.

## How to use this repository?

- Clone the repository: ```git clone git@github.com:javieralmonacid/zeroth-order-operators.git```
- Open Matlab and add the ```tools``` folder to the path: ```addpath tools```
- Run the scripts (most of the codes can be run independently from each other). Each script produces a .mat files to avoid doing the main computation every time the script is executed.

## How to cite this repository?

Just cite the main article:

> J. A. Almonacid and N. Nigam. *Characterization of singular flows of zeroth-order pseudo-differential operators via elliptic eigenfunctions: a numerical study.* Journal of Computational and Applied Mathematics (2023), doi: [https://doi.org/10.1016/j.cam.2023.115510](https://doi.org/10.1016/j.cam.2023.115510).

This reference will be updated once the article has been published in its final form.

## External tools used

- ColorBrewer: Attractive and Distinctive Colormaps ([GitHub](https://github.com/DrosteEffect/BrewerMap), [Mathworks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps)),
- ```sptoeplitz``` [(Mathworks File Exchange)](https://www.mathworks.com/matlabcentral/fileexchange/13353-sparse-toeplitz-matrix-construction),
- ```textprogressbar``` [(Mathworks File Exchange)](https://www.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar).
