#### Noiseless simulation of the QAOA ansatz circuit for the MaxCut problem

See `demo.ipynb` for an example (nbviewer [link](https://nbviewer.org/github/PratikSathe/QAOA_approxratio_landscape/blob/main/demo.ipynb)).

This repository also contains `angles_regular_graphs.json`, a file containing the values of fixed angles from the [fixed angle conjecture](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.104.052419%5D), which shows that at these angles, the QAOA ansatz results in a high approximation ratio. This data file was copied over from the (fixed angles repository)[[https://github.com/danlkv/fixed-angle-QAOA]](https://github.com/danlkv/fixed-angle-QAOA%5D).

Some parts of this code have been adapted from the MaxCut implementation in the [QED-C benchmarking framework](https://github.com/PratikSathe/QC-App-Oriented-Benchmarks), which the author contributed to. As of Nov 29 2022, this code is in the `add_maxcut` branch of the repository.



The (unweighted) graph for which the maxcut problem is to be solved should be specified in the form of a list of tuples, with each tuple denoting an edge of the graph. Currently, the code in this repository works only for unweighted maxcut problems.

Once the `methods.qaoa_anstaz` is instantiated, the graph can be visualized using `ansatz.draw_graph()`. The graph used in the demo is:



![](.\figures\graph.png)



The ansatz circuit can be printed by calling `ansatz.draw_ansatz()`

The empirically obtained distribution of cut sizes can be plotted using `ansatz.plot_cutsize_dist()`. The plot that is generated will also show the distribution of cut sizes for a random uniform sampling over cuts.

![cutsize](.\figures\Distribution%20of%20Cut%20sizes.png)

It is straightforward to create a wrapper function to compute the approximation ratio for any given $\beta$ and $\gamma$ values. For example, the approximation ratio landscape for the graph above, obtained using the statevector simulation looks like below:

![landscapear](.\figures\Sweep-Size=4-Rounds=1.png)

Here, the star represents the location of the fixed angle from the fixed angle conjecture. The square marks the location of the grid point for which the highest value of the approximation ratio was obtained.



As an aside, we also plot the landscape of the best cut ratio with 5000 shots. (Here, best cut ratio equals the largest cut size from 5000 shots divided by the known maxcut value.) It is clear from the landscape that it is entirely featureless, implying that it is not a good objective function for iterative QAOA.

![bestcutlandscape](.\figures\Sweep-Size=4-Rounds=1_bestcutratio.png)