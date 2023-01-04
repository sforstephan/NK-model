
# Introduction

The project is an extended version of Kauffman's $NK$ model (see, e.g., [[1]](#1) and [[2]](#2)). The model
captures a genome consisting of N (binary) genes, and a genome can be represented in the form of an $N$- dimensional bitstring. Every gene contributes the 
$0 \leq c_i \leq 1$ to the genome's fitness. However, the fitness contribution might not depend on a single gene, but there are at most $K = N-1$ 
interdependencies with other genes. This means that the fitness contribution $c_i$ of gene $n_i$
is not only affected by its binary value but also by the values that the interdependent (up to $K$) other genes take.
The fitness contributions are computed according to:

$$ c_i = f(n_i, n_{i_1}, \dots, n_{i_K}) \sim U(0, 1),$$

where $n_i \in \\{0, 1\\}$ indicates gene $i \in \\{1, \dots, N\\}$ and $\\{i_1, \dots, i_K\\} \subseteq \\{1, \dots, i-1, i+1, \dots, N\\}$.
The genome's fitness is the mean of the fitness contributions:

$$ C = \frac{1}{N} \sum_{i=1}^{N} c_i .$$

The fitness contributions are captured in so-called $NK$-fitness landscapes. Their complexity can be tuned via the parameter $K$. The higher (lower) $K$, more more (less) rugged is the resulting fitness landscape, and the more difficult (easier) it is to find the global maximum in the landscape. The fitness landscape can best be imagined as a hypercube that becomes more complex with increases in $K$. 

The genome's evolution is captured as follows: An agent moves in the fitness landscape by performing hill climbing-based search for genomes with a higher fitness. However, to avoid long jumps, the space in which the agent can move is restricted to the neighbourhood of the current genome (the current bit-string), whereby the neighbourhood is defined as the Hamming Distance. This means, the neighbourhood are all genomes in which one bit is flipped compared to the current genome. 

The proposed model is an extension of Kauffman's $NK$ model: In the original model, it is assumed that evolution happens without noise. In the implementation included here, I add an optional normally distributed error term to the process of evolution, meaning that the fitness of genomes if uncertain during hill climbing (but, of course, certain once the genome is selected). 

# Distinctiveness and Complexity

Kauffman's $NK$ model is a well-established model that is applied in research in a multiplicity of fields (also beyond the field of biology). Initializing the landscape is a computationally complex problem that is solved following an object-oriented approach. Also, the project does not only implement but also extend the original $NK$ model. 


# Features

The model simulated the process of evolution for an $N$-dimensional genome with $K$ interdependencies between genes. 
The genome's fitness is captured in a $NK$-fitness landscape, and the process of evolution is implemented as an agent that
moves on the landscape. 

For explanations on how to run the model and the parameters see [here](#run): The pattern of dependency can either be selected from a pre-defined dependency matrix (`main`, `block`, `ring`) or randomly 
generated using parameters `n` and `k`. The evolutionary process is simulated `repeat` times for `period` time steps.  
There is uncertainty in the evolutionary process, meaning that the agent might make a normally distributed error with mean `mean` and standard deviation `std` when estimating the fitness during the process of evolution. The model stores and returns simulated data, 
statistics, and a plot containing the mean fitness and confidence intervals with confidence level `confidence`. The computation of confidence intervals is contingent on the number of repetitions (i.e., observations). For less than 30 observations per timee step, the computation of the confidence intervals
relies on the t-Distribution, and on the Normal Distribution otherwise. 

The model stores the following files:

- `parameters.json`: A JSON file containing the entered parameters.
- `fitness.json`: A JSON file containing the genome's fitness for all time steps and all repetitions.
- `statistics.json`: A JSON file containing the mean fitness for each time steps (over all repetitions), and the lower and upper boundaries for the confidence interval. The confidence-level is defined by the parameter `confidence` 
- `fitness.jpg`: A plot that includes the genome's mean fitness for all time steps and confidence intervals. The confidence-level is defined by the parameter `confidence` 

# <a id="run"></a>How to run the program

### Start the program via the command line 
Start the program via the command line (`python project.py`). 


### Optional parameters
The following arguments are **optional**:
  - `-n`: Number of genes in the genome, integer, default value = `5`
  - `-k`: Number of interdependencies between genes, integer (0 <= k <= N-1), default value = `2`
  - `-matrix`: Pre-defined dependency matrix (see [below](#matrices)), type string. Allowed keys are: `main`, `block`, `ring`, `random`, default value = `random`. If random, a dependency matrix with parameter `n` and `k` is randomly generated, otherwise one of the pre-defined matrices is used. 
  - `-time`: Time steps in one simulation round, i.e., time steps the genome can evolve, integer, default value = `100`
  - `-repeat`: Number of simulation rounds, i.e., number of times the entire process of evolution is observed, integer, default value = `8`
  - `-mean`: Mean of the prediction error during evolution, float, efault value = `0.0`
  - `-std`: Standard deviation of the error during evolution, float, efault value = `0.0`
  - `-confidence`: Confidence interval used for the plot, float, default value = `0.9`

To avoid long jumps (and chaotic behavior), the neighbourhood during evolution is fixed at a Hamming Distance of $1$. 

### <a id="matrices"></a>Pre-defined dependency matrices

The following pre-defined dependency matrices can be used as optional parameter (or a random dependency matrix can be generated, see above). Performance contributions are in rows, genes are in columns. `1` indicates that a dependency between genes exists, `0` indicates that no dependency exists.  

- Key `main`: A unity matrix (only 1 at the main diagonal) that indicates that no dependencies between genes exist. The resulting fitness landscape only has one peak (that is relatively easy to find).
```
[
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
]
```

- Key `block`: Block-like interdependencies along the main diagonal. 

```
[
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
]
```

- Key `ring`: Ring-like interdependencies, extension of the block structure. The first block is interdependent with the second block, the second block is interdependent with the third block, etc. 

```
[
    [1, 1, 0, 0, 0, 0, 0, 0, 1, 1],
    [1, 1, 0, 0, 0, 0, 0, 0, 1, 1],
    [1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
    [0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 1, 1, 1, 0, 0],
    [0, 0, 0, 0, 1, 1, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
]
```

# Required packages

- numpy
- argparse
- SciPy
- matplotlib
- json
- pytest


# Project structure

The model is organized in three files.
- `NK.py` contains the class definitions for `NK`(the parend class) and `Landscape` and `Agent` (the two child classes)
- `project.py` includes the main simulation loop and functions to create and check random dependency matrices and data analysis
- `matrices.py` includes the pre-defined dependency matrices

```bash
.
|-- NK.py
|   |-- class NK (parent class)
|   |   |-- properties: n, k 
|   |   |-- get()
|   |   |-- get_n()
|   |   |-- get_k()
|   |   |-- convert_number_to_bin_nparray()
|   |-- class Landscape (child class)
|   |   |-- properties: global_max_position, contributions, lookuptables, alleles, dependencymap
|   |   |-- get()
|   |   |-- check_dependencymap()
|   |   |-- get_required_fitness_contributions()
|   |   |-- check_genome()
|   |   |-- get_fitness_gene()
|   |   |-- get_fitness_genome()
|   |   |-- get_global_max()
|   |-- class Agent (child class)
|   |   |-- properties: position, error_mean, error_std
|   |   |-- get()
|   |   |-- flip_bit()
|   |   |-- get_random_position()
|   |   |-- get_alternative_position()
|   |   |-- evolve()
|   |   |-- update_position()
|-- project.py
|   |-- main()
|   |-- get_random_matrix()
|   |-- check_matrix()
|   |-- compute_statistics()
|-- test_project.py
|   |-- main()
|   |-- test_get_random_matrix()
|   |-- test_check_matrix()
|   |-- test_compute_statistics()
|-- matrices.py
|-- readme.md
```




# References
<a id="1">[1]</a> 
Kauffman, S.A. and Weinberger, E.D., 1989. The NK model of rugged fitness landscapes and its application to maturation of the immune response. Journal of Theoretical Biology, 141(2), pp.211-245.

<a id="2">[2]</a> 
Weinberger, E.D., 1991. Local properties of Kauffman?s N-k model: A tunably rugged energy landscape. Physical review A, 44(10), p.6399.
