
# Final project

# Introduction

The final project is an extended version of Kauffman's $NK$-model (see, e.g., [[1]](#1) and [[2]](#2)). The model
captures a genome consisting of N (binary) genes, and a genome can be represented in the form of an $N$- dimensional bitstring. Every gene contributes the 
$0 \leq c_i \leq 1$ to the genome's fitness. However, the fitness contribution might not depend on a single gene, but there are at most $K = N-1$ 
interdependencies with other genes. This means that the fitness contribution $c_i$ of gene $n_i$
is not only affected by its binary value but also by the values that the interdependent (up to $K$) other genes take.
The fitness contributions are computed according to:

$$ c_i = f(n_i, n_{i_1}, \dots, n_{i_K}) \sim U(0, 1),$$

where $n_i \in \\{0, 1\\}$ indicates gene $i \in \\{1, \dots, N\\}$ and $\\{i_1, \dots, i_K\\} \subseteq \\{1, \dots, i-1, i+1, \dots, N\\}$.
The genomes fitness is the mean of the fitness contributions:

$$ C = \frac{1}{N} \sum_{i=1}^{N} c_i .$$

The fitness contributions are captured in so-called $NK$-fitness landscapes. Their complexity can be tuned via the parameter $K$. The higher (lower) $K$, more more (less) rugged is the resulting fitness landscape, and the more difficult (easier) it is to find the global maximum in the landscape. The fitness landscape can best be imagined as a hypercube that becomes more complex with increses in $K$. 

The genome's evolution is captured as follows: An agent moves in the fitness lanscape by performing hill climbing-based search for genomes with a higher fitness. However, to avoid long jumps, the space in which the agent can move is restricted to the neighbourhood of the current genome (the current bit-string), whereby the neighourhood is defined as the Hamming Distance (and is often set to 1). This means, the neighbourhood are all genomes in which one bit is flipped compared to the current genome. 

The proposed model is an extension of Kauffman's $NK$ model: In the original model, it is assumed that evolution happens without noise. In the implementation included here, I add an optional error term to the process of evolution, meaning that the fitness of genomes if uncertain during hill climbing (but, of course, certain once the genome is selected). 

# Distinctiveness and Complexity

Kauffman's $NK$-model is a well-established model that is applied in research in a multiplicity of fields (also beyond the field of biology). Initializing the landscape is a computationally complex problem that is solved following an object-oriented approach. Also, the project does not only implment but also extend the original $NK$-model. 

# How to run the program

How to run your application.

# Features

What can the program do


# Project structure

What’s contained in each file you created.



# References
<a id="1">[1]</a> 
Kauffman, S.A. and Weinberger, E.D., 1989. The NK model of rugged fitness landscapes and its application to maturation of the immune response. Journal of Theoretical Biology, 141(2), pp.211-245.

<a id="2">[2]</a> 
Weinberger, E.D., 1991. Local properties of Kauffman’s N-k model: A tunably rugged energy landscape. Physical review A, 44(10), p.6399.





Any other additional information the staff should know about your project.
