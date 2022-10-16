
# Final project in course CS50P

# Introduction

The final project is an extended version of Kauffman's $NK$-model (see, e.g., [[1]](#1) and [[2]](#2)). The model
captures a genome consisting of N (binary) genes, and a genome can be represented in the form of an $N$- dimensional bitstring. Every gene contributes the 
$0 \leq c_i \leq 1$ to the genome's fitness. However, the fitness contribution might not depend on a single gene, but there are at most $K = N-1$ 
interdependencies with other genes. This means that the fitness contribution $c_i$ of gene $n_i$
is not only affected by its binary value but also by the values that the interdependent (up to $K$) other genes take.
The fitness contributions are computed according to:

$$ c_i = f(n_k, n_{i_1}, \dots, n_{i_K}) \sim U(0, 1),$$

where $n_i \in \\{0, 1\\}$ indicates gene $i \in \\{1, \dots, N\\}$ and $\\{i_1, \dots, i_K\\} \subseteq \\{1, \dots, i-1, i+1, \dots, N\\}$.
The genomes fitness is the mean of the fitness contributions:

$$ C = \frac{1}{N} \sum_{i=1}^{N} c_i .$$



# Distinctiveness and Complexity
Why you believe your project satisfies the distinctiveness and complexity requirements, mentioned above.

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
