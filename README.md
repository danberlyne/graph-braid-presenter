# About this tool
This tool combines computational methods for graph braid groups developed in the [graph-braid-splitter](https://github.com/danberlyne/graph-braid-splitter) and [graph-morse](https://github.com/tmaciazek/graph-morse) packages by Daniel Berlyne and Tomasz Maciazek (the authors of this package). 

Given a graph braid group, this tool will compute a presentation for that group. It does this by first attempting to split the group as a free product, where the factors of this splitting are expressed as either free groups or simpler graph braid groups. Next, the tool expresses each of these free factors as a direct product if it is directly decomposable. The tool will then attempt to further split these direct factors, alternating between free products and direct products until a splitting with freely indecomposable factors is acquired. The tool will then compute the presentations of any graph braid group factors that appear in the splitting, finally combining the presentations into a presentation for the group as a whole. 

# Instructions
1. Enter the data for your graph braid group in the `gbg_data.txt` file. Detailed instructions are given in the file, but briefly, this should include the following.
- The number of particles/strands.
- The adjacency matrix for your graph. Note that [this handy online graph calculator](https://www.mas.ncl.ac.uk/graph-curvature/)[^1] can easily produce an adjacency matrix from a hand-drawn graph, which can then be copy-pasted into `gbg_data.txt`.
- An initial configuration of the particles, if required.

[^1]: *The Graph Curvature Calculator and the curvatures of cubic graphs, Experimental Mathematics, 2019
(arXiv:1712.03033 [math.CO])*

2. Run `presenter.py` and follow the on-screen instructions.

# How to read output data in `presentation.txt` and `splitting.txt`
If a presentation is found, the generators and relators will be displayed in the terminal. The relators are expressed as lists of 2-tuples consisting of a generator and a power. Multiplying these powers of generators gives the relator. This presentation data is also saved to `presentation.txt`.

If the original graph braid group is split before the presentation is computed, this splitting data will be saved to `splitting.txt`. Factors in the splitting will either be: 
- free groups, displayed in the form Z or F_k;
- other known groups, displayed in the form given in `known_gbgs.txt`;
- graph braid groups, displayed in the form B_n(Gamma_i) or RB_n(Gamma_i).

The following further information about B_n(Gamma_i) and RB_n(Gamma_i) can be found in `splitting.txt`.
- The adjacency matrix of Gamma_i.
- The presentation of B_n(Gamma_i) or RB_n(Gamma_i) used in the computation of the presentation of the original graph braid group.

Tip: Click 'load' in the bottom-right corner of [this handy online graph calculator](https://www.mas.ncl.ac.uk/graph-curvature/)[^1] and paste the adjacency matrix to quickly obtain a picture of the graph.

# Known graph braid groups
Below is a list of known graph braid groups that are non-cyclic and directly and freely indecomposable, with citations in the literature. These are used to aid in computations. The author would appreciate any further contributions to this list. The data for these is included in `known_gbgs.txt` (note: in order to make these graphs canonical, any vertices of degree 2 must be removed before adding the data to this list).
1. B_2(K_5) = pi_1(N_6), where K_5 is the complete graph on 5 vertices and N_6 is the closed non-orientable surface of genus 6.
- Example 3.18 of ['Graph of groups decompositions of graph braid groups'](https://arxiv.org/pdf/2209.03860.pdf).
- Example 5.1 of ['Configuration spaces and braid groups of graphs'](https://www.proquest.com/docview/304583880).
2. B_2(K_{3,3}) = pi_1(N_4), where K_{3,3} is the complete bipartite graph on two sets of 3 vertices and N_4 is the closed non-orientable surface of genus 4.
- Example 3.19 of ['Graph of groups decompositions of graph braid groups'](https://arxiv.org/pdf/2209.03860.pdf).
- Example 5.2 of ['Configuration spaces and braid groups of graphs'](https://www.proquest.com/docview/304583880).
3. B_3(\theta_4) = pi_1(S_3), where \theta_n is a generalised theta graph and S_3 is the closed orientable surface of genus 3.
- Example 4.3 of ['Characteristics of graph braid groups'](https://arxiv.org/pdf/1101.2648.pdf).