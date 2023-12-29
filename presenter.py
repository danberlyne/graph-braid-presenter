#!/usr/bin/env python3
# presenter.py - Computes presentations of graph braid groups

from splitting_tools.graph import Graph
from splitting_tools.graph_braid_group import GraphBraidGroup, is_same
from splitting_tools.splitter import enter_data_manually, get_data_from_file, VertexException, start_exit_sequence, stringify_factors, combine_strings
from morse_tools.morse_utils import graph_braid_group

import sys
import itertools
import copy
import numpy as np

# Takes as input an object of class `Graph`.
# Returns a graph in the correct format for use in `graph_braid_group()` in `morse_utils.py`.
def reformat_graph_to_morse(graph, num_particles):
    # First identify a root for the spanning tree.
    root = get_root(graph)
    # Next, pick a path of length `num_particles` starting at `root`.
    rooted_path = get_rooted_path(graph, num_particles, root)
    # Generate a spanning tree, starting at the end of the rooted path.
    spanning_tree = graph.get_spanning_tree(rooted_path, rooted_path[0][-1])
    # Edit adjacency matrix by adding '-1' for edges not in the spanning tree, reordering rows and columns to match the order of the vertices in `spanning_tree[0]`.
    # Reordering in this way is important to make sure the order agrees with [Farley-Sabalka's vertex ordering](https://arxiv.org/abs/math/0410539).
    spanning_matrix = [[graph.adj_matrix[spanning_tree[0][i]][spanning_tree[0][j]] for j in range(len(graph.adj_matrix))] for i in range(len(graph.adj_matrix))]
    for edge in [e for e in graph.edges if e not in spanning_tree[1]]:
        spanning_matrix[spanning_tree[0].index(edge[0])][spanning_tree[0].index(edge[1])] = -1
        spanning_matrix[spanning_tree[0].index(edge[1])][spanning_tree[0].index(edge[0])] = -1
    print(spanning_matrix)
    # Then convert matrix to a NumPy array.
    return np.array(spanning_matrix)

# Returns a vertex of the Graph object `graph` to be used as a root for a spanning tree.
def get_root(graph):
    # If all vertices have valence 2, then the graph is a cycle and any vertex can be chosen as the root.
    if graph.essential_vertices == []:
        return 0
    # Otherwise, pick an essential vertex, prioritising vertices of valence 1.
    else:
        for v in graph.essential_vertices:
            if graph.get_degree(v) == 1:
                return v
        return graph.essential_vertices[0]

# Returns a path in `graph` of length `num_paticles`, starting at vertex `root`.
def get_rooted_path(graph, num_particles, root):
    rooted_path = ([root], [])
    while len(rooted_path[0]) < num_particles:
        end_vertex = rooted_path[0][-1]
        for edge in [e for e in graph.edges if end_vertex in e]: 
            if edge not in rooted_path[1]:
                rooted_path[0].append(edge[edge.index(end_vertex) - 1])
                rooted_path[1].append(edge)
                break
    return rooted_path

def print_presentation(gens, rels):
    print("Generators:")
    for i in range(len(gens)):
        print(f'g{i} = {gens[i]}')
    print('')   
    print("Relators:")
    for rel in rels:
        print(rel)

def print_presentation_without_gen_info(gens, rels):
    print("Generators:")
    for i in range(len(gens)):
        print(gens[i])
    print('') 
    print("Relators:")
    for rel in rels:
        print(rel)

def save_presentation(gens, rels):
    presentation_string = 'Generators:\n'
    for i in range(len(gens)):
        presentation_string += f'g{i} = {gens[i]}\n' 
    presentation_string += '\nRelators:\n'
    for rel in rels:
        presentation_string += f'{rel}\n'
    file = open('presentation.txt', 'w')
    file.write(presentation_string)
    file.close()

def save_presentation_without_gen_info(gens, rels):
    presentation_string = 'Generators:\n'
    for i in range(len(gens)):
        presentation_string += f'{gens[i]}\n'
    presentation_string += '\nRelators:\n'
    for rel in rels:
        presentation_string += f'{rel}\n'
    file = open('presentation.txt', 'w')
    file.write(presentation_string)
    file.close()

# Computes presentations of all GraphBraidGroup factors in the splitting and 
# returns a version of `splitting` with the GraphBraidGroup factors replaced with (gens, rels) 2-tuples.
# Note, `previous_gbgs` should be initially input as an empty dictionary. Its keys will be graph braid groups we have already encountered 
# and its values will be 2-tuples consisting of the presentation and the output of `is_reduced()`.
def compute_presentations(splitting, previous_gbgs):
    for i, factor in enumerate(splitting):
        if type(factor) == GraphBraidGroup:
            if factor.is_reduced():
                essential_graph = factor.graph.make_essential()
                is_new_gbg = True
                # Loop over all B_n(Gamma) we have already encountered.
                for previous_gbg in [gbg for gbg in previous_gbgs if previous_gbgs[gbg][1]]:
                    # If we already encountered this graph braid group, use the presentation we already computed.
                    if is_same(essential_graph.adj_matrix, previous_gbg.graph.make_essential().adj_matrix) and factor.num_particles == previous_gbg.num_particles:
                        splitting[i] = previous_gbgs[previous_gbg][0]
                        is_new_gbg = False
                        break
                if is_new_gbg:
                    arr = reformat_graph_to_morse(factor.graph, factor.num_particles)
                    gens, rels = graph_braid_group(arr, factor.num_particles)
                    splitting[i] = (gens, rels)
                    previous_gbgs.update({factor: ((gens, rels), True)})
            else:
                is_new_gbg = True
                # Loop over all RB_n(Gamma) we have already encountered.
                for previous_gbg in [gbg for gbg in previous_gbgs if not previous_gbgs[gbg][1]]:
                    if is_same(factor.adj_matrix, previous_gbg.graph.adj_matrix) and factor.num_particles == previous_gbg.num_particles:
                        splitting[i] = previous_gbgs[previous_gbg][0]
                        is_new_gbg = False
                        break
                if is_new_gbg:
                    arr = reformat_graph_to_morse(factor.graph, factor.num_particles)
                    gens, rels = graph_braid_group(arr, factor.num_particles)
                    splitting[i] = (gens, rels)
                    previous_gbgs.update({factor: ((gens, rels), False)})
        elif isinstance(factor, list):
            compute_presentations(factor, previous_gbgs)
    return splitting

# For each direct factor of the graph braid group, combines presentations into a single list of generators and a single list of relators.
# Edits `presentations` in place so that it is a list of (generators, relators) 2-tuples.
def combine_presentations(presentations, is_direct = True):
    for i, factor in enumerate(presentations):
        # If a factor is a list of non-lists, then replace the factor in `presentations` with its combined presentation.
        if isinstance(factor, list):
            for subfactor in factor:
                # If we find a list in factor's subfactors, then feed factor back into `combine_strings`.
                if isinstance(subfactor, list):
                    combine_presentations(factor, not is_direct)
                    break
            # If `presentations` is a direct splitting, then `factor` is a free splitting.
            if is_direct:
                # Collect together F_m terms and Z terms.
                free_rank = improved_sum(int(subfactor[2:]) for subfactor in factor if isinstance(subfactor, str) and subfactor.startswith('F_')) + len([subfactor for subfactor in factor if subfactor == 'Z'])
                # Replace F_m and Z terms with their presentations.
                factor = [([j for j in range(k)], []) for k in [free_rank] if free_rank >= 1] + [subfactor for subfactor in factor if not (isinstance(subfactor, str) and subfactor.startswith('F_')) and not subfactor == 'Z']
                if len(factor) > 1:
                    # Concatenate the lists of generators, reindexing so we don't have repeats.
                    gens = [f'g{factor[j][0].index(gen) + improved_sum(len(factor[k][0]) for k in range(j))}' for j in range(len(factor)) for gen in factor[j][0]]
                    # Concatenate the relators, replacing any generator substrings with the reindexed versions.
                    rels = [[(rel[k][0][0] + str(int(rel[k][0][1:]) + improved_sum(len(factor[m][0]) for m in range(j))), rel[k][1]) for k in range(len(rel))] for j in range(len(factor)) for rel in factor[j][1]]
                    presentations[i] = (gens, rels)
                else:
                    presentations[i] = factor[0]
            else:
                # Replace F_m and Z terms with their presentations
                factor = [([j for j in range(int(subfactor[2:]))], []) for subfactor in factor if isinstance(subfactor, str) and subfactor.startswith('F_')] + [([0],[]) for subfactor in factor if subfactor == 'Z'] + [subfactor for subfactor in factor if not (isinstance(subfactor, str) and subfactor.startswith('F_')) and not subfactor == 'Z']
                if len(factor) > 1:
                    # Concatenate the lists of generators, reindexing so we don't have repeats.
                    for j in range(len(factor)):
                        factor[j] = ([f'g{factor[j][0].index(gen) + improved_sum(len(factor[k][0]) for k in range(j))}' for gen in factor[j][0]], factor[j][1])
                    gens = [gen for j in range(len(factor)) for gen in factor[j][0]]
                    # Concatenate the relators, replacing any generator substrings with the reindexed versions, then add extra commutation relators.
                    rels = [[(rel[k][0][0] + str(int(rel[k][0][1:]) + improved_sum(len(factor[m][0]) for m in range(j))), rel[k][1]) for k in range(len(rel))] for j in range(len(factor)) for rel in factor[j][1]]
                    commutators = [[(a, 1), (b, 1), (a, -1), (b, -1)] for (j, k) in itertools.combinations(range(len(factor)), 2) for a in factor[j][0] for b in factor[k][0]]
                    rels += commutators
                    presentations[i] = (gens, rels)
                else:
                    presentations[i] = factor[0]

# Has same functionality as `sum`, except returns 0 when the generator is empty.
def improved_sum(generator):
    return sum([0] + [gen for gen in generator])

#############
# Main code #
#############

def main():
    print('Checking gbg_data.txt...')
    file = open('gbg_data.txt')
    file_as_list = file.readlines()
    file.close()
    for i, line in enumerate(file_as_list):
        file_as_list[i] = line.strip()
    non_empty_lines = [line for line in file_as_list if line != '']
    dotted_line = non_empty_lines.index('-' * 50)
    gbg_data = non_empty_lines[dotted_line + 1:]

    # If data has not been entered in `gbg_data.txt`, prompt user to enter it manually.
    if len(gbg_data) == 3:
        (num_particles, adj_matrix, initial_config) = enter_data_manually()
    # Otherwise, get data from `gbg_data.txt` and verify it is formatted correctly.
    else: 
        (num_particles, adj_matrix, initial_config) = get_data_from_file(gbg_data)

    gbg = GraphBraidGroup(Graph(adj_matrix), num_particles, initial_config)
    if not gbg.is_reduced():
        print('''WARNING: In order to perform computations for B_n(\Gamma), the graph \Gamma must satisfy the following conditions:
1. All cycles must have length at least n+1.
2. All paths between vertices of degree not equal to 2 must have length at least n-1.
At least one of these conditions is not satisfied by your graph.
If you choose to continue, any results obtained will only be true for the reduced braid group RB_n(\Gamma).
Do you wish to continue? (y/n)''')
        while True:
            response = input().lower()
            if response == 'n':
                while True:
                    exit_sequence = input('Please amend your data then run the script again. Press [Enter] to exit.')
                    if not exit_sequence:
                        sys.exit()
            elif response == 'y':
                for comp in gbg.num_initial_particles_per_component:
                    if len(comp[0]) < gbg.num_initial_particles_per_component[comp]:
                        raise VertexException
                break

    print('Verified successfully.\nSearching for free splittings...')

    if gbg.is_trivial():
        if gbg.is_reduced():
            print(f'B_{num_particles} = 1')
            start_exit_sequence()
        else:
            print(f'RB_{num_particles} = 1')
            start_exit_sequence()

    # Returns a nested list, where the odd levels of the list correspond to direct factors and the even levels correspond to free factors.
    # We avoid splittings with non-trivial graph of groups factors.
    splitting = gbg.get_splitting_without_gogs()

    print('Search complete.')

    # Makes fresh copies of `presentation.txt` and `splitting.txt`, the files containing presentation information and splitting information.
    file = open('presentation.txt', 'w')
    file.write('')
    file.close()
    file = open('splitting.txt', 'w')
    file.write('')
    file.close()

    # If no splittings are found, compute the presentation of the graph braid group as-is.
    if len(splitting) == 1:
        if not isinstance(splitting[0], str):
            if not isinstance(splitting[0], list):
                print('No splittings found. Computing presentation...')
                spanning_array = reformat_graph_to_morse(gbg.graph, gbg.num_particles)
                gens, rels = graph_braid_group(spanning_array, gbg.num_particles)
                print('Presentation found.')
                print_presentation(gens, rels)
                print('The relators above are expressed as lists of 2-tuples consisting of a generator and a power. Multiplying these powers of generators gives the relator.')
                save_presentation(gens, rels)
                print('Saved to presentation.txt.')
                start_exit_sequence()
            elif len(splitting[0]) == 1 and not isinstance(splitting[0][0], str):
                print('No splittings found. Computing presentation...')
                spanning_array = reformat_graph_to_morse(gbg.graph, gbg.num_particles)
                gens, rels = graph_braid_group(spanning_array, gbg.num_particles)
                print('Presentation found.')
                print_presentation(gens, rels)
                print('The relators above are expressed as lists of 2-tuples consisting of a generator and a power. Multiplying these powers of generators gives the relator.')
                save_presentation(gens, rels)
                print('Saved to presentation.txt.')
                start_exit_sequence()

    # Otherwise, compute presentations of factors and create a copy of `splitting` containing these presentations instead of the factors.
    print('Splitting found. Saved to splitting.txt.\nComputing presentations of factors...')
    presentations = copy.deepcopy(splitting)
    presentations = compute_presentations(presentations, {})

    # Turns all factors in the splitting into strings and adds detailed data to `splitting.txt` where appropriate.
    # This includes saving the presentations of each GraphBraidGroup factor to `splitting.txt`.
    stringify_factors(splitting, 1, 1, True, presentations)

    combine_strings(splitting)
    final_string = ' x '.join(splitting)

    # Adds the splitting to the beginning of `splitting.txt`.
    with open('splitting.txt','r') as contents:
        save = contents.read()
    if gbg.is_reduced():
        with open('splitting.txt','w') as contents:
            contents.write(f'B_{num_particles} = ' + final_string + '\n\n')
        with open('splitting.txt','a') as contents:
            contents.write(save)
    else:
        with open('splitting.txt','w') as contents:
            contents.write(f'RB_{num_particles} = ' + final_string + '\n\n')
        with open('splitting.txt','a') as contents:
            contents.write(save)

    # For each direct factor of the graph braid group, combine the presentations appearing in that factor.
    combine_presentations(presentations)
    # Finally, combine the presentations of the direct factors to give a presentation of the whole group.
    # First replace F_m and Z terms with their presentations
    presentations = [([f'g{j}' for j in range(int(factor[2:]))], []) for factor in presentations if isinstance(factor, str) and factor.startswith('F_')] + [([f'g{0}'],[]) for factor in presentations if factor == 'Z'] + [factor for factor in presentations if not (isinstance(factor, str) and factor.startswith('F_')) and not factor == 'Z']
    if len(presentations) > 1:
        # Concatenate the lists of generators, reindexing so we don't have repeats.
        for j in range(len(presentations)):
            presentations[j] = ([f'g{presentations[j][0].index(gen) + improved_sum(len(presentations[k][0]) for k in range(j))}' for gen in presentations[j][0]], presentations[j][1])
        final_gens = [gen for j in range(len(presentations)) for gen in presentations[j][0]]
        # Concatenate the relators, replacing any generator substrings with the reindexed versions, then add extra commutation relators.
        final_rels = [[(rel[k][0][0] + str(int(rel[k][0][1:]) + improved_sum(len(presentations[m][0]) for m in range(j))), rel[k][1]) for k in range(len(rel))] for j in range(len(presentations)) for rel in presentations[j][1]]
        commutators = [[(a, 1), (b, 1), (a, -1), (b, -1)] for (j, k) in itertools.combinations(range(len(presentations)), 2) for a in presentations[j][0] for b in presentations[k][0]]
        final_rels += commutators
    else:
        (final_gens, final_rels) = presentations[0]
    print('Presentation found.')
    print_presentation_without_gen_info(final_gens, final_rels)
    print('The relators above are expressed as lists of 2-tuples consisting of a generator and a power. Multiplying these powers of generators gives the relator.')
    save_presentation_without_gen_info(final_gens, final_rels)
    print('Saved to presentation.txt.')
    start_exit_sequence()

if __name__ == '__main__':
    main()