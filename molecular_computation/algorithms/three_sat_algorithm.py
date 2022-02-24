import math

from molecules.dna_sequence import DNASequence
from molecules.dna_molecule import DNAMolecule
from enzymes.ligase import Ligase


class ThreeSATAlgorithm:

    __literal_solution_size = 20
    __edge_sequence_split_size = 10
    __literals_repr_dict = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    __clause_size = 3

    def __init__(self, clauses):
        """
        Initialize the Three SAT Algorithm with formula to solve.

        :param clauses: The formula to solve, in form of:<br>
        [C_1, C_2, ..., C_n]<br>
        Where C_i is the clause if the formula but in form of:<br>
        [literal_num_1, literal_num_2, literal_num_3]<br>
        Where literal_num_i is the number of the literal as part of the whole literals count.<br>
        E.g [[1,2,3], [-1, 2, 3]] represents (A∨B∨C)∧(A'∨B∨C) as the number of literals are 3.
        """
        if not (
            type(clauses) == list and
            all([type(clause) == list and len(clause) == self.__clause_size for clause in clauses])
        ):
            raise ValueError(f'Given bad clause format. Each clause should contain 3 literals.')

        self.__clauses = clauses
        self.__num_literals = max(abs(literal_index) for clause in self.__clauses for literal_index in clause)

        if self.__num_literals > len(self.__literals_repr_dict):
            raise ValueError(f'The number of literals is above expectation.\n'
                             f'Should be at most {len(self.__literals_repr_dict)}')

    def solve(self):
        """
        Solves the Three SAT Algorithm for the given formula.
        """
        self.__print_formula_representation()

        print('Generating solution sequences...')
        edges, literals, negated_literals, \
            vertices, starting_vertex, ending_vertex, \
            comp_literals, comp_negated_literals, comp_vertices, \
            comp_starting_vertex, comp_ending_vertex = self.__generate_solution_sequences()
        self.__print_representations_info(
            edges,
            literals,
            negated_literals,
            vertices,
            starting_vertex,
            ending_vertex,
            comp_literals,
            comp_negated_literals,
            comp_vertices,
            comp_starting_vertex,
            comp_ending_vertex
        )

        print('Mixing solution values in test tube with ligase...')
        generated_paths = self.__mix_components(
            edges,
            comp_literals,
            comp_negated_literals,
            comp_vertices,
            comp_starting_vertex,
            comp_ending_vertex
        )

        solution_identification_map, len_path_id = self.__create_solution_identification_map(generated_paths)
        self.__print_paths_by_solution_identification_map(solution_identification_map, len_path_id)

        print('Filtering the satisfying paths...')
        satisfying_paths = self.__filter_satisfying_paths(generated_paths, comp_literals, comp_negated_literals)
        self.__print_paths_by_solution_identification_map(solution_identification_map, len_path_id, satisfying_paths)

        results = \
            "There is combination of values that satisfies the given formula." \
            if len(satisfying_paths) > 0 \
            else "There is no combination of values that satisfies the given formula."
        print(f'Results: {results}')

    def __generate_solution_sequences(self):
        """
        Generates the whole representations for the algorithm.

        :return: List containing the representations for:<br>
            [<br>
            edges - The edges of the possible path solution graph template,<br>
            literals - Representation of each literal,<br>
            negated_literals - Representation of each negated literal,<br>
            vertices - Representation of the vertices in the solution graph template,<br>
            starting_vertex - Representation of the starting vertex in the solution graph template,<br>
            ending_vertex - Representation of the ending vertex in the solution graph template,<br>
            comp_literals - Representation of the complement literals,<br>
            comp_negated_literals - Representation of the complement negated literals,<br>
            comp_vertices - Representation of the complement vertices in the solution graph template,<br>
            comp_starting_vertex - Representation of the complement starting vertex in the solution graph template,<br>
            comp_ending_vertex - Representation of the complement ending vertex in the solution graph template<br>
            ]
        """
        literals = []
        negated_literals = []
        existing_literals = {}
        comp_literals = []
        comp_negated_literals = []

        starting_vertex = None
        ending_vertex = None
        comp_starting_vertex = None
        comp_ending_vertex = None

        vertices = []
        comp_vertices = []

        edges = []

        vertices_len = self.__num_literals - 1
        literals_created = 0

        # Create representations for the literals and their complements
        while literals_created != self.__num_literals:
            dna_seq = DNASequence.create_random_sequence(self.__literal_solution_size)
            raw_seq = ''.join(dna_seq.bases)

            if existing_literals.get(raw_seq) is None:
                literals.append(dna_seq)
                comp_dna_seq = dna_seq.get_complement()
                comp_literals.append(comp_dna_seq)

                existing_literals[raw_seq] = True
                existing_literals[''.join(comp_dna_seq.bases)] = True

                literals_created += 1

        # Create representations for the negated literals and their complements
        literals_created = 0
        while literals_created != self.__num_literals:
            dna_seq = DNASequence.create_random_sequence(self.__literal_solution_size)
            raw_seq = ''.join(dna_seq.bases)

            if existing_literals.get(raw_seq) is None:
                negated_literals.append(dna_seq)
                comp_dna_seq = dna_seq.get_complement()
                comp_negated_literals.append(comp_dna_seq)

                existing_literals[raw_seq] = True
                existing_literals[''.join(comp_dna_seq.bases)] = True

                literals_created += 1

        # Create representations for the vertices and their complements
        literals_created = 0
        while literals_created != self.__num_literals - 1:
            dna_seq = DNASequence.create_random_sequence(self.__literal_solution_size)
            raw_seq = ''.join(dna_seq.bases)

            if existing_literals.get(raw_seq) is None:
                vertices.append(dna_seq)
                comp_dna_seq = dna_seq.get_complement()
                comp_vertices.append(comp_dna_seq)

                existing_literals[raw_seq] = True
                existing_literals[''.join(comp_dna_seq.bases)] = True

                literals_created += 1

        # Create representation for the starting vertex and its complement
        found = False
        while not found:
            dna_seq = DNASequence.create_random_sequence(self.__literal_solution_size)
            comp_dna_seq = dna_seq.get_complement()
            raw_seq = ''.join(dna_seq.bases)
            raw_comp_dna_seq = ''.join(comp_dna_seq.bases)

            if existing_literals.get(raw_seq) is None and existing_literals.get(raw_comp_dna_seq) is None:
                starting_vertex = dna_seq
                comp_starting_vertex = comp_dna_seq
                found = True

        # Create representation for the ending vertex and its complement
        found = False
        while not found:
            dna_seq = DNASequence.create_random_sequence(self.__literal_solution_size)
            comp_dna_seq = dna_seq.get_complement()
            raw_seq = ''.join(dna_seq.bases)
            raw_comp_dna_seq = ''.join(comp_dna_seq.bases)

            if existing_literals.get(raw_seq) is None and existing_literals.get(raw_comp_dna_seq) is None:
                ending_vertex = dna_seq
                comp_ending_vertex = comp_dna_seq
                found = True

        # Create the edges of the template solution graph, starting from the starting vertex
        current_literal_vertex = 0
        edges.append(
            DNASequence(
                starting_vertex.bases +
                literals[current_literal_vertex].bases[:self.__edge_sequence_split_size]
            )
        )
        edges.append(
            DNASequence(
                starting_vertex.bases +
                negated_literals[current_literal_vertex].bases[:self.__edge_sequence_split_size]
            )
        )

        # Create the edges of the template solution graph,
        # using each time pair of literal and negated literal which connects to vertex
        # and after it vertex which connect to literal and negated literal, leaving the last
        # literal and negated literal to connect to the ending vertex
        for index in range(vertices_len):
            # Edges from literals to vertex
            edges.append(
                DNASequence(
                    literals[current_literal_vertex].bases[self.__edge_sequence_split_size:] +
                    vertices[index].bases[:self.__edge_sequence_split_size]
                )
            )
            edges.append(
                DNASequence(
                    negated_literals[current_literal_vertex].bases[self.__edge_sequence_split_size:] +
                    vertices[index].bases[:self.__edge_sequence_split_size]
                )
            )

            current_literal_vertex += 1

            # Edges from vertex to literals
            edges.append(
                DNASequence(
                    vertices[index].bases[self.__edge_sequence_split_size:] +
                    literals[current_literal_vertex].bases[:self.__edge_sequence_split_size]
                )
            )
            edges.append(
                DNASequence(
                    vertices[index].bases[self.__edge_sequence_split_size:] +
                    negated_literals[current_literal_vertex].bases[:self.__edge_sequence_split_size]
                )
            )

        # At last, connect the last literal and negated literal with the ending vertex
        edges.append(
            DNASequence(
                literals[current_literal_vertex].bases[self.__edge_sequence_split_size:] +
                ending_vertex.bases
            )
        )
        edges.append(
            DNASequence(
                negated_literals[current_literal_vertex].bases[self.__edge_sequence_split_size:] +
                ending_vertex.bases
            )
        )

        # Return the whole representations
        return [
            edges,
            literals,
            negated_literals,
            vertices,
            starting_vertex,
            ending_vertex,
            comp_literals,
            comp_negated_literals,
            comp_vertices,
            comp_starting_vertex,
            comp_ending_vertex
        ]

    @staticmethod
    def __mix_components(
        edges,
        comp_literals,
        comp_negated_literals,
        comp_vertices,
        comp_starting_vertex,
        comp_ending_vertex
    ):
        """
        Mix the whole components to create all the possible paths in the template solution graph as
        mentioned, with the enzyme ligase.

        :param edges: The edges of the template solution graph.
        :param comp_literals: Complement representation of the literals.
        :param comp_negated_literals: Complement representation of the negated literals.
        :param comp_vertices: Complement representation of the vertices of the template solution graph.
        :param comp_starting_vertex: Complement representation of the starting vertex of the template solution graph.
        :param comp_ending_vertex: Complement representation of the ending vertex of the template solution graph.
        :return: List of all possible paths as DNA molecules
        """
        # Number of possibilities.
        possibilities = 2 ** len(comp_literals)

        # The current path selection as binary representation, where
        # 0 means negated literal and 1 regular literal.
        # The 0b1 remains as placeholder for representing binary number as the size
        # of the number of literals.
        # E.g. if we have 3 literals, the binary number will be 0b1000 which saves 3 placeholders
        # after the 0b1 for all the possibilities of selections.
        current_selection = 0b1 << len(comp_literals)
        selection_extractor = lambda selection: bin(selection).replace('0b1', '')

        # Will contain all possible paths as DNA molecules.
        paths = []

        # For each possibility, create new path and add it into the paths array
        for possibility in range(possibilities):

            # Extract the current selection from the binary representation
            actual_selection = selection_extractor(current_selection)

            # Select the starting of the first edge (by the selection of the first literal)
            selected_edge_literal = int(actual_selection[0])

            # Building the top DNA molecule strand of the curr path
            curr_path = edges[selected_edge_literal]

            # Continue to the next edge
            selected_edge_literal += 2
            curr_path = Ligase.combine_molecules(curr_path, edges[selected_edge_literal])

            # Continue to the next selection of literal as the 2 edges of the previous one selected,
            # need to determine the next edge by the position of the current edge - if it in odd index
            # means the negated literal is selected and we need to increase the edge index by 1.
            # Otherwise, it is in even index - means the literal is selected and need to increase
            # it by 2.
            selected_edge_literal += 2 if selected_edge_literal % 2 == 0 else 1

            # Building the bottom DNA molecule strand of the curr path (as comp. for the top strand)
            comp_curr_path = comp_starting_vertex
            comp_curr_path = Ligase.combine_molecules(
                comp_curr_path,
                comp_negated_literals[0] if actual_selection[0] == '1' else comp_literals[0]
            )

            # Add the other vertices and literals based on the actual selection
            for index in range(0, len(comp_literals) - 1):

                # Select the correct edge as the result of the current selection literal
                selected_edge_literal += int(actual_selection[index + 1])
                curr_path = Ligase.combine_molecules(
                    curr_path,
                    edges[selected_edge_literal]
                )

                # Continue to the next edge of the current selection
                selected_edge_literal += 2
                curr_path = Ligase.combine_molecules(
                    curr_path,
                    edges[selected_edge_literal]
                )

                # Continue to the next selection of literal as the 2 edges of the previous one selected,
                # need to determine the next edge by the position of the current edge - if it in odd index
                # means the negated literal is selected and we need to increase the edge index by 1.
                # Otherwise, it is in even index - means the literal is selected and need to increase
                # it by 2.
                selected_edge_literal += 2 if selected_edge_literal % 2 == 0 else 1

                # Select and combine the path using the current vertex and the actual selection of the literal
                comp_curr_path = Ligase.combine_molecules(
                    comp_curr_path,
                    comp_vertices[index]
                )
                comp_curr_path = Ligase.combine_molecules(
                    comp_curr_path,
                    comp_negated_literals[index + 1] if actual_selection[index + 1] == '1' else comp_literals[index + 1]
                )

            # Combine the end vertex to mark the current path as complete
            comp_curr_path = Ligase.combine_molecules(comp_curr_path, comp_ending_vertex)

            # Add the current generated path combination as DNA molecule
            paths.append(
                DNAMolecule(
                    curr_path,
                    comp_curr_path
                )
            )

            # Continue to the next possible selection
            current_selection += 1

        return paths

    def __filter_satisfying_paths(self, paths, comp_literals, comp_negated_literals):
        """
        Filter the satisfying paths from all possible paths given.<br>
        This done by iterating each clause, and filter only the paths which satisfies the clause (by finding the
        part of a literal in the DNA molecule at the second strand, which satisfies the clause - as can be
        performed with magnetic beads).<br>
        By doing so, after all the clauses are satisfied, only the paths which satisfy the whole formula left over.

        :param paths: All possible paths of the template solution graph.
        :param comp_literals: Complement representation of the literals.
        :param comp_negated_literals: Complement representation of the negated literals.
        :return: List of paths that satisfy the whole formula
        """
        def contains_literal_bases(literal_bases, path_bases):
            for i in range(len(path_bases) - len(literal_bases) + 1):
                for j in range(len(literal_bases)):
                    if path_bases[i + j] != literal_bases[j]:
                        break
                else:
                    return True
            return False

        current_satisfying_test_tube = [*paths]

        # Filtering each clause with the satisfying paths according to each literal
        for clause in self.__clauses:

            # Previous satisfying test tube for tracking the elimination of paths each clause
            prev_satisfying_test_tube = current_satisfying_test_tube[:]

            # Current satisfying test tube for that clause (And the whole clause one after one)
            current_satisfying_test_tube = []

            # Remaining paths tube after filter each literal satisfaction
            remaining_paths_tube = prev_satisfying_test_tube[:]

            for curr_literal in clause:

                # The index of the literal in both comp_literals and comp_negated_literals lists
                literal_index = abs(curr_literal) - 1

                # If the literal has positive value - should search when the literal is 1
                # Otherwise the literal has negative value - should search when the literal is 0
                filtered_literal = \
                    comp_literals[literal_index] \
                    if curr_literal > 0 \
                    else comp_negated_literals[literal_index]

                # The new remaining paths tube as the results of the previous paths filtered
                new_remaining_paths_tube = []

                for path in remaining_paths_tube:

                    # Check if the filtered literal is part of the path -
                    # Which means it included in the path and therefore, satisfies the clause.
                    contains_literal = contains_literal_bases(filtered_literal.bases, path.second_strand.bases)

                    # Path contain the filtered literal - means we can save the path as satisfying
                    # Otherwise - throw it to the new remaining paths tube
                    if contains_literal:
                        current_satisfying_test_tube.append(path)
                    else:
                        new_remaining_paths_tube.append(path)

                # Setting the new remaining paths tube as the remaining paths tube for the next round
                remaining_paths_tube = new_remaining_paths_tube

        return current_satisfying_test_tube

    def __print_formula_representation(self):
        """
        Prints the formula in friendly format (In mathematical format).
        """
        formula = ''
        and_sign = '∧'
        or_sign = '∨'
        negated_sign = "'"

        for clause in self.__clauses:
            formula += '('

            for literal in clause:
                current_literal = self.__literals_repr_dict[abs(literal) - 1]
                formula += current_literal

                # If the literal is negated, need to add the negated sign
                if literal < 0:
                    formula += negated_sign

                formula += or_sign

            # Remove the last or sign
            formula = formula[:len(formula) - 1]

            formula += f'){and_sign}'

        # Remove the last and sign
        formula = formula[:len(formula) - 1]

        print('Given formula to solve:')
        print(f'{formula}\n')

    def __print_representations_info(
        self,
        edges,
        literals,
        negated_literals,
        vertices,
        starting_vertex,
        ending_vertex,
        comp_literals,
        comp_negated_literals,
        comp_vertices,
        comp_starting_vertex,
        comp_ending_vertex
    ):
        """
        Prints the representation of the components for the algorithm.

        :param edges - The edges of the possible path solution graph template.
        :param literals - Representation of each literal.
        :param negated_literals - Representation of each negated literal.
        :param vertices - Representation of the vertices in the solution graph template.
        :param starting_vertex - Representation of the starting vertex in the solution graph template.
        :param ending_vertex - Representation of the ending vertex in the solution graph template.
        :param comp_literals - Representation of the complement literals.
        :param comp_negated_literals - Representation of the complement negated literals.
        :param comp_vertices - Representation of the complement vertices in the solution graph template.
        :param comp_starting_vertex - Representation of the complement starting vertex in the solution graph template.
        :param comp_ending_vertex - Representation of the complement ending vertex in the solution graph template.
        """
        print(f'\nLiterals Representations:')
        for index, literal in enumerate(literals):
            print(f'Literal {self.__literals_repr_dict[index]}:            {literal}')
            print(f'Literal {self.__literals_repr_dict[index]} Complement: {comp_literals[index]}')

        print(f'\nNegated - Literals Representations:')
        for index, negated_literal in enumerate(negated_literals):
            print(f"Negated-Literal {self.__literals_repr_dict[index]}':            {negated_literal}")
            print(f"Negated-Literal {self.__literals_repr_dict[index]}' Complement: {comp_negated_literals[index]}")

        print(f'\nVertices Representations:')
        for index, vertex in enumerate(vertices):
            print(f'Vertex {index}:            {vertex}')
            print(f'Vertex {index} Complement: {comp_vertices[index]}')

        print(f'\nStarting Vertex:            {starting_vertex}')
        print(f'Complement Starting Vertex: {comp_starting_vertex}')
        print(f'Ending Vertex:              {ending_vertex}')
        print(f'Complement Ending Vertex:   {comp_ending_vertex}')

        print(f'\nEdges ({len(edges)}) Representations:')
        for index in range(0, len(edges), 4):
            literal_repr = self.__literals_repr_dict[index // 4]
            print(f"From vertex a{index // 4}  --->  {literal_repr}:  {edges[index]}")
            print(f"From vertex a{index // 4}  --->  {literal_repr}': {edges[index + 1]}")
            print(f'From literal {literal_repr}  --->  a{(index // 4) + 1}: {edges[index + 2]}')
            print(f"From literal {literal_repr}' --->  a{(index // 4) + 1}: {edges[index + 3]}")
            print()

    @staticmethod
    def __print_paths_by_solution_identification_map(solution_identification_map, len_path_id, partial_paths=None):
        """
        Prints the paths by the solution identification map.

        :param solution_identification_map: Dictionary which maps each path to representational string (which identifies
        the solution as binary number which represent the values for each literal (0 for False, 1 for True))
        and backwards.
        :param len_path_id: The length of the representational string of the solution identification map.
        :param partial_paths: If given, showing only the partial paths as part of the solution identification map,
        using the backward mapping to show the representational string.
        """

        # If partial paths presented, means we need to print only those paths and their selection representation
        if partial_paths is not None:
            if len(partial_paths) > 0:
                print('\nPaths:')
                for partial_path in partial_paths:
                    print(f'Path representing selection {solution_identification_map[partial_path.__repr__()]}: '
                          f'\n{partial_path.__repr__()}\n')

        # Otherwise, printing all the paths as usual by their path id
        else:
            print('\nPaths:')
            for path_id in solution_identification_map.keys():

                # Filter only the representational string format instead printing duplicates values
                if len(path_id) == len_path_id:
                    print(f'Path representing selection {path_id}: \n{solution_identification_map[path_id]}\n')

    @staticmethod
    def __create_solution_identification_map(paths):
        """
        Creating a solution identification map which maps each path to representational string and the opposite.<br>
        The representational string identifies the solution as binary number which represent the values set
        for each literal (0 for False, 1 for True)).<br>
        For example, given A,B,C literals 101 represents A=True, B=False, C=True.

        :param paths: The paths to map into solution identification map.
        :return: Dictionary of solution identification map (as described above).
        """
        path_selection = (0b1 << (int(math.log2(len(paths))) + 1)) - 1
        solution_identification_map = {}

        for path in paths:
            solution_identification_map[bin(path_selection).replace("0b1", "")] = path
            solution_identification_map[path.__repr__()] = bin(path_selection).replace("0b1", "")
            path_selection -= 1

        return solution_identification_map, int(math.log2(len(paths)))


if __name__ == '__main__':
    three_sat = ThreeSATAlgorithm([[1, 2, 3], [1, 2, -3], [-1, -2, 4], [1, -2, -4]])
    three_sat.solve()
