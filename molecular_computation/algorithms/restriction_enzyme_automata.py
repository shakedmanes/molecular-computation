from molecules.dna_sequence import DNASequence
from molecules.dna_molecule import DNAMolecule
from enzymes.fok1 import Fok1
from enzymes.ligase import Ligase


class RestrictionEnzymeAutomata:

    # Alphabet of the automaton (^ represent the start, $ represent the terminator)
    __alphabet = {
        '^': DNASequence('GGATG'),
        'a': DNASequence('CTGGCT'),
        'b': DNASequence('CGCAGC'),
        '$': DNASequence('TGTCGC')
    }

    # States of the automaton
    __states = {
        ''.join(__alphabet['a'].bases[2:]): 'q0',
        ''.join(__alphabet['b'].bases[2:]): 'q0',
        ''.join(__alphabet['$'].bases[2:]): 'q0',

        ''.join(__alphabet['a'].bases[:4]): 'q1',
        ''.join(__alphabet['b'].bases[:4]): 'q1',
        ''.join(__alphabet['$'].bases[:4]): 'q1'
    }

    # Accepted states of the automaton
    __accepting_states = ['q0']

    # Represent the transition rules in the automaton
    __transition_rules = [

        # q0 --(a)--> q0
        {
            'from': 'q0',
            'to': 'q0',
            'molecule':
                DNAMolecule(
                    DNASequence('GGATGTACSSSS'),
                    DNASequence('CCTACATGCCGA')
                )
        },

        # q0 --(b)--> q1
        {
            'from': 'q0',
            'to': 'q1',
            'molecule':
                DNAMolecule(
                    DNASequence('GGATGACGACSSSS'),
                    DNASequence('CCTACTGCTGGTCG')
                )
        },

        # q1 --(a)--> q1
        {
            'from': 'q1',
            'to': 'q1',
            'molecule':
                DNAMolecule(
                    DNASequence('GGATGACGSSSS'),
                    DNASequence('CCTACTGCGACC')
                )
        },

        # q1 --(b)--> q0
        {
            'from': 'q1',
            'to': 'q0',
            'molecule':
                DNAMolecule(
                    DNASequence('GGATGGSSSS'),
                    DNASequence('CCTACCGCGT')
                )
        }
    ]

    def __init__(self, automaton_input):
        """
        Initializes the restriction enzyme automaton by input given to the automaton.

        :param automaton_input: String representing the input the automaton should process.
        """
        self.__original_input = automaton_input
        self.__input = self.__build_input_representation('^' + automaton_input + '$')

    def process_automaton(self):
        """
        Process the input given to the automaton and determine if the automaton accept / reject the input.
        """
        print(f'Given original input:\n{self.__original_input}')
        print(f'\nInput in automaton format:\n{self.__input}')

        # Copy the automaton representation input to work on copy rather than the original representation input.
        copy_input = DNAMolecule(self.__input.first_strand, self.__input.second_strand)

        # Performing the cleavage using Fok1 enzyme
        print('\nApplying Fok1 enzyme to cleavage input:')
        restriction_output = Fok1.cleavage(copy_input)

        print('Results of cleavage:')
        print(restriction_output)
        print(f'\nSticky part: {restriction_output.sticky_part}')

        # Marking the raw original input index and current and next states.
        # The current state is the start state q0, and the next state is resolved from
        # the sticky part of the cleavage results.
        original_input_index = 0
        current_state = 'q0'
        next_state = self.__states[restriction_output.sticky_part.__repr__()]
        print(f'\nState Transition: {current_state}'
              f' --{self.__original_input[original_input_index]}--> {next_state}')

        # Trying to fine compatible transition rule from the current state reached which found at the
        # sticky part of the cleavage output results.
        compatible_transition_rule_molecule = self.__find_compatible_transition_rule(restriction_output)

        # Also, the sticky part of the cleavage results belongs to specific letter representation and therefore
        # connects to compatible transition rule which knows the letter received and applying the transition
        # to the appropriate next state.
        if compatible_transition_rule_molecule is not None:
            print('\nAppropriate transition rule molecule:')
            print(compatible_transition_rule_molecule['molecule'])
            print(f'\nSticky part: {compatible_transition_rule_molecule["molecule"].sticky_part}')

            print('\nCombining cleavage output and transition rule molecules:')
            copy_input = DNAMolecule(
                Ligase.combine_molecules(
                    compatible_transition_rule_molecule['molecule'].first_strand,
                    restriction_output.first_strand
                ),
                Ligase.combine_molecules(
                    compatible_transition_rule_molecule['molecule'].second_strand,
                    restriction_output.second_strand
                )
            )
            print(copy_input)

        # If no transition rule found, means we reach end of input.
        else:
            print('Reached end of input.')
            copy_input = restriction_output
            restriction_output = None

        while restriction_output is not None:
            print('\nApplying Fok1 enzyme to cleavage input:')
            restriction_output = Fok1.cleavage(copy_input)

            print('Results of cleavage:')
            print(restriction_output)
            print(f'\nSticky part: {restriction_output.sticky_part}')

            # Move forward in the original input and set the current state to the previous next state.
            # The new next state should be resolved also from the cleavage output same as done before.
            original_input_index += 1
            current_state = next_state
            next_state = self.__states[restriction_output.sticky_part.__repr__()]
            print(f'\nState Transition: {current_state} '
                  f'--{self.__original_input[original_input_index]}--> {next_state}')

            # Trying to fine compatible transition rule from the current state reached which found at the
            # sticky part of the cleavage output results.
            compatible_transition_rule_molecule = self.__find_compatible_transition_rule(restriction_output)

            # Also, the sticky part of the cleavage results belongs to specific letter representation and therefore
            # connects to compatible transition rule which knows the letter received and applying the transition
            # to the appropriate next state.
            if compatible_transition_rule_molecule is not None:
                print('\nAppropriate transition rule molecule:')
                print(compatible_transition_rule_molecule['molecule'])
                print(f'\nSticky part: {compatible_transition_rule_molecule["molecule"].sticky_part}')

                print('\nCombining cleavage output and transition rule molecules:')
                copy_input = DNAMolecule(
                    Ligase.combine_molecules(
                        compatible_transition_rule_molecule['molecule'].first_strand,
                        restriction_output.first_strand
                    ),
                    Ligase.combine_molecules(
                        compatible_transition_rule_molecule['molecule'].second_strand,
                        restriction_output.second_strand
                    )
                )
                print(copy_input)

            # If no transition rule found, means we reach end of input.
            else:
                print('\nReached end of input.')
                copy_input = restriction_output
                restriction_output = None

        print(f'\nLast state is {next_state} - means the automaton '
              f'{"accept" if next_state in self.__accepting_states else "reject"} this input.')

    @classmethod
    def __build_input_representation(cls, automaton_input):
        """
        Build the input representation for the automaton.

        :param automaton_input: Raw string input the automaton should process.
        :return: DNA Molecule representation of the automaton input.
        """
        first_strand_seq_raw = []

        # Attaching first start letter
        first_strand_seq_raw += cls.__alphabet[automaton_input[0]].bases

        # Adding spacer for the first letter
        first_strand_seq_raw +=\
            DNASequence.create_random_sequence(1).bases \
            if automaton_input[1] == 'a' \
            else DNASequence.create_random_sequence(3).bases

        first_strand_seq_raw += cls.__alphabet[automaton_input[1]].bases

        for index in range(2, len(automaton_input)):
            first_strand_seq_raw += cls.__alphabet[automaton_input[index]].bases

        first_strand_seq = DNASequence(first_strand_seq_raw)
        second_strand_seq = first_strand_seq.get_complement()

        return DNAMolecule(first_strand_seq, second_strand_seq)

    @classmethod
    def __find_compatible_transition_rule(cls, current_sticky_dna):
        """
        Finds compatible transition rule based on the given DNA molecule sticky part.

        :param current_sticky_dna: DNA molecule with sticky part.
        :return: The compatible transition rule if found, Otherwise None.
        """
        for transition_rule in cls.__transition_rules:
            if transition_rule['molecule'].sticky_part.complement_of(current_sticky_dna.sticky_part):
                return transition_rule

        return None


if __name__ == '__main__':
    automaton = RestrictionEnzymeAutomata('abbabab')
    automaton.process_automaton()
