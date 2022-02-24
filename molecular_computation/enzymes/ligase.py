from molecules.dna_sequence import DNASequence


class Ligase:
    """
    Represent the Ligase enzyme which attach 2 DNA sequences to each other.
    """

    @staticmethod
    def combine_molecules(dna_sequence_one, dna_sequence_two):
        """
        Combining two dna sequence into one sequence.
        (The first dna sequence can be dna molecule, so the other dna sequence can stick to its sticky part)

        :param dna_sequence_one:
        :param dna_sequence_two:
        :return: DNA molecule which connected by the two dna sequences given.
        """
        if isinstance(dna_sequence_one, DNASequence) and\
           isinstance(dna_sequence_two, DNASequence):

            # Find the sticky part to remove from the first sequence at end
            remove_from_end_index = dna_sequence_one.length
            for index in range(dna_sequence_one.length - 1, 0, -1):
                if dna_sequence_one.bases[index] in DNASequence.sticky_bases:
                    remove_from_end_index -= 1
                else:
                    break

            # Find the sticky part to remove from the second sequence at start
            remove_from_start_index = 0
            for index in range(dna_sequence_two.length):
                if dna_sequence_two.bases[index] in DNASequence.sticky_bases:
                    remove_from_start_index += 1
                else:
                    break

            return DNASequence(
                dna_sequence_one.bases[:remove_from_end_index] +
                dna_sequence_two.bases[remove_from_start_index:]
            )
