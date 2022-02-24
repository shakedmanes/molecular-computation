

class Polymerase:
    """
    Represent the polymerase enzyme which completes a single strand of DNA,
    which is used as a template for a double-stranded DNA molecule.
    """

    @staticmethod
    def create_complement_sequence(dna_sequence):
        """
        Creates complement DNA sequence for a given DNA sequence.

        :param dna_sequence: DNA sequence to get complement DNA sequence for.
        :return: Complement DNA sequence for given DNA sequence.
        """
        return dna_sequence.get_complement()
