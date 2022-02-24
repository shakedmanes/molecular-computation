from molecules.dna_sequence import DNASequence


class Primer:
    """
    Represent a Primer DNA Sequence which can be attached to DNA sequence in order to mark specific
    segment in the DNA sequence.
    """

    def __init__(self, searched_sequence):
        """
        Initialize Primer sequence.

        :param searched_sequence: DNA Sequence to search.
        """
        self.__searched_sequence = searched_sequence

    @property
    def searched_sequence(self):
        return self.__searched_sequence

    def try_attach_to_dna_sequence(self, dna_sequence):
        """
        Try attach to a given DNA Sequence.

        :param dna_sequence: DNA Sequence to try attach to.
        :return: The index in the DNA Sequence where the attachment occurred, Otherwise None.
        """
        if not isinstance(dna_sequence, DNASequence):
            raise ValueError('DNA Sequence given is not instance of DNASequence.')

        # Mark if attachment point found
        found_attachment_point = False

        # Search the complement sequence of the search sequence in the DNA sequence
        complement_search_sequence = [DNASequence.complement_table.get(base) for base in self.__searched_sequence]

        # Current index navigated
        current_index = 0

        # Iterate all sub sequences in the DNA sequence to find the complement sequence of the searched sequence
        while current_index <= dna_sequence.length - len(self.__searched_sequence):

            # If the complement sequence of the search sequence is found, break and return the index
            if dna_sequence.bases[current_index: current_index + len(self.__searched_sequence)] ==\
               complement_search_sequence:
                found_attachment_point = True
                break

            current_index += 1

        # If found attachment point, return it
        if found_attachment_point:
            return current_index

        return None

    def __eq__(self, other):
        """
        Equality checking of 2 Primers.

        :param other: Primer object.
        :return: True if both Primers represents the same search sequence, Otherwise False.
        """
        return isinstance(other, self.__class__) and self.searched_sequence == other.searched_sequence
