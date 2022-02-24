from molecules.dna_sequence import DNASequence


class DNAMolecule:
    """
    Represent DNA Molecule containing two DNA Sequence strands.
    """

    def __init__(self, dna_seq_one, dna_seq_two):
        if isinstance(dna_seq_one, DNASequence) and \
           isinstance(dna_seq_two, DNASequence) and \
           (complement_pos := dna_seq_one.complement_from(dna_seq_two)) is not None:
            self.__first_strand = dna_seq_one
            self.__second_strand = dna_seq_two
            self.__complement_pos = complement_pos
            self.__sticky_part_strand = None
            self.__sticky_part_pos = None
            self.__sticky_part_size = None

            first_strand_sticky_index, first_strand_sticky_size = self.first_strand.sticky_part_index_and_size
            second_strand_sticky_index, second_strand_sticky_size = self.second_strand.sticky_part_index_and_size

            if first_strand_sticky_index is not None and first_strand_sticky_size is not None:
                self.__sticky_part_strand = 'second_strand'
                self.__sticky_part_pos = first_strand_sticky_index
                self.__sticky_part_size = first_strand_sticky_size
            elif second_strand_sticky_index is not None and second_strand_sticky_size is not None:
                self.__sticky_part_strand = 'first_strand'
                self.__sticky_part_pos = second_strand_sticky_index
                self.__sticky_part_size = second_strand_sticky_size

        else:
            raise ValueError(
                f'The two given DNA sequences are not complement to each other or not even a DNA sequences.\n'
                f'DNA sequence one: {dna_seq_one}\n'
                f'DNA sequence two: {dna_seq_two}\n'
            )

    @property
    def first_strand(self):
        return self.__first_strand

    @property
    def second_strand(self):
        return self.__second_strand

    @property
    def sticky_part_strand(self):
        return self.__sticky_part_strand

    @property
    def sticky_part_index(self):
        return self.__sticky_part_pos

    @property
    def sticky_part_size(self):
        return self.__sticky_part_size

    @property
    def sticky_part(self):
        return DNASequence(
            self.__getattribute__(
                self.sticky_part_strand
            ).bases[self.sticky_part_index: self.sticky_part_index + self.sticky_part_size]
        )

    @property
    def length(self):
        return max(self.first_strand.length, self.second_strand.length)

    def __eq__(self, other):
        """
        Equality checking of 2 DNA molecules.

        :param other: DNA molecule object.
        :return: True if both DNA molecules represents the same sequence, Otherwise False.
        """
        return (
                isinstance(other, self.__class__) and
                (
                        (self.first_strand == other.first_strand and self.second_strand == other.second_strand) or
                        (self.first_strand == other.second_strand and self.second_strand == other.first_strand)
                )
        )

    def __repr__(self):
        """
        Return string representation of DNA molecule

        :return: String representation of DNA molecule
        """
        return f'{self.first_strand.__repr__()}\n{self.second_strand.__repr__()}'
