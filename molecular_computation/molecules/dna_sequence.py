from re import compile
from random import choice


class DNASequence:
    """
    Represent a DNA Sequence strand which constructed from ACGT bases.
    """

    __dna_bases_regex = compile(r'[AaTtCcGgSs]+')

    __dna_bases = 'ACTG'

    sticky_bases = ['S', 's']

    complement_table = {
        'a': 't',
        'A': 'T',
        't': 'a',
        'T': 'A',
        'c': 'g',
        'C': 'G',
        'g': 'c',
        'G': 'C'
    }

    def __init__(self, sequence):
        """
        Creates a DNA sequence by given sequence of DNA bases.

        :param sequence: String/List/DNA Sequence object represents DNA sequence.
        """
        if self.__validate_sequence(sequence):
            self.__bases = self.__parse_sequence(sequence)
            self.__attached_primer = None
            self.__attached_primer_point = None
        else:
            raise ValueError('DNA sequence must contains ACGT bases and be in format of String / List / DNASeqeunce.')

    @property
    def bases(self):
        return self.__bases

    @property
    def length(self):
        return len(self.__bases)

    @property
    def attached_primer(self):
        return self.__attached_primer

    @property
    def attached_primer_point(self):
        return self.__attached_primer_point

    @property
    def sticky_part_index_and_size(self):
        size = 0
        first_index = None

        # Searching for sticky part in the DNA sequence
        for index in range(len(self.bases)):

            # If reach sticky part, need to increase the size of it
            if self.bases[index] in self.sticky_bases:
                size += 1

                # If it is the first time we reach sticky part, mark it start point
                if first_index is None:
                    first_index = index

            # If the sequence of sticky parts is over, stop searching for another sticky parts
            if first_index is not None and self.bases[index] not in self.sticky_bases:
                break

        # If found sticky part in the DNA sequence, return it starting point and size
        if first_index is not None:
            return first_index, size

        # Otherwise return None for both
        return None, None

    def get_complement(self):
        """
        Returns the complement DNA Sequence.

        :return: Complement DNA Sequence
        """
        return DNASequence([self.complement_table.get(base) for base in self.bases])

    def complement_of(self, other_dna_seq):
        """
        Determine if other DNA Sequence is the complement DNA Sequence.

        :param other_dna_seq: DNA Sequence object.
        :return: True if the DNA Sequence given is complement, Otherwise False.
        """
        if isinstance(other_dna_seq, self.__class__):
            reverse_comp = [
                self.complement_table[base]
                if base not in self.sticky_bases
                else base
                for base in other_dna_seq.bases
            ]

            return reverse_comp == self.bases

        return False

    def complement_from(self, other_dna_seq):
        """
        Find the complement point on the current sequence by other sequence.

        :param other_dna_seq: Other DNA Sequence.
        :return: Index of the complement point or None if does not have complement point.
        """
        if not isinstance(other_dna_seq, self.__class__):
            return None

        curr_pos_seq = 0
        curr_pos_other = 0

        while curr_pos_seq < len(self.bases) and curr_pos_other < len(other_dna_seq.bases):

            # Matching point, the DNA sequences complement each other
            if self.bases[curr_pos_seq] in self.sticky_bases or \
               other_dna_seq.bases[curr_pos_other] in self.sticky_bases or \
               self.complement_table.get(self.bases[curr_pos_seq]) == other_dna_seq.bases[curr_pos_other]:
                curr_pos_seq += 1
                curr_pos_other += 1

            # No matching point, but continue with the DNA sequence and back to first base in other sequence.
            else:
                curr_pos_seq += 1
                curr_pos_other = 0

        # If the whole other sequence matched, return the difference between the last location
        # visited with the size of the other sequence to reach first complement point.
        if curr_pos_other == len(other_dna_seq.bases):
            return curr_pos_seq - len(other_dna_seq.bases)

        # If we reach end of current sequence but the other sequence was match in part of it's bases,
        # return the difference between the current reach point at the other sequence between the length of
        # self bases.
        if curr_pos_seq == len(self.bases) and curr_pos_other > 0:
            return len(self.bases) - curr_pos_other

        return None

    def attach_primer(self, primer):
        """
        Trying to attach primer to the DNA sequence.

        :param primer: Primer to try to attach the DNA sequence.
        """
        if self.__attached_primer is None:
            attachment_point = primer.try_attach_to_dna_sequence(self)
            if attachment_point is not None:
                self.__attached_primer = primer
                self.__attached_primer_point = attachment_point

    def __validate_sequence(self, sequence):
        """
        Validates if a given sequence is valid as DNA sequence.<br>
        It must be string or list or actual DNA Sequence object which contains the DNA bases.

        :param sequence: String/List/DNA Sequence object represents DNA sequence.
        :return: True if the sequence passed is DNA sequence, Otherwise False.
        """
        return (
            (type(sequence) == str and self.__dna_bases_regex.fullmatch(sequence)) or
            (type(sequence) == list and self.__dna_bases_regex.fullmatch(''.join(sequence))) or
            isinstance(sequence, self.__class__)
        )

    def __parse_sequence(self, sequence):
        """
        Parses a valid sequence into DNA sequence of bases (list of bases).

        :param sequence: String/List/DNA Sequence object represents DNA sequence.
        :return: List containing the DNA sequence given.
        """
        if type(sequence) == str or type(sequence) == list:
            return list(sequence)
        elif isinstance(sequence, self.__class__):
            return list(sequence.bases)
        else:
            raise ValueError('Cannot parse sequence - must be type of String / List / DNA sequence.')

    @classmethod
    def create_random_sequence(cls, size=20):
        """
        Create random DNA sequence by given size.

        :param size: Size of the DNA sequence to create.
        :return: Random DNA sequence of the given size.
        """
        return cls([choice(cls.__dna_bases) for _ in range(size)])

    @classmethod
    def create_sticky_part(cls, size=4):
        """
        Creates DNA sequence of sticky part by given size.

        :param size: Size of the DNA sequence of sticky part to create.
        :return: DNA sequence contains only sticky parts of the given size.
        """
        return cls(cls.sticky_bases[0] * size)

    def __eq__(self, other):
        """
        Equality checking of 2 DNA sequences.

        :param other: DNA sequence object.
        :return: True if both DNA sequences represents the same sequence, Otherwise False.
        """
        return (
            isinstance(other, self.__class__) and
            self.bases == other.bases
        )

    def __repr__(self):
        """
        Return string representation of DNA sequence

        :return: String representation of DNA sequence
        """
        if self.attached_primer:
            return ' ' * self.attached_primer_point + \
                   self.attached_primer.searched_sequence + \
                   '' * (self.length - self.attached_primer_point) + \
                   '\n' +\
                   ''.join(self.bases)
        return ''.join(self.bases)
