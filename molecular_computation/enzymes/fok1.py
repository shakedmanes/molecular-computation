from molecules.dna_sequence import DNASequence
from molecules.dna_molecule import DNAMolecule


class Fok1:

    __interaction_sequence = DNASequence('GGATG')
    __first_strand_cut = 9
    __second_strand_cut = 13
    __sticky_part_size = 4

    @classmethod
    def cleavage(cls, dna_molecule):
        """
        Cleavage part of DNA molecule due to interaction sequence of Fok1 enzyme.

        :param dna_molecule: DNA molecule to operate on.
        :return: New DNA molecule after cleavage by the Fok1 enzyme if the interaction point found, Otherwise None.
        """
        if not isinstance(dna_molecule, DNAMolecule):
            raise ValueError(f'Fok1 only operates on DNA molecules. Given molecule: {DNAMolecule}')

        index_at_interaction = 0
        index_at_molecule = 0

        while index_at_molecule < dna_molecule.first_strand.length and \
              index_at_interaction < cls.__interaction_sequence.length:

            # If part of interaction point found, continue this path
            if dna_molecule.first_strand.bases[index_at_molecule] == cls.__interaction_sequence.bases[index_at_interaction]:
                index_at_interaction += 1

            # Otherwise, back to the starting point
            else:
                index_at_interaction = 0

            # Continue the path in the DNA molecule
            index_at_molecule += 1

        try:
            # Found the interaction sequence in the molecule
            if index_at_interaction == cls.__interaction_sequence.length:
                first_strand_cut_point = index_at_molecule + cls.__first_strand_cut
                second_strand_cut_point = index_at_molecule + cls.__second_strand_cut

                # Returning new DNa molecule after removing the necessary parts by the interaction
                return DNAMolecule(
                    DNASequence(dna_molecule.first_strand.bases[first_strand_cut_point:]),
                    DNASequence(
                        DNASequence.create_sticky_part(cls.__sticky_part_size).bases +
                        dna_molecule.second_strand.bases[second_strand_cut_point:]
                    )
                )
        except Exception as exp:
            # Cleavage failed by lower DNA sequence
            return None

        # Not found any interaction sequence, returning none cleavage
        return None
