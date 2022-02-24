from molecules.dna_sequence import DNASequence
from molecules.dna_molecule import DNAMolecule
from enzymes.polymerase import Polymerase
from enzymes.primer import Primer
from utils.probabilites import chance


class PCR:
    """
    Represents PCR machine and performs the equivalent procedure the PCR machine does with DNA molecules.
    """

    __chance_primer_attach = 0.95
    __chance_polymerase_complement = 0.95

    def __init__(self, dna_molecules, first_strand_primer, second_strand_primer):
        """
        Initialize PCR machine with DNA molecules and primers.

        :param dna_molecules: DNA molecules to operate on.
        :param first_strand_primer: First strands primer which should attach in the 5' direction.
        :param second_strand_primer: Second strands primer which should attach in the 3' direction.
        """
        self.__current_dna_molecules = dna_molecules
        self.__first_strand_primer = first_strand_primer
        self.__second_strand_primer = second_strand_primer

    @property
    def dna_molecules(self):
        return self.__current_dna_molecules

    def perform_operation(self):
        """
        Performing the PCR operation.
        Storing the results DNA molecules in the DNA molecules property.
        """
        print('\nCurrent DNA molecules:')
        self.__print_dna_parts(self.__current_dna_molecules)

        print('Separating the DNA molecules strands into DNA sequences')
        first_strands, second_strands = self.__separate_dna_molecules()

        print('\nSeparated Strands:')
        self.__print_dna_parts(first_strands, second_strands)

        print('Attaching the primers')
        self.__attach_primers(first_strands, second_strands)

        print('\nAfter attaching the primers, DNA strands:')
        self.__print_dna_parts(first_strands, second_strands)

        print('Applying polymerases to duplicate selection DNA portions')
        self.__current_dna_molecules = self.__apply_polymerases(first_strands, second_strands)

        print('\nResults DNA molecules after polymerases:')
        self.__print_dna_parts(self.__current_dna_molecules)

    def __separate_dna_molecules(self):
        """
        Separates the DNA molecules strands by heat as the PCR machine should do.

        :return: Separated strands divided into 2 groups:
        List of first strands and list of second strands.
        """
        first_strands = []
        second_strands = []

        for dna_molecule in self.__current_dna_molecules:
            first_strands.append(dna_molecule.first_strand)
            second_strands.append(dna_molecule.second_strand)

        return first_strands, second_strands

    def __attach_primers(self, first_strands, second_strands):
        """
        Attaching primers to mark specific gene in the given DNA strands.

        :param first_strands: First strands of DNA sequences.
        :param second_strands: Second strands of DNA sequences.
        """
        # Trying to attach primer to each first strands
        for strand in first_strands:

            # Chance of primer attached without errors
            if chance(self.__chance_primer_attach):
                strand.attach_primer(self.__first_strand_primer)

        # Trying to attach primer to each second strands
        for strand in second_strands:

            # Chance of primer attached without errors (by the chance of success)
            if chance(self.__chance_primer_attach):
                strand.attach_primer(self.__second_strand_primer)

    def __apply_polymerases(self, first_strands, second_strands):
        """
        Applying polymerases on the first and second strands by the point of the primers attachments.

        :param first_strands: First strands of DNA sequences.
        :param second_strands: Second strands of DNA sequences.
        :return: New DNA molecules created by the polymerase process.
        """
        dna_molecules = []

        # For each of the first strands, which are in the 3' direction of the DNA molecule
        # Apply the polymerase enzyme
        for strand in first_strands:

            # Checking if primer is attached and there's no error (by the chance of success)
            if strand.attached_primer and chance(self.__chance_polymerase_complement):

                # Because the 3' direction of the strand, the polymerase will copy from the primer
                # attachment point, but in the reverse direction of array indices.
                # Here, we save time for that and start from the first index until we reach the primer
                # attachment point, leading to the same results.
                dna_molecules.append(
                    DNAMolecule(
                        DNASequence(
                            strand.bases[
                                0:strand.attached_primer_point + len(strand.attached_primer.searched_sequence)
                            ]
                        ),
                        DNASequence(
                            Polymerase.create_complement_sequence(
                                DNASequence(
                                    strand.bases[
                                        0:strand.attached_primer_point + len(strand.attached_primer.searched_sequence)
                                    ]
                                )
                            )
                        )
                    )
                )

        # For each of the second strands, which are in the 5' direction of the DNA molecule
        # Apply the polymerase enzyme
        for strand in second_strands:

            # Checking if primer is attached and there's no error (by the chance of success)
            if strand.attached_primer and chance(self.__chance_polymerase_complement):

                # Because the 5' direction of the strand, it is easy to complete the strand as
                # the default direction of the list indices from the attached primer point.
                dna_molecules.append(
                    DNAMolecule(
                        DNASequence(
                            Polymerase.create_complement_sequence(
                                DNASequence(strand.bases[strand.attached_primer_point:])
                            )
                        ),
                        DNASequence(strand.bases[strand.attached_primer_point:])
                    )
                )

        return dna_molecules

    @staticmethod
    def __print_dna_parts(dna_parts, dna_parts_other=None):
        """
        Pretty print DNA parts.

        :param dna_parts: DNA molecules or sequences to print.
        :param dna_parts_other: DNA molecules or sequences to print next to the other DNA parts.
        """
        if dna_parts_other is None:
            for dna_part in dna_parts:
                print(f'{dna_part}\n')
        else:
            for index in range(len(dna_parts)):
                print(f'{dna_parts[index]}\n<--\n{dna_parts_other[index]}\n-->\n')


if __name__ == '__main__':
    sequence1 = DNASequence('AAAATTTTGGGGCCCC')
    comp_sequence1 = sequence1.get_complement()

    sequence2 = DNASequence('AAAATGGGATATATATCCCC')
    comp_sequence2 = sequence2.get_complement()

    first_strand_primer_selection = Primer('GGGG')
    second_strand_primer_selection = Primer('TGGG')

    molecules = [
        DNAMolecule(sequence1, comp_sequence1),
        DNAMolecule(sequence2, comp_sequence2)
    ]

    pcr = PCR(
        molecules,
        first_strand_primer_selection,
        second_strand_primer_selection
    )

    pcr.perform_operation()
