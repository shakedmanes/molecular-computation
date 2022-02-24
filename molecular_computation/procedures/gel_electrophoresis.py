from molecules.dna_molecule import DNAMolecule
from molecules.dna_sequence import DNASequence


class GelElectrophoresis:
    """
    Produce the Gel Electrophoresis procedure to sort DNA molecules by their size.
    """

    @staticmethod
    def run_gel(dna_molecules):
        """
        Runs the Gel Electrophoresis procedure to sort DNA molecules by their size.

        :param dna_molecules: DNA molecules.
        :return: Sorted list of the DNA molecules given.
        """
        molecules = list(dna_molecules)
        molecules.sort(key=lambda mol: mol.length)
        return molecules


if __name__ == '__main__':
    dna_sequences = [
        DNASequence.create_random_sequence(size=20),
        DNASequence.create_random_sequence(size=10),
        DNASequence.create_random_sequence(size=30),
        DNASequence.create_random_sequence(size=5),
        DNASequence.create_random_sequence(size=15)
    ]

    ex_dna_molecules = [
        DNAMolecule(dna_sequences[index], dna_sequences[index].get_complement())
        for index in range(len(dna_sequences))
    ]

    print('DNA molecules:')
    for molecule in ex_dna_molecules:
        print(f'{molecule}\n')

    print('\nRun Gel Electrophoresis on DNA molecules:')
    gel_dna_molecules = GelElectrophoresis.run_gel(ex_dna_molecules)

    print('Results DNA molecules:')
    for molecule in gel_dna_molecules:
        print(f'{molecule}\n')


