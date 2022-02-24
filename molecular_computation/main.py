from os import system

options = {
    1: {
        'display': '3-SAT Algorithm',
        'path': './algorithms/three_sat_algorithm.py',
    },
    2: {
        'display': 'Restriction Enzyme Automata',
        'path': './algorithms/restriction_enzyme_automata.py',
    },
    3: {
        'display': 'Gel Electrophoresis',
        'path': './procedures/gel_electrophoresis.py',
    },
    4: {
        'display': 'PCR',
        'path': './procedures/pcr.py',
    }
}


def print_selection_options():
    print('Choose which virtual computation would you want to run:')
    for option in options.keys():
        print(f'({option}) {options[option]["display"]}')


if __name__ == '__main__':
    option_range = max(options.keys())
    selection = 0

    while selection < 1 or selection > 4:
        print_selection_options()
        try:
            selection = int(input('Your choice: '))
        except Exception as exp:
            selection = 0

    print(f'Running {options[selection]["display"]}...\n\n')

    system(f'python {options[selection]["path"]}')
