from random import random


def chance(percentage_chance):
    """
    Determine percentage chance to occur

    :param percentage_chance: Percentage of the chance to occur
    :return: True if chance occur, Otherwise False
    """
    return random() < percentage_chance
