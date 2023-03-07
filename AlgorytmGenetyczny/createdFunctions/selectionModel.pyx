#!python
#cython: language_level=3

def elite_selection_model(generation, selection_percent):
    max_selected = int(len(generation) * selection_percent)
    sorted_by_assess = sorted(generation, key=lambda x: x.fitness, reverse=True)
    return sorted_by_assess[:max_selected]