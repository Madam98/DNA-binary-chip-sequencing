#!python
#cython: language_level=3
from abc import abstractmethod, ABC


class Element(ABC):
    @abstractmethod
    def __init__(self):
        self.found_max_cover = []
        self.fitness = self.evaluate_function()  # <-----------zwracamy wynik dopasowanych obiektów

    def mutation(self):
        self._perform_mutation()

    @abstractmethod
    def _perform_mutation(self):  # <------------ sposob mutacji
        pass

    @abstractmethod
    def crossover(self, element2: 'Element', number_of_crossover) -> 'Element':  # <------------ sposob krzyzowania
        pass

    @abstractmethod
    def evaluate_function(self):  # <------------ sposob oceny uzyskanego wyniku
        pass

    '''
    @abstractmethod
    def global_evaluate_function(self, element2: 'Element') -> 'Element':
        pass
    '''
