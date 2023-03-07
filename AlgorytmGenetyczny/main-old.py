#**********************************************************************************

from createdClasses.colors        import Bcolors
from createdClasses.element       import Element
from createdClasses.genethicClass import GeneticAlgorithm
from createdClasses.text          import Text
from createdClasses.queue         import Queue

#**********************************************************************************

from createdFunctions.read             import parse_file
from createdFunctions.firstPopulation  import first_sequence_k_minus_1_for_both_poll, \
                                              look_for_indexes_max
from createdFunctions.log              import log, DEBUG
from createdFunctions.stringOperations import find_same_last_letter, \
                                              duplicate, \
                                              endOverlap, \
                                              last_element_decode, \
                                              decode, \
                                              scale, \
                                              merge
from createdFunctions.selectionModel   import elite_selection_model
from createdFunctions.stopCondition    import stop_condition

#**********************************************************************************

from createdThread.threadPopulation import

#**********************************************************************************
import itertools
import sys
import collections
import time
import threading
from random import random, randint, choice, sample, randrange, shuffle
import pandas as pd
import time



# --------------------------------------------------------------
# ----- MAIN PROGRAM -------------------------------------------
# --------------------------------------------------------------
population_number = 30
mutation_probability = 0.1

pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_rows', None)
start, seq_len, probe1_pattern, probe2_pattern, ols1, ols2, s1df, s2df = parse_file("sample_cases/bio500.xml")
Z_set_first, P_set_first = decode(start, True)
k_minus_1 = len(start) - 1

log(''.join((Bcolors.OKGREEN, 'Informacje początkowe:', Bcolors.ENDC)), DEBUG)
log(''.join(('Start:\t\t\t\t\t\t|\t', start)), DEBUG)
log(''.join(('Długość szukanej sekwencji:\t|\t', str(seq_len))), DEBUG)
log(''.join(('Proba1:\t\t\t\t\t\t|\t', str(probe1_pattern))), DEBUG)
log(''.join(('Proba2:\t\t\t\t\t\t|\t', str(probe2_pattern))), DEBUG)
log(''.join(('Pierwsza sekwencja Z:\t\t|\t', Z_set_first)), DEBUG)
log(''.join(('Pierwsza sekwencja P:\t\t|\t', P_set_first)), DEBUG)
log(''.join(('Pokrycie k minus jeden:\t\t|\t', str(k_minus_1))), DEBUG)
log(''.join(('S1 lista:\n', str(s1df.T))), DEBUG)
log(''.join(('S2 lista:\n', str(s2df.T))), DEBUG)
log(''.join(('Powstały słownik wykorzystanych sekwencji dla Z:\t', str(ols1))), DEBUG)
log(''.join(('Powstały słownik wykorzystanych sekwencji dla P:\t', str(ols2))), DEBUG)
log('\n\n\n', DEBUG)
log(''.join((Bcolors.OKCYAN, 'Uruchamiam klasę genethicAlghorithm...\n', Bcolors.ENDC)), DEBUG)

ga = GeneticAlgorithm(first_sequence_k_minus_1_for_both_poll,
                      elite_selection_model,
                      stop_condition,
                      start,
                      seq_len,
                      k_minus_1,
                      s1df,
                      s2df,
                      Z_set_first,
                      P_set_first,
                      ols1,
                      ols2,
                      population_number,
                      mutation_probability)
ga.run()


#Queue.size()

#queue = Queue()
#queue2 = Queue()



