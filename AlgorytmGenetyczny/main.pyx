#!python
#cython: language_level=3

#**********************************************************************************
import pyximport
pyximport.install()

from createdClasses.colors        import Bcolors
from createdClasses.element       import Element
from createdClasses.genethicClass import GeneticAlgorithm
from createdClasses.text          import Text
from createdClasses.queue         import Queue

#**********************************************************************************

from createdFunctions.read             import parse_file
from createdFunctions.firstPopulation  import first_sequence_k_minus_1_for_both_poll, \
                                              look_for_indexes_max, \
                                              combine_string
from createdFunctions.log              import log, DEBUG, end_time
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

from createdThread.threadPopulation import run_thread_population

#**********************************************************************************
import itertools
import sys
import collections
import time
import threading
from random import random, randint, choice, sample, randrange, shuffle
import pandas as pd
import time

from main_iter import fun_iter

# --------------------------------------------------------------
# ----- MAIN PROGRAM -------------------------------------------
# --------------------------------------------------------------

#***************************************************************
population_number    = 10
mutation_probability = 0.01
size_sequence        = 500
cross_number         = 3                    #<---- nie zmieniac!
selection_percent    = 0.5
TIMEOUT              = 1000
number_of_average    = 2
run_tests            = False
#***************************************************************

pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_rows', None)

if not(run_tests):
    start, seq_len, probe1_pattern, probe2_pattern, ols1, ols2, s1df, s2df, s1, s2 = parse_file(
        "sample_cases/positive_numbers_max_k_max_sqt_max_pos/" + str(size_sequence) + ".xml")
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
                          size_sequence,
                          cross_number,
                          selection_percent,
                          mutation_probability,
                          TIMEOUT)
    ga.run()
    final_string, final_time = ga.return_class_string_time()
    print(final_string)
    print(final_time)

######################################################
#
# T E S T
#
######################################################

else:
    item_size_sequence   = [30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]
    item_mutation_prob   = [0.00, 0.05, 0.1, 0.15]
    item_population_size = [40]

    #item_size_sequence = [30, 40]
    #item_mutation_prob = [0.00]
    #item_population_size = [10, 20]

    time_result_gene = []

    #-----------------------------------------------------------------

    print("Testy algorytm genetyczny")
    for i_1 in range(len(item_population_size)):
        population_number = item_population_size[i_1]
        for i_2 in range(len(item_mutation_prob)):
            mutation_probability = item_mutation_prob[i_2]
            end_list_times = []
            for i_3 in range(len(item_size_sequence)):
                size_sequence = item_size_sequence[i_3]
                start, seq_len, probe1_pattern, probe2_pattern, ols1, ols2, s1df, s2df, s1, s2 = parse_file("sample_cases/positive_numbers_max_k_max_sqt_max_pos/" + str(size_sequence) + ".xml")
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
                log(''.join((Bcolors.OKCYAN, 'Uruchamiam klasę genethicAlghorithm...\n', Bcolors.ENDC)), DEBUG)
                sum_of_time = 0
                for j in range(number_of_average):
                    string_ga = "ga"
                    string_temp = string_ga + str(j)
                    ols1 = {x: False for x in s1}
                    ols2 = {x: False for x in s2}
                    string_temp = GeneticAlgorithm(first_sequence_k_minus_1_for_both_poll,
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
                                          size_sequence,
                                          cross_number,
                                          selection_percent,
                                          mutation_probability,
                                          TIMEOUT)
                    string_temp.run()
                    final_string_genethic, time_gen = string_temp.return_class_string_time()
                    sum_of_time += time_gen
                end_time_gen = sum_of_time / number_of_average
                end_list_times.append(end_time_gen)
            print()
            print("**************************************")
            print(Bcolors.LIGHT_WHITE + "Rozmiar populacji:          " + Bcolors.LIGHT_PURPLE + str(item_population_size[i_1]) + Bcolors.ENDC)
            print(Bcolors.LIGHT_WHITE + "Prawdopodobieństwo mutacji: " + Bcolors.LIGHT_PURPLE + str(item_mutation_prob[i_2]) + Bcolors.ENDC)
            for i_4 in range(len(item_size_sequence)):
                print(Bcolors.LIGHT_WHITE + "* Sekwencja n: " + Bcolors.LIGHT_PURPLE + str(
                    item_size_sequence[i_4]) + Bcolors.LIGHT_WHITE + "  Czas: " + Bcolors.RED + str(
                    round(end_list_times[i_4], 3)) + Bcolors.ENDC)
            print("**************************************\n")
    
    #-----------------------------------------------------------------



    #-----------------------------------------------------------------
    print("\n\n")
    print("ILOŚĆ TESTÓW: " + str(number_of_average))
    print("\n\nAlgorytm dokładny")
    print("*************************************")
    for grow_sequence in range(len(item_size_sequence)):
        sum_time = 0
        for number_of_test in range(number_of_average):
            result, result_time = fun_iter(item_size_sequence[grow_sequence])
            sum_time += result_time
        sum_time = sum_time / number_of_average
        print(Bcolors.LIGHT_WHITE + "* Sekwencja n: " + Bcolors.LIGHT_PURPLE + str(item_size_sequence[grow_sequence]) + Bcolors.LIGHT_WHITE + "  Czas: " + Bcolors.RED + str(round(sum_time, 3)) + Bcolors.ENDC)
    print("*************************************\n\n")
    #-----------------------------------------------------------------