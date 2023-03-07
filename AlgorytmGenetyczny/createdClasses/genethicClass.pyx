#!python
# cython: language_level=3
from createdClasses.colors              import Bcolors
from createdClasses.text                import Text
from createdClasses.element             import Element
from createdFunctions.log               import log, \
                                               DEBUG, \
                                               st, \
                                               end_time
from createdThread.threadPopulation     import run_thread_population, \
                                               threading
from createdFunctions.firstPopulation   import create_permutation, \
                                               combine_string
from createdFunctions.stringOperations  import last_element_decode, \
                                               endOverlap, \
                                               decode
from random                             import random
'''
from multiprocessing                    import Pool, \
                                               Process, \
                                                Value, Manager
'''
#from multiprocessing import *
import multiprocessing as mp
import pandas as pd
import time
import copy
import numpy
from functools import partial
from itertools import repeat

pd.set_option("display.max_colwidth", None)

# --------------------------------------------------------------
# ----- GENETIC CLASS ------------------------------------------
# --------------------------------------------------------------
class GeneticAlgorithm:
    def __init__(self, first_sequence: callable,
                 selection_model: callable,
                 stop_condition: callable,
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
                 selection_percent: float = 0.2,
                 mutation_probability: float = 0.05,
                 TIMEOUT = 60):
        self.first_generation_func = first_sequence
        self.selection_model = selection_model
        self.stop_condition = stop_condition
        self.start = start
        self.seq_len = seq_len
        self.k_minus_1 = k_minus_1
        self.s1df = s1df
        self.s2df = s2df
        self.Z_set_first = Z_set_first
        self.P_set_first = P_set_first
        self.ols1 = ols1
        self.ols2 = ols2
        self.population_number = population_number
        self.size_sequence = size_sequence
        self.cross_number = cross_number
        self.selection_percent = selection_percent
        self.mutation_probability = mutation_probability
        self.TIMEOUT = TIMEOUT

    def run(self):
        FIRST_POP_DF, SECOND_POP_DF, self.first_can, self.second_can, self.cross_number = self.first_generation_func(self.s1df, self.s2df,
                                                                                                  self.ols1, self.ols2,
                                                                                                  self.k_minus_1,
                                                                                                  self.population_number,
                                                                                                  self.cross_number,
                                                                                                  self.size_sequence)

        dont_start_first, dont_start_second = False, False
        if len(self.first_can) == 1:
            print("Znaleziona sekwencja przed uruchomieniem algorytmu genetycznego!")
            print(str(self.first_can['SEQUENCE']))
            dont_start_first = True
        if len(self.second_can) == 1:
            print("Znaleziona sekwencja przed uruchomieniem algorytmu genetycznego!")
            print(str(self.second_can['SEQUENCE']))
            dont_start_second = True

        if True:

            final_string, time = population_threads(FIRST_POP_DF, SECOND_POP_DF, self)
            #print("Wynik zakodowany:     " + str(final_string))
            self.final_string = final_string
            self.time = time

        else:
            populationSW, len_populationSW, best_SW, length_best_SW = init_population(FIRST_POP_DF, self.first_can,
                                                                                      self)
            populationRY, len_populationRY, best_RY, length_best_RY = init_population(SECOND_POP_DF, self.second_can,
                                                                                      self)
            i_generation, counter_SW, counter_RY = 0, 0, 0
            len_seq_str_SW = 1000000
            len_seq_str_RY = 1000000

            while True:
                new_populationSW, new_populationRY = select_model_SW_RY(populationSW, populationRY, self.selection_model)
                while len(new_populationSW) != len_populationSW:
                    while_mutation(populationSW, populationRY, new_populationSW, new_populationRY, self.cross_number,
                                   self.mutation_probability)
                best_SW, best_RY = max(new_populationSW, key=lambda x: x.fitness), max(new_populationRY,
                                                                                       key=lambda x: x.fitness)
                if best_SW.found_max_cover:
                    populationSW, self.first_can, best_SW, counter_SW = if_found_max_cover(best_SW.found_max_cover,
                                                                                           self.first_can, self.k_minus_1,
                                                                                           self.population_number,
                                                                                           self.Z_set_first,
                                                                                           self.P_set_first,
                                                                                           self.size_sequence)
                elif len_seq_str_SW > len(self.first_can):
                    show_generation_info(i_generation, best_SW, best_SW.fitness, best_SW.seq_parts_str, Bcolors.FAIL)
                    len_seq_str_SW = len(self.first_can)
                    if best_SW.fitness > length_best_SW:
                        populationSW = new_populationSW
                if best_RY.found_max_cover:
                    populationRY, self.second_can, best_RY, counter_RY = if_found_max_cover(best_RY.found_max_cover,
                                                                                            self.second_can, self.k_minus_1,
                                                                                            self.population_number,
                                                                                            self.Z_set_first,
                                                                                            self.P_set_first,
                                                                                            self.size_sequence)
                elif len_seq_str_RY > len(self.second_can):
                    show_generation_info(i_generation, best_RY, best_RY.fitness, best_RY.seq_parts_str, Bcolors.OKGREEN)
                    len_seq_str_RY = len(self.second_can)
                    if best_RY.fitness > length_best_RY:
                        populationRY = new_populationRY

                if counter_SW > 10:
                    first_population_df = pd.DataFrame(
                        list(create_permutation(self.first_can["SEQUENCE"], self.population_number)),
                        columns=['SEQUENCE', 'PATH', 'NEW_STRING_START', 'COVER'])
                    populationSW = [Text(first_population_df.iloc[i]['SEQUENCE'], first_population_df.iloc[i]['PATH'],
                                         self.first_can["SEQUENCE"], self.first_can["PATH"], self.Z_set_first,
                                         self.P_set_first, self.k_minus_1) for i in range(len(first_population_df))]
                    counter_SW = 0

                if counter_RY > 10:
                    second_population_df = pd.DataFrame(
                        list(create_permutation(self.second_can["SEQUENCE"], self.population_number)),
                        columns=['SEQUENCE', 'PATH', 'NEW_STRING_START', 'COVER'])
                    populationRY = [Text(second_population_df.iloc[i]['SEQUENCE'], second_population_df.iloc[i]['PATH'],
                                         self.second_can["SEQUENCE"], self.second_can["PATH"], self.Z_set_first,
                                         self.P_set_first, self.k_minus_1) for i in range(len(second_population_df))]
                    counter_RY = 0

                i_generation += 1
                counter_SW += 1
                counter_RY += 1
                '''
                for j in range(len(new_populationSW)):
                    new_populationSW[j].check_if_correct()
                for j in range(len(new_populationRY)):
                    new_populationRY[j].check_if_correct()
                '''
    def return_class_string_time(self):
        return self.final_string, self.time


def population_threads(FIRST_POP_DF, SECOND_POP_DF, self):

    # Creates two processes
    FINISH_TIME = self.TIMEOUT

    start_time = time.time()
    queue = mp.Queue()
    items = [(FIRST_POP_DF, self.first_can, self, Bcolors.FAIL, FINISH_TIME, queue, ), (SECOND_POP_DF, self.second_can, self, Bcolors.OKGREEN, FINISH_TIME, queue, )]


    vs = range(2)
    '''
    pool = mp.Pool()
    result = [pool.apply_async(run_thread_population, args = items[v]) for v in vs]
    pool.close()
    pool.join()
    '''
    result = [mp.Process(target=run_thread_population, args = items[v]) for v in vs]
    for p in result:
        p.start()

    #print()
    #print(Bcolors.BOLD + Bcolors.LIGHT_BLUE + "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" + Bcolors.ENDC)
    #print(Bcolors.BOLD + Bcolors.LIGHT_BLUE + "ALGORYTM GENETYCZNY SZUKANIE DLA SOND SW I RY" + Bcolors.ENDC)
    for p in result:
        p.join()
    #print(Bcolors.BOLD + Bcolors.LIGHT_BLUE + "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" + Bcolors.ENDC)

    unsorted_result = [queue.get() for _ in vs]
    final_string = encode(unsorted_result[0], unsorted_result[1])
    end_time = time.time() - start_time
    return final_string, end_time

def init_population(POP_DF, can, self):
    population = [
        Text(POP_DF.iloc[i]['SEQUENCE'], POP_DF.iloc[i]['PATH'], can["SEQUENCE"], can["PATH"], self.Z_set_first,
             self.P_set_first, self.k_minus_1) for i in range(len(POP_DF))]
    population.sort(key=lambda x: x.fitness, reverse=True)
    population_len = len(population)
    temp_best = 0
    temp_length = 10000000
    return population, population_len, temp_best, temp_length


def select_model_SW_RY(populationSW, populationRY, selection_model):
    selectedSW = selection_model(populationSW)
    selectedRY = selection_model(populationRY)
    return selectedSW.copy(), selectedRY.copy()


def while_mutation(populationSW, populationRY, new_populationSW, new_populationRY, cross_number, mutation_probability):
    temp1 = numpy.random.randint(0, len(populationSW) - 1)
    temp2 = numpy.random.randint(0, len(populationSW) - 1)
    childSW = populationSW[temp1].crossover(populationSW[temp2], cross_number)
    childRY = populationRY[temp1].crossover(populationRY[temp2], cross_number)
    new_populationSW = while_mutation_append(childSW, mutation_probability, new_populationSW)
    new_populationRY = while_mutation_append(childRY, mutation_probability, new_populationRY)


def while_mutation_append(child, mutation_probability, new_population):
    if child != None:
        if random() <= mutation_probability:
            child.mutation()
        if not any(x.PATH == child.PATH for x in new_population):
            new_population.append(child)
    return new_population


def if_found_max_cover(best, can, k_minus_1, population_number, Z_set_first, P_set_first, size_sequence):
    can = combine_string(best, can, k_minus_1, size_sequence)
    can = can.reset_index(drop=True)
    if len(can['SEQUENCE']) < 3:
        exit()
    population_df = pd.DataFrame(list(create_permutation(can["SEQUENCE"], population_number)),
                                 columns=['SEQUENCE', 'PATH', 'NEW_STRING_START', 'COVER'])
    populationNEW = [Text(population_df.iloc[i]['SEQUENCE'], population_df.iloc[i]['PATH'],
                          can["SEQUENCE"], can["PATH"], Z_set_first, P_set_first,
                          k_minus_1) for i in range(len(population_df))]
    temp_best = 0
    counter = 0
    return populationNEW, can, temp_best, counter


def show_generation_info(i, the_best_match, fitness, seq_parts_str, color_of_string):
    log(''.join((color_of_string, "-------------------------------------")), DEBUG)
    log(''.join("Generation:               {}".format(i)), DEBUG)
    log(''.join("The best match SEQUENCE:  {}".format(the_best_match)), DEBUG)
    log(''.join("FITNESS:                  {}".format(fitness)), DEBUG)
    log(''.join("Sequence candidates SIZE: {}".format(len(seq_parts_str))), DEBUG)
    log(''.join(("-------------------------------------", Bcolors.ENDC)), DEBUG)


def check_if_str_prt_correct(self, s1df):
    for i in range(len(self.seq_parts_str)):
        temp = self.seq_parts_path
        to_read_path = temp.tolist()
        to_read_str = s1df['SW poll'].tolist()
        seq = to_read_str[to_read_path[i][0]]
        for j in range(len(to_read_path[i]) - 1):
            cover = endOverlap(decode(to_read_str[to_read_path[i][j]]), to_read_str[to_read_path[i][j + 1]])
            seq = ''.join((last_element_decode(seq), to_read_str[to_read_path[i][j + 1]][cover:]))
        if seq != self.seq_parts_str[i]:
            print("NOT OK")
            exit()

def encode(Z_set, P_set):
    result = ""
    if Z_set[0] == "R" or Z_set[0] == "Y":
        temp = Z_set
        Z_set = P_set
        P_set = temp


    for i in range (len(Z_set)):
        if Z_set[i] == "W" and P_set[i] == "R":
            result += "A"
        elif Z_set[i] == "W" and P_set[i] == "Y":
            result += "T"
        elif Z_set[i] == "S" and P_set[i] == "R":
            result += "G"
        elif Z_set[i] == "S" and P_set[i] == "Y":
            result += "C"
        else:
            result += Z_set[i]
    return result
