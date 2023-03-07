#!python
#cython: language_level=3
from createdClasses.text import Text
from createdFunctions.log import end_time, DEBUG, DEBUG_GENERATION, st, log
from createdFunctions.firstPopulation import create_permutation, combine_string
import threading
from multiprocessing import Pool, Queue
import time
import numpy
from random import random
import pandas as pd
from createdClasses.colors              import Bcolors


pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_rows', None)

def run_thread_population(POP_DF, can, self, color, TIMEOUT, queue):
    time_population = time.time()

    i_generation, counter, population, population_len, best_match, best_length, best_can = init_population(POP_DF, can, self)

    while True:
        if time.time() - time_population > TIMEOUT:
            print(Bcolors.BOLD + Bcolors.BRIGHT_MAGNETA + "|-----------------------|" + Bcolors.ENDC)
            print(Bcolors.BOLD + Bcolors.BRIGHT_MAGNETA + "| TIMEOUT PRZEKROCZONY! |\n| KONCZE PRACE          |" + Bcolors.ENDC)
            print(Bcolors.BOLD + Bcolors.BRIGHT_MAGNETA + "|-----------------------|" + Bcolors.ENDC)
            print(color + "********************")
            end_time()

            print(color + "ZNALEZIONA SEKWENCJA: " + str(max(new_population, key=lambda x: x.fitness)))
            end_time()
            print("********************" + Bcolors.ENDC)
            queue(max(new_population, key=lambda x: x.fitness))
            return max(new_population, key=lambda x: x.fitness)

        new_population, self.cross_number, population_len = select_model_SW_RY(population, self.selection_model, population_len, self)
        while len(new_population) != population_len:
            while_mutation(population, new_population, self.cross_number, self.mutation_probability)
        best_match = max(new_population, key=lambda x: x.fitness)

        if best_match.found_max_cover:
            population, can, best_match, counter = if_found_max_cover(best_match.found_max_cover, can, self.k_minus_1, self.population_number, self.Z_set_first, self.P_set_first, self.size_sequence, self.cross_number, color)
            if population == True:
                queue.put(can)
                return
        elif best_length > len(can):

            show_generation_info(i_generation, best_match, best_match.fitness, best_match.seq_parts_str, color)
            best_length = len(can)
            if best_match.fitness > best_length:
                population = new_population

        if counter > 10:
            population_df = pd.DataFrame(
                list(create_permutation(can["SEQUENCE"], self.population_number, self.cross_number)),
                columns=['SEQUENCE', 'PATH', 'NEW_STRING_START', 'COVER'])
            population = [Text(population_df.iloc[i]['SEQUENCE'], population_df.iloc[i]['PATH'],
                                 can["SEQUENCE"], can["PATH"], self.Z_set_first,
                                 self.P_set_first, self.k_minus_1) for i in range(len(population_df))]
            counter = 0
        i_generation += 1
        counter += 1

def init_population(POP_DF, can, self):
    population = [
        Text(POP_DF.iloc[i]['SEQUENCE'], POP_DF.iloc[i]['PATH'], can["SEQUENCE"], can["PATH"], self.Z_set_first,
             self.P_set_first, self.k_minus_1) for i in range(len(POP_DF))]
    population.sort(key=lambda x: x.fitness, reverse=True)
    population_len = len(population)
    i, counter, temp_best, best_length, best_can = 0, 0, 1, 10000000000, 1000000000
    return i, counter, population, population_len, temp_best, best_length, best_can

def select_model_SW_RY(population, selection_model, population_len, self):
    #strong_check, number_strong = strong_number_check(self.cross_number, self.population_number)
    #if not strong_check:
        #self.cross_number = 1
        #population_len = number_strong / 2
        #population.sort(key=lambda x: x.fitness, reverse=True)
        #selected = [population[0], population[1]].copy()
    #else:
    selected = selection_model(population, self.selection_percent)
    return selected.copy(), self.cross_number, population_len

def while_mutation(population, new_population, cross_number, mutation_probability):
    temp1 = numpy.random.randint(0, len(population) - 1)
    temp2 = numpy.random.randint(0, len(population) - 1)
    child = population[temp1].crossover(population[temp2], cross_number)
    new_population = while_mutation_append(child, mutation_probability, new_population)

def while_mutation_append(child, mutation_probability, new_population):
    if child != None:
        if random() <= mutation_probability:
            child.mutation()
        if not any(x.PATH == child.PATH for x in new_population):
            new_population.append(child)
    return new_population

def if_found_max_cover(best_match_found_max_cover, can, k_minus_1, population_number, Z_set_first, P_set_first, size_sequence, cross_number, color):
    can = combine_string(best_match_found_max_cover, can, k_minus_1, size_sequence)

    if type(can) is str:
        print(color + "********************")
        print("ZNALEZIONA SEKWENCJA: " + str(can))
        end_time()
        print("********************"  + Bcolors.ENDC)
        return True, can, 1, 0

    can = can.reset_index(drop=True)

    #print(can)
    #print(len(can))
    #print(population_number)
    #print(cross_number)


    population_df = pd.DataFrame(list(create_permutation(can["SEQUENCE"], population_number, cross_number)),
                                 columns=['SEQUENCE', 'PATH', 'NEW_STRING_START', 'COVER'])
    population = [Text(population_df.iloc[i]['SEQUENCE'], population_df.iloc[i]['PATH'],
                       can["SEQUENCE"], can["PATH"], Z_set_first, P_set_first,
                       k_minus_1) for i in range(len(population_df))]
    temp_best = 1
    counter = 0
    return population, can, temp_best, counter

def show_generation_info(i, the_best_match, fitness, seq_parts_str, color_of_string):
    log(''.join((color_of_string, "-------------------------------------")), DEBUG_GENERATION)
    log(''.join("Generation:               {}".format(i)), DEBUG_GENERATION)
    log(''.join("The best match SEQUENCE:  {}".format(the_best_match)), DEBUG_GENERATION)
    log(''.join("FITNESS:                  {}".format(fitness)), DEBUG_GENERATION)
    log(''.join("Sequence LENGTH:          {}".format(the_best_match.get_length())), DEBUG_GENERATION)
    log(''.join("Sequence candidates SIZE: {}".format(len(seq_parts_str))), DEBUG_GENERATION)
    log(''.join(("-------------------------------------", Bcolors.ENDC)), DEBUG_GENERATION)

'''
def strong_number_check(cross_number, population_number):
    Temp, i, n = cross_number + 1, 1, 1
    while i <= Temp:
        n = n * i
        if population_number < n: return True, 0
        i += 1
    return False, n
'''
