#!python
#cython: language_level=3

from createdClasses.element import Element
from createdFunctions.stringOperations import last_element_decode, endOverlap, decode
from createdFunctions.log              import log, DEBUG, DEBUG_CROSSOVER, st, end_time
import pandas as pd
import numpy
import time
from random import randint, sample
from createdClasses.colors import *

class Text(Element):
    #def __init__(self, population_df_str, population_df_path, population_df_new_string_start, population_df_cover, seq_parts_str, seq_parts_path, Z_set_first, P_set_first, k_minus_1):
    def __init__(self, population_df_str, population_df_path, seq_parts_str, seq_parts_path, Z_set_first, P_set_first, k_minus_1):
        self.SEQ = population_df_str
        self.PATH = population_df_path
        #self.NEW_STRING_START = population_df_new_string_start
        #self.COVER = population_df_cover
        self.seq_parts_str = seq_parts_str
        self.seq_parts_path = seq_parts_path
        self.Z_set_first = Z_set_first
        self.P_set_first = P_set_first
        self.k_minus_1 = k_minus_1
        self.charsSW = ('S', 'W')
        self.charsRY = ('R', 'Y')
        self.evaluate_value = 0
        #self.found_max_cover = []
        super().__init__()

    '''
    Do zmiennej PATH bierzemy jej INDEKSY nie wartosci! self.PATH.index[0]
    '''
    def evaluate_function(self):
        counter, result, created_seq = 0, 0, last_element_decode(self.seq_parts_str[self.PATH[0]])
        if self.SEQ.startswith(self.charsSW): start_check = self.Z_set_first
        else: start_check = self.P_set_first
        found_max_cover = []
        for i in range(len(self.PATH) - 1):
            max_cover = endOverlap(last_element_decode(self.seq_parts_str[self.PATH[i]]), self.seq_parts_str[self.PATH[i + 1]])
            if max_cover == self.k_minus_1:
                found_max_cover.append((self.PATH[i], self.PATH[i + 1]))
                result = result + 20
            result = result + max_cover
            created_seq = ''.join((created_seq, last_element_decode(self.seq_parts_str[self.PATH[i + 1]])[max_cover:]))
            counter += 1
        self.evaluate_value = result
        if found_max_cover:
            self.found_max_cover = found_max_cover
        return self.evaluate_value

    def crossover(self, element2: 'Element', number_of_crossover) -> 'Element':
        first_list, second_list = self.PATH, element2.PATH

        #print("ŚCIEŻKI PRZED CROSSOVER")
        #print(self.PATH)
        #print(element2.PATH)
        #print()

        to_cut_list_first, to_cut_list_second = cut(number_of_crossover, first_list, second_list)
        if True:
            i, result = 0, []
            while i != number_of_crossover + 1:
                how_to_combine = define_combine(to_cut_list_first, to_cut_list_second, i)
                if how_to_combine == 0:
                    if i == number_of_crossover: result.append(first_list[to_cut_list_first[i]:])
                    else: result.append(first_list[to_cut_list_first[i]:to_cut_list_first[i + 1]])
                else:
                    if i == number_of_crossover: result.append(second_list[to_cut_list_second[i]:])
                    else: result.append(second_list[to_cut_list_second[i]:to_cut_list_second[i + 1]])
                i += 1


            flat_list = [x for xs in result for x in xs]

            #print("ŚCIEŻKA PO CROSSOVER")
            #print(flat_list)
            #exit()
            flat_list, counter, yes_duplicate = find_duplicate(flat_list)
            list_example = numpy.arange(0, len(self.PATH), 1)
            result_diffrence = list(set(list_example) - set(flat_list))
            if counter > len(result_diffrence): return None
            gen_sample = sample(result_diffrence, len(result_diffrence))
            for i in range(len(flat_list)):
                if flat_list[i] == -1: flat_list[i] = gen_sample.pop()
            seq = last_element_decode(self.seq_parts_str[flat_list[0]])
            for j in range(len(flat_list) - 1):
                cover = endOverlap(decode(self.seq_parts_str[flat_list[j]]), self.seq_parts_str[flat_list[j + 1]])
                seq = ''.join((last_element_decode(seq), self.seq_parts_str[flat_list[j + 1]][cover:]))
        else:
            #pierwszy wyraz podstawienie
            path, new_string_list, cover_list, i,  = [], [], [], 0
            how_to_combine = define_combine(to_cut_list_first, to_cut_list_second, i)
            #how_to_combine = 1
            if how_to_combine == 0:
                if i == number_of_crossover: path.append(first_list[to_cut_list_first[i]:])
                else: path.append(first_list[to_cut_list_first[i]:to_cut_list_first[i + 1]])
                seq = self.SEQ[0:self.NEW_STRING_START[len(path[0])]]
            else:
                if i == number_of_crossover: path.append(second_list[to_cut_list_second[i]:])
                else: path.append(second_list[to_cut_list_second[i]:to_cut_list_second[i + 1]])
                seq = element2.SEQ[0:element2.NEW_STRING_START[len(path[0])]]
            crossover_show_info_first(how_to_combine, self.seq_parts_str, self.SEQ, element2.SEQ, seq, self.PATH, element2.PATH, to_cut_list_first, to_cut_list_second, self.NEW_STRING_START, element2.NEW_STRING_START)

            to_add_first = to_cut_list_first[i + 1]
            to_add_second = to_cut_list_second[i + 1]

            if how_to_combine == 0:
                new_string_list.append(self.NEW_STRING_START[to_cut_list_first[0]:to_cut_list_first[1]])
                cover_list.append(self.COVER[to_cut_list_first[0]:to_cut_list_first[1] - 1])
            else:
                new_string_list.append(element2.NEW_STRING_START[to_cut_list_second[0]:to_cut_list_second[1]])
                cover_list.append(element2.COVER[to_cut_list_second[0]:to_cut_list_second[1] - 1])
            i = 1
            #-------------------------------------------------------------------------------------
            while i != number_of_crossover + 1:
                previous = how_to_combine
                how_to_combine = define_combine(to_cut_list_first, to_cut_list_second, i)
                #how_to_combine = 1
                cross_len_path_first, cross_len_path_second,  offset_first, offset_second = set_start_while_cross(to_cut_list_first, to_cut_list_second, self, element2, i)
                crossover_show_info(how_to_combine, self.SEQ, element2.SEQ, offset_first, offset_second)
                to_add_first = to_cut_list_first[i + 1]
                to_add_second = to_cut_list_second[i + 1]
                if how_to_combine == 0:
                    if i == number_of_crossover:
                        path.append(first_list[to_cut_list_first[i]:])
                        new_string_list.append(self.NEW_STRING_START[to_cut_list_first[i]:])
                    else:
                        path.append(first_list[to_cut_list_first[i]:to_cut_list_first[i + 1]])
                        new_string_list.append(self.NEW_STRING_START[to_cut_list_first[i]:to_cut_list_first[i + 1]])
                    cover_list.append(self.COVER[to_cut_list_first[i]:to_cut_list_first[i + 1] - 1])
                else:
                    if i == number_of_crossover:
                        path.append(second_list[to_cut_list_second[i]:])
                        new_string_list.append(element2.NEW_STRING_START[to_cut_list_second[i]:])
                    else:
                        path.append(second_list[to_cut_list_second[i]:to_cut_list_second[i + 1]])
                        new_string_list.append(element2.NEW_STRING_START[to_cut_list_second[i]:to_cut_list_second[i + 1]])
                    cover_list.append(element2.COVER[to_cut_list_second[i]:to_cut_list_second[i + 1] - 1])

                if how_to_combine == 0:
                    cover_cut_list = to_cut_list_first[i] - 1
                    new_part_string_begin = self.NEW_STRING_START[cross_len_path_first]
                    new_part_string_end   = self.NEW_STRING_START[to_add_first]
                    created_string = self.seq_parts_str[first_list[to_cut_list_first[i]]][:self.COVER[cover_cut_list]] + self.SEQ[new_part_string_begin:new_part_string_end]
                else:
                    cover_cut_list = to_cut_list_second[i] - 1
                    new_part_string_begin = element2.NEW_STRING_START[cross_len_path_second]
                    new_part_string_end   = element2.NEW_STRING_START[to_add_second]
                    created_string = element2.seq_parts_str[second_list[to_cut_list_second[i]]][:element2.COVER[cover_cut_list]] + element2.SEQ[new_part_string_begin:new_part_string_end]
                cover = endOverlap(seq, created_string)
                cover_list[-1].insert(0, cover)

                while_crossover_show_info(previous, path, seq, self.SEQ, self.NEW_STRING_START, self.COVER, element2.SEQ, created_string, self.seq_parts_str, offset_first, offset_second, to_add_first, first_list, to_cut_list_first, i, cover_cut_list, cover, element2.COVER)
                new_string_list[i][0] = len(seq)
                seq = ''.join((seq, created_string[cover:]))
                cross_len_path_first = to_add_first
                cross_len_path_second = to_add_second
                i += 1

            flat_list            = [x for xs in path for x in xs]
            flat_list_new_string = [x for xs in new_string_list for x in xs]
            flat_list_cover      = [x for xs in cover_list for x in xs]
            print(flat_list)
            print(flat_list_new_string)
            print(flat_list_cover)
            flat_list, counter, yes_duplicate = find_duplicate(flat_list)
            if yes_duplicate:
                print(flat_list, counter)
                list_example = numpy.arange(0, len(self.PATH), 1)
                result_diffrence = list(set(list_example) - set(flat_list))
                if counter > len(result_diffrence): return None
                gen_sample = sample(result_diffrence, len(result_diffrence))
                for i in range(len(flat_list)):
                    if flat_list[i] == -1:
                        print(i)
                        print()
                        seq_dup_1 = seq[0:flat_list_new_string[i]]
                        seq_dup_2 = seq[flat_list_new_string[i]:]
                        to_change = gen_sample.pop()
                        flat_list[i] = to_change
                        print(to_change)
                        print(seq)
                        print(seq_dup_1)
                        print(seq_dup_2)
                        print(flat_list)
                        #musimy zmodyfikowac z wielu stron
                        cover = endOverlap(decode(seq_dup_1), decode(self.seq_parts_str[flat_list[i]]))
                        print(self.seq_parts_str[flat_list[i]])
                        print(cover)
                        #print(self.seq_parts_str[flat_list[i]][cover:])
                        seq = ''.join((seq_dup_1, self.seq_parts_str[flat_list[i]][cover:]))
                        seq = last_element_decode(seq)
                        print(flat_list_cover)
                        print(flat_list_new_string)
                        flat_list_cover[i - 1] = cover
                        flat_list_new_string[i + 1] = len(seq)
                        print(flat_list_cover)
                        print(flat_list_new_string)
                        #cover = endOverlap(decode(seq_dup), seq[flat_list_new_string[i + 1]:flat_list_new_string[len(flat_list) - 1]]))
                        #seq = ''.join((seq_dup, self.seq_parts_str[flat_list[i + 1]][cover:]))

                        exit()
        tt1 = self.seq_parts_str.to_frame('SEQUENCE')
        tt2 = self.seq_parts_path.to_frame('PATH')
        return Text(seq, flat_list, tt1["SEQUENCE"], tt2["PATH"], self.Z_set_first, self.P_set_first, self.k_minus_1)

    def _perform_mutation(self):
        temp_list = numpy.arange(0, len(self.seq_parts_path), 1)
        list_diffrences = list(set(temp_list) - set(self.PATH))
        if len(list_diffrences) > 0:
                temp_random = sample(list_diffrences, 1)
                score, created_seq, created_path = [], [], []
                old_path, old_string, old_fitness = self.PATH.copy(), self.SEQ, self.fitness,
                for i in range(len(old_path)):
                    temp_path = old_path.copy()
                    temp_path[i] = temp_random[0]
                    seq = last_element_decode(self.seq_parts_str[temp_path[0]])
                    for j in range(len(temp_path) - 1):
                        cover = endOverlap(decode(self.seq_parts_str[temp_path[j]]), self.seq_parts_str[temp_path[j + 1]])
                        seq = ''.join((last_element_decode(seq), self.seq_parts_str[temp_path[j + 1]][cover:]))
                    self.PATH = temp_path
                    self.SEQ = seq
                    self.fitness = self.evaluate_function()
                    created_path.append(temp_path)
                    created_seq.append(seq)
                    score.append(self.fitness)
                max_value = max(score)
                if max_value > old_fitness:
                    max_index = score.index(max_value)
                    self.seq_path = created_path[max_index]
                    self.seq_string = created_seq[max_index]
                    self.fitness = score[max_index]
                else:
                    self.seq_path = old_path
                    self.seq_string = old_string
                    self.fitness = old_fitness

    def check_if_correct(self):
        seq = last_element_decode(self.seq_parts_str[self.PATH[0]])
        cover_list = []
        for j in range(len(self.PATH) - 1):
            #print(j)
            cover = endOverlap(decode(self.seq_parts_str[self.PATH[j]]), self.seq_parts_str[self.PATH[j + 1]])
            cover_list.append(cover)
            seq = ''.join((last_element_decode(seq), self.seq_parts_str[self.PATH[j + 1]][cover:]))
        if seq!= self.SEQ:
            print("seq:         " + seq)
            print("cover:       " + str(cover_list))
            print("self.SEQ:    " + self.SEQ)
            print("self.COVER:  " + str(self.COVER))
            print("SEKWENCJA NIE ROWNA!")
            print()
            print(self.PATH)
            exit()

    def get_text(self):
        return self.PATH

    def __repr__(self):
        return self.SEQ

    def get_length(self):
        return len(self.SEQ)

    def get_max_covers(self):
        return self.found_max_cover

    def return_path_length(self):
        return len(self.seq_parts_path)



def define_combine(to_cut_list_first, to_cut_list_second, i):
    if i <= max(len(to_cut_list_first), len(to_cut_list_second)):
        if i < len(to_cut_list_first) and i < len(to_cut_list_second): how_to_combine = randint(0, 1) # obie listy nie osiagnely maksa, losujemy ktory nastepny element bierzemy
        elif i >= len(to_cut_list_first) and i < len(to_cut_list_second): how_to_combine = 1
        elif i < len(to_cut_list_first) and i >= len(to_cut_list_second): how_to_combine = 0
    return how_to_combine

def cut(number_of_crossover, first_list, second_list):
    if number_of_crossover < len(first_list):
        to_cut_list_first = sample([*range(1, len(first_list), 1)], number_of_crossover)
        to_cut_list_first.sort()
    else:
        to_cut_list_first = [*range(0, len(first_list), 1)]
    if number_of_crossover < len(second_list):
        to_cut_list_second = sample([*range(1, len(second_list), 1)], number_of_crossover)
        to_cut_list_second.sort()
    else:
        to_cut_list_second = [*range(0, len(second_list), 1)]
    to_cut_list_first.insert(0, 0)
    to_cut_list_first.append(len(first_list))
    to_cut_list_second.insert(0, 0)
    to_cut_list_second.append(len(second_list))
    return to_cut_list_first, to_cut_list_second

def find_duplicate(flat_list):
    counter, seen = 0, set()
    yes_duplicate = False
    for i, e in enumerate(flat_list):
        if e in seen:
            flat_list[i] = -1
            counter += 1
            yes_duplicate = True
        else:
            seen.add(e)
    return flat_list, counter, yes_duplicate

def crossover_show_info_first(how_to_combine, seq_parts_str, first_SEQ, second_SEQ, seq, first_PATH, second_PATH, to_cut_list_first, to_cut_list_second, first_NEW_STRING, second_NEW_STRING):
    log(''.join('\n\n--------------------------------------------------------------------------------------'), DEBUG_CROSSOVER)
    #log('\n'.join((seq_parts_str)), DEBUG)
    log(''.join(("HOT_TO_COMBINE: ", str(how_to_combine))), DEBUG_CROSSOVER)
    log(''.join((Bcolors.OKGREEN, "\n" "first.SEQ:   ", first_SEQ, Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.FAIL, "second.SEQ:  ", second_SEQ, Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join(("seq:         ", seq)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.OKGREEN, "\n", "first.PATH:  ", str(first_PATH), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.FAIL, "second.PATH: ", str(second_PATH), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.OKGREEN, "\n", "to_cut_list_first:  ", str(to_cut_list_first), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.FAIL, "to_cut_list_second: ", str(to_cut_list_second), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.OKGREEN, "\n", "cross_len_path_first:     ", str(to_cut_list_first[1]), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.FAIL, "cross_len_path_second:    ", str(to_cut_list_second[1]), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.OKGREEN, "\n", "first.NEW_STRING_START:  ", str(first_NEW_STRING), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.FAIL, "second.NEW_STRING_START: ", str(second_NEW_STRING), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join('--------------------------------------------------------------------------------------'), DEBUG_CROSSOVER)

def crossover_show_info(how_to_combine, first_SEQ, second_SEQ, offset_first, offset_second):
    log(''.join('\n\n*************************************************************************************'), DEBUG_CROSSOVER)
    log(''.join(("HOW_TO_COMBINE: ", str(how_to_combine))), DEBUG_CROSSOVER)
    log(''.join((Bcolors.OKGREEN, "\n", first_SEQ, Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.FAIL, second_SEQ, Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.OKGREEN, "\n", "offset_first:  ", str(offset_first), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join((Bcolors.FAIL, "offset_second: ", str(offset_second), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join('*************************************************************************************'), DEBUG_CROSSOVER)

def while_crossover_show_info(previous, path, seq, first_SEQ, first_NEW_STRING_START, first_COVER, second_SEQ, created_string, self_seq_parts_str, offset_first, offset_second, to_add_first, first_list, to_cut_list_first, i, te, cover, second_COVER):
    log(''.join(("\n\n.....................................................................................")), DEBUG_CROSSOVER)
    log(''.join(("| PATH:              ", str(path))), DEBUG_CROSSOVER)
    log(''.join(("| Before cross SEQ:  ", str(seq))), DEBUG_CROSSOVER)
    log(''.join(("|", Bcolors.FAIL, " FIRST SEQ:         ", first_SEQ, Bcolors.ENDC)), DEBUG_CROSSOVER)
    printer = ''
    if previous == 0: printer += ' ' * offset_first
    else:             printer += ' ' * offset_second
    log(''.join(("| ", Bcolors.OKGREEN, "SECOND SEQ:        ", printer, first_SEQ[offset_first:first_NEW_STRING_START[to_add_first]], Bcolors.ENDC)), DEBUG_CROSSOVER)
    printer = ''
    #print(previous)
    #print(te)
    #print(len(first_COVER))
    #print(len(second_COVER))
    '''
    if previous == 0:
        if te < first_COVER:  offset = offset_first - first_COVER[te]
        else:
            print("aa")
    else:
        if te < second_COVER: offset = offset_second - second_COVER[te]
        else:
            print("bb")
    '''
    printer += ''# * offset
    log(''.join(("| ", Bcolors.OKBLUE, "END SEQ:           ", printer, created_string, Bcolors.ENDC)), DEBUG_CROSSOVER)
    #log(''.join((self_seq_parts_str[first_list[to_cut_list_first[i]]])), DEBUG)
    #log(''.join((str(first_list[to_cut_list_first[i]]))), DEBUG)
    if te < len(first_COVER):  log(''.join(("| ", Bcolors.FAIL, "COVER FIRST:       ", str(first_COVER[te]), Bcolors.ENDC)), DEBUG_CROSSOVER)
    if te < len(second_COVER): log(''.join(("| ", Bcolors.OKGREEN, "COVER SECOND:      ", str(second_COVER[te]), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join(("| ", Bcolors.OKBLUE, "COVER SEQ-CROSS:   ", str(cover), Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join(("| ", Bcolors.FAIL, "FIRST SEQ:         ", first_SEQ, Bcolors.ENDC)), DEBUG_CROSSOVER)
    log(''.join(("| ", Bcolors.OKGREEN, "SECOND SEQ:        ", second_SEQ, Bcolors.ENDC)), DEBUG_CROSSOVER)
    seq = ''.join((seq, created_string[cover:]))
    log(''.join(("| ", "After cross SEQ:   ", str(seq))), DEBUG_CROSSOVER)
    log(''.join((".....................................................................................")), DEBUG_CROSSOVER)

def set_start_while_cross(to_cut_list_first, to_cut_list_second, self, element2, i):
    return to_cut_list_first[i], to_cut_list_second[i], self.NEW_STRING_START[to_cut_list_first[i]], element2.NEW_STRING_START[to_cut_list_second[i]]