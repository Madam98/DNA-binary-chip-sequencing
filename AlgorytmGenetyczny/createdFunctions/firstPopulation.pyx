#!python
#cython: language_level=3
import pandas as pd
from createdClasses.colors             import Bcolors
from createdFunctions.stringOperations import find_same_last_letter, \
                                              last_element_decode, \
                                              endOverlap, \
                                              decode
from createdFunctions.log              import log, DEBUG, st, end_time
import time
import numpy
from random import randint, sample
pd.set_option("display.max_colwidth", None)
'''
WARUNKIEM ZNALEZIENIA MAKSYMALNEGO POKRYCIA NIECH BEDZIE 
ZNALEZIENIE MAKSYMALNEGO POKRYCIA Z TA SAMA KONCOWA LITERA 
DLA OBU SOND
'''
def first_sequence_k_minus_1_for_both_poll(s1df, s2df, ols1, ols2, k_minus_1, population_number, cross_number, size_sequence):
    log(''.join((Bcolors.FAIL, 'Ustalam największe pokrycie równe k-1 dla\nobu sond z taką samą ostatnią literą...\n', Bcolors.ENDC)), DEBUG)
    gen_path_1, gen_path_2, end_path_1, end_path_2, end_string_1, end_string_2 = [],[],[],[],[],[]
    if max(len(s1df), len(s2df)) == len(s1df):
        for i in range(len(s1df)):
            check, i2 = find_same_last_letter(s1df._get_value(i, 'SW poll'), s2df, ols2, 'RY poll')
            if check == True and ols1[s1df._get_value(i, 'SW poll')] == False:
                gen_path_1.append(i)
                gen_path_2.append(i2)
                saved_path_1, saved_path_2, saved_string_1, saved_string_2 = look_for_indexes_max(i, i2, s1df, s2df, ols1, ols2, s1df._get_value(i, 'SW poll'), s2df._get_value(i2, 'RY poll'), k_minus_1, gen_path_1, gen_path_2, [], [], [], [])
                append_end_lists(saved_path_1, saved_path_2, saved_string_1, saved_string_2, end_path_1, end_path_2, end_string_1, end_string_2)
    else:
        for i in range(len(s2df)):
            check, i2 = find_same_last_letter(s2df._get_value(i, 'RY poll'), s1df, ols1, 'SW poll')
            if check == True and ols2[s2df._get_value(i, 'RY poll')] == False:
                gen_path_1.append(i2)
                gen_path_2.append(i)
                saved_path_1, saved_path_2, saved_string_1, saved_string_2 = look_for_indexes_max(i, i2, s2df, s1df, ols2, ols1, s2df._get_value(i, 'RY poll'), s1df._get_value(i2, 'SW poll'), k_minus_1, gen_path_2, gen_path_1, [], [], [], [])
                append_end_lists(saved_path_2, saved_path_1, saved_string_2, saved_string_1, end_path_1, end_path_2, end_string_1, end_string_2)
    [add_not_used_seq(ols1, s1df, end_path_1, end_string_1, 'SW poll', i) if ols1[s1df._get_value(i, 'SW poll')] == False else '' for i in range(len(s1df))]
    [add_not_used_seq(ols2, s2df, end_path_2, end_string_2, 'RY poll', i) if ols2[s2df._get_value(i, 'RY poll')] == False else '' for i in range(len(s2df))]
    print_generated_parts_info(end_string_1, end_string_2, end_path_1, end_path_2, DEBUG)
    # USE IF YOU WANT TO CHECK IF GENERATED CANDIDATES IS DONE RIGHT
    check_candidates(end_path_1, end_path_2, end_string_1, end_string_2, k_minus_1)
    log('\n'.join((Bcolors.OKBLUE, "Wyznaczanie losowych permutacji\nWyznaczam pierwszą przykładową populację...", Bcolors.ENDC)), DEBUG)

    check_end_string_1, check_end_string_2 = False, False
    if (len(end_string_1) == 1):
        print("Wielkosc wygenerowanych kandydatow rowna 1")
        #print(end_string_1[0])
        check_end_string_1 = True

    if (len(end_string_2) == 1):
        print("Wielkosc wygenerowanych kandydatow rowna 1")
        #print(end_string_2[0])
        check_end_string_2 = True



    if check_end_string_1:
        first_population_df = pd.DataFrame(list(zip(end_string_1, end_path_1, [0], [0])), columns=['SEQUENCE', 'PATH', 'NEW_STRING_START', 'COVER'])
    else:
        first_population_df = pd.DataFrame(list(create_permutation(end_string_1, population_number, cross_number)), columns=['SEQUENCE', 'PATH', 'NEW_STRING_START', 'COVER'])
    first_candidates_df = pd.DataFrame(list(zip(end_string_1, end_path_1)), columns=['SEQUENCE', 'PATH'])
    first_candidates_df.set_index('SEQUENCE')

    if check_end_string_2:
        second_population_df = pd.DataFrame(list(zip(end_string_2, end_path_2, [0], [0])),
                                        columns=['SEQUENCE', 'PATH', 'NEW_STRING_START', 'COVER'])
    else:
        second_population_df = pd.DataFrame(list(create_permutation(end_string_2, population_number, cross_number)), columns=['SEQUENCE', 'PATH', 'NEW_STRING_START', 'COVER'])

    #print(first_population_df)
    #print(second_population_df)
    second_candidates_df = pd.DataFrame(list(zip(end_string_2, end_path_2)), columns=['SEQUENCE', 'PATH'])
    second_candidates_df.set_index('SEQUENCE')


    log(''.join((Bcolors.OKBLUE, "Wyznaczyłem pierwszą populację!\n\n\n", Bcolors.ENDC)), DEBUG)

    return first_population_df, second_population_df, first_candidates_df, second_candidates_df, cross_number

def look_for_indexes_max(i1, i2, s1df, s2df, ols1, ols2, gen_string_1, gen_string_2, k_minus_1, created_path_1, created_path_2, saved_path_1, saved_path_2, saved_string_1, saved_string_2):
    found_indexes_1, found_indexes_2 = [], []
    head1, head2 = s1df.head().columns[0], s2df.head().columns[0]
    t_string = last_element_decode(s1df._get_value(i1, head1))
    found_s1_element, found_s2_element, index_1, index_2 = find_next_pair(i1, i2, s1df, s2df, ols1, ols2, k_minus_1, head1, head2, t_string)
    found_indexes_1.append(index_1)
    found_indexes_2.append(index_2)
    if index_2 == None:
        if len(gen_string_1) == (k_minus_1 + 1):
            ols1[s1df._get_value(i1, head1)] = True
            ols2[s2df._get_value(i2, head2)] = True
            saved_path_1, saved_path_2 = append_path(created_path_1, created_path_2, saved_path_1, saved_path_2)
            list_clear(created_path_1, created_path_2)
            return saved_path_1, saved_path_2, gen_string_1, gen_string_2
        elif ols1[s1df._get_value(i1, head1)] == True and ols2[s2df._get_value(i2, head2)] == True and len(gen_string_1) != (k_minus_1 + 1):
            saved_path_1, saved_path_2 = append_path(created_path_1, created_path_2, saved_path_1, saved_path_2)
            list_clear(created_path_1, created_path_2)
            return saved_path_1, saved_path_2, gen_string_1, gen_string_2
    else:
        if len(found_indexes_1) > 1:
            print("ZNALEZIONO PODWOJNE DOPASOWANIE MAKSYMALNE!\n KONIEC PROGRAMU")
            exit()
        else:
            created_path_1, created_path_2 = append_path(found_indexes_1, found_indexes_2, created_path_1, created_path_2)
            set_all_on_true(s1df, s2df, ols1, ols2, i1, i2, found_indexes_1, found_indexes_2, head1, head2)
            gen_string_1 = ''.join((last_element_decode(gen_string_1), s1df._get_value(found_indexes_1[0], head1)[-1]))
            gen_string_2 = ''.join((last_element_decode(gen_string_2), s2df._get_value(found_indexes_2[0], head2)[-1]))
            sav_1, sav_2, gen_string_1, gen_string_2 = look_for_indexes_max(index_1, index_2, s1df, s2df, ols1, ols2, gen_string_1, gen_string_2, k_minus_1, created_path_1, created_path_2, saved_path_1, saved_path_2, saved_string_1, saved_string_2)
            return sav_1, sav_2, gen_string_1, gen_string_2

def create_permutation(gen_str_list, population_number, cross_number):
    result, path_list, list_new_seq_start, list_cover = [], [], [], []

    for i in range(population_number):

        #print("CREATE PERMUATATION")
        #print(cross_number)
        #print(len(gen_str_list))
        #print((range(len(gen_str_list)))[:numpy.random.randint(cross_number + 1, len(gen_str_list))])
        range_gen = (range(len(gen_str_list)))
        #print(range_gen)
        #print(numpy.random.randint(cross_number + 1, size=len(gen_str_list)))

        generated_path = list(numpy.random.permutation(range_gen)[:numpy.random.randint(cross_number, high=len(gen_str_list))])
        #print(len(gen_str_list))
        #print((range(len(gen_str_list)))[:numpy.random.randint(cross_number + 1, len(gen_str_list))])
        #print(generated_path)
        #print()
        new_seq_start = [0]
        cover_to_add = []
        if not generated_path in path_list:
            path_list.append(generated_path)
            seq = decode(gen_str_list[generated_path[0]])
            for j in range(len(generated_path) - 1):
                cre_new_string = decode(gen_str_list[generated_path[j]])
                cover = endOverlap(cre_new_string, gen_str_list[generated_path[j + 1]])
                new_seq_start.append(len(seq))
                seq = ''.join((last_element_decode(seq), gen_str_list[generated_path[j + 1]][cover:]))
                cover_to_add.append(cover)
            result.append(seq)
            new_seq_start.append(len(seq))
            list_cover.append(cover_to_add)
            list_new_seq_start.append(new_seq_start)
        else:
            i -= 1

    return zip(result, path_list, list_new_seq_start, list_cover)

def find_next_pair(i1, i2, s1df, s2df, ols1, ols2, k_minus_1, head1, head2, t_string):
    for i in range(len(s1df)):
        if s1df._get_value(i1, head1) != s1df._get_value(i, head1) and ols1[s1df._get_value(i, head1)] == False:
            if endOverlap(t_string, s1df._get_value(i, head1)) == k_minus_1:
                for j in range(len(s2df)):
                    if s2df._get_value(i2, head2) != s2df._get_value(j, head2) and ols2[s2df._get_value(j, head2)] == False:
                        if endOverlap(last_element_decode(s2df._get_value(i2, head2)), s2df._get_value(j, head2)) == k_minus_1:
                            return s1df._get_value(i, head1), s2df._get_value(j, head2), i, j
    return None, None, None, None

def set_all_on_true(s1df, s2df, ols1, ols2, i1, i2, found_indexes_1, found_indexes_2, head1, head2):
    ols1[s1df._get_value(i1, head1)], ols2[s2df._get_value(i2, head2)], ols1[s1df._get_value(found_indexes_1[0], head1)], ols2[s2df._get_value(found_indexes_2[0], head2)] = True, True, True, True

def list_clear(created_path_1, created_path_2):
    created_path_1 *= 0  # <-- fast clear
    created_path_2 *= 0

def append_path(created_path_1, created_path_2, saved_path_1, saved_path_2):
    if len(created_path_2) == 1 : saved_path_1.append(created_path_1[0]); saved_path_2.append(created_path_2[0])
    else: saved_path_1, saved_path_2 = created_path_1.copy(), created_path_2.copy()
    return saved_path_1, saved_path_2

def append_end_lists(saved_path_1, saved_path_2, saved_string_1, saved_string_2, end_path_1, end_path_2, end_string_1, end_string_2):
    end_path_1.append(saved_path_1); end_path_2.append(saved_path_2); end_string_1.append(saved_string_1); end_string_2.append(saved_string_2)

def print_generated_parts_info(end_string_1, end_string_2, end_path_1, end_path_2, DEBUG):
    log(''.join((Bcolors.FAIL, "Uzyskane elementy które posłużą do budowania następnych populacji", Bcolors.ENDC)), DEBUG)
    log(''.join(("Dla sondy SW:")), DEBUG)
    log('\n'.join([''.join((str(i), ":\t", end_string_1[i])) for i in range(len(end_string_1))]), DEBUG)
    log(end_path_1, DEBUG)
    log(''.join(("Dla sondy RY:")), DEBUG)
    log('\n'.join([''.join((str(i), ":\t", end_string_2[i])) for i in range(len(end_string_2))]), DEBUG)
    log(end_path_2, DEBUG)

def add_not_used_seq(ols1, s1df, end_path_1, end_string_1, name, i):
    ols1[s1df._get_value(i, name)] = True; end_path_1.append([i]); end_string_1.append(s1df._get_value(i, name))

def check_candidates(end_path_1, end_path_2, end_string_1, end_string_2, k_minus_1):
    log(''.join((Bcolors.FAIL, "Sprawdzam czy stworzeni kandydaci zostali poprawnie wygenerowani...", Bcolors.ENDC)), DEBUG)
    flat_list_1, flat_list_2 = numpy.sort([x for xs in end_path_1 for x in xs]), numpy.sort([x for xs in end_path_2 for x in xs])
    check_list = [i for i in range(len(flat_list_1)) if not i in flat_list_1]
    if check_list != []:
        print(''.join((Bcolors.FAIL, "ZOSTAŁ NIEPRZYDZIELONY KANDYDAT!\nKOŃCZĘ PRACĘ PROGRAMU!", Bcolors.ENDC))); exit()
    sum1, sum2 = 0, 0
    for i in range(len(end_string_1)): sum1 = sum1 + 1 + (len(end_string_1[i]) - (k_minus_1 + 1))
    for i in range(len(end_string_2)): sum2 = sum2 + 1 + (len(end_string_2[i]) - (k_minus_1 + 1))
    if len(flat_list_1) == sum1 and len(flat_list_2) == sum2: log(''.join((Bcolors.OKGREEN, "OK", Bcolors.ENDC)), DEBUG)
    else:
        log(''.join((Bcolors.FAIL, "Uzyskany różny wynik dla kandydatów, kończenie programu...", Bcolors.ENDC))); exit()

def combine_string(max_cover, can_df, k_minus_1, size_sequence):
        temp_list = []
        for i in reversed(range(len(max_cover))):
            left = max_cover[i][0]
            right = max_cover[i][1]
            temp_list.append(right)
            can_df.iloc[left]["PATH"].extend(can_df.iloc[right]["PATH"])
            string_len = len(can_df.iloc[left]["SEQUENCE"])
            string_len_left = len(can_df.iloc[left]["SEQUENCE"])
            string_len_right = len(can_df.iloc[right]["SEQUENCE"])
            can_df.iloc[left]["SEQUENCE"] = can_df.iloc[left]["SEQUENCE"][0:string_len - k_minus_1] + can_df.iloc[right]["SEQUENCE"]

            if (string_len_left - k_minus_1) + (string_len_right - k_minus_1) + k_minus_1 == len(
                    can_df.iloc[left]["SEQUENCE"]):
                if len(can_df.iloc[left]["SEQUENCE"]) == size_sequence:
                    return can_df.iloc[left]["SEQUENCE"]
            else:
                print("KRYTYCZNY BLAD PODCZAS ŁĄCZENIA SEKWENCJI!!!\nKONIEC PROGRAMU")
                exit()
        can_df = can_df.drop(temp_list)
        return can_df

def strong_number_check(cross_number, population_number):
    Temp, i, n = cross_number + 1, 1, 1
    while i <= Temp:
        n = n * i
        if population_number < n: return True
        i += 1
    return False