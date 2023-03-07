import xml.etree.ElementTree as et
import time
import sys

sys.setrecursionlimit(3000)
result_remember = ""

def parse_file(file_path):
    try: tree = et.parse(file_path)
    except IOError:
        print(f"Problem reading file {file_path}")
        raise
    root = tree.getroot()
    start = root.attrib["start"]
    seq_len = int(root.attrib["length"])
    probe1_pattern, probe2_pattern = root[0].attrib["pattern"], root[1].attrib["pattern"]
    s1, s2 = [], []
    for oligon in root[0]: s1.append(oligon.text)
    for oligon in root[1]: s2.append(oligon.text)
    return start, seq_len, probe1_pattern, probe2_pattern, s1, s2

def replace_last_letter(first, second):
    first  = first.replace("A", "W").replace("T", "W").replace("G", "S").replace("C", "S")
    second = second.replace("A", "R").replace("T", "Y").replace("G", "R").replace("C", "Y")
    return first, second

def encode(Z_set, P_set):
    result = ""
    for i in range (len(Z_set)):
        if    Z_set[i] == "S" and P_set[i] == "R": result = ''.join((result, "G"))
        elif  Z_set[i] == "W" and P_set[i] == "R": result = ''.join((result, "A"))
        elif  Z_set[i] == "W" and P_set[i] == "Y": result = ''.join((result, "T"))
        elif  Z_set[i] == "S" and P_set[i] == "Y": result = ''.join((result, "C"))
        else: result += Z_set[i]
    return result

def decode(string):
    first, second = '', ''
    dlugosc = len(string)
    for i in range(dlugosc):
        if i < (dlugosc - 1):
            if   string[i] == 'G': first = ''.join((first, 'S')); second = ''.join((second, 'R'))
            elif string[i] == 'A': first = ''.join((first, 'W')); second = ''.join((second, 'R'))
            elif string[i] == 'T': first = ''.join((first, 'W')); second = ''.join((second, 'Y'))
            elif string[i] == 'C': first = ''.join((first, 'S')); second = ''.join((second, 'Y'))
        else:                      first = ''.join((first, string[i])); second = ''.join((second, string[i]))
    return first, second

def endOverlap(a, b):
    max = 0
    for i in range(len(a)):
        if b.startswith(a[-i:]): max = i
    if max != 0: return max
    return max

def found_element(string_seq, sx, ols):
    levels = {}
    for element in sx:
        pos = sx.index(element)
        if endOverlap(string_seq, element) > 0 and ols[element] == False:
            x = endOverlap(string_seq, element)
            if x in list(levels.keys()):
                temp_dictionary_variable = levels[x]
                temp_dictionary_variable.append([element, pos])
                dictionary_variable = {x: temp_dictionary_variable}
            else:
                dictionary_variable = {x: [[element, pos]]}
            levels.update(dictionary_variable)
    keys = list(levels.keys())
    res, val = [], []
    for i in reversed(sorted(keys)):
        for j in levels[i]:
            res.append(j); val.append(i)
    return res, val, levels
    #res - [string pokrycia, numer w xml]
    #val - [lista posortowana pokryć zaczynając od największego pokrycia]
    #levels - znalezione pokrycia

def search(first, second, ols1, ols2, s1, s2):
    first, second = replace_last_letter(first, second)
    res1, val1, level1 = found_element(first, s1, ols1)
    res2, val2, level2 = found_element(second, s2, ols2)
    if len(level1) == 0 or len(level2) == 0: return (0,0), (0,0), (0,0)
    else: return (first, second), (res1, res2), (val1, val2)

def bin_chip_recursion(end_re, cre_nuc_len, seq_len, set_first, paths, ols, sx, result_found, no_more_solutions):
    if result_found == True: return end_re, cre_nuc_len, seq_len, set_first, paths, ols, sx, result_found, no_more_solutions
    if (cre_nuc_len == seq_len):
        end_re = encode(set_first[0], set_first[1])
        return end_re, cre_nuc_len, seq_len, set_first, paths, ols, sx, True, no_more_solutions
    set_first, gen_str, pos = search(set_first[0], set_first[1], ols[0], ols[1], sx[0], sx[1])
    if len(gen_str[0]) == 0 or len(gen_str[1]) == 0:
        return end_re, cre_nuc_len, seq_len, set_first, paths, ols, result_found, sx, True, no_more_solutions
    if (no_more_solutions == True and cre_nuc_len != seq_len) or cre_nuc_len > seq_len:
        last_wrong_Z_element = paths[0][-1]
        ols[0][str(list(ols[0].keys())[paths[0][-1]])] = False
        to_delete_1 = paths[3].pop()
        temp_len = len(paths[0]) - (len(list(ols[0].keys())[paths[0][-1]]) - to_delete_1)
        paths[0].pop()
        set_first[0] = set_first[0][:temp_len]
        ols[1][str(list(ols[1].keys())[paths[1][-1]])] = False
        to_delete_2 = paths[4].pop()
        temp_len = len(paths[1]) - (len(list(ols[1].keys())[paths[1][-1]]) - to_delete_2)
        paths[1].pop()
        set_first[1] = set_first[1][:temp_len]
        cre_nuc_len -= (len(list(ols[0].keys())[last_wrong_Z_element]) - to_delete_1)
        return end_re, cre_nuc_len, seq_len, set_first, paths, ols, sx, result_found, no_more_solutions
    else:
        return addOutgoingVerticesToList(end_re, cre_nuc_len, seq_len, set_first, paths, ols, sx, gen_str, pos, result_found, False)

def addOutgoingVerticesToList(end_re, cre_nuc_len, end_seq_len, set_first, paths, ols, sx, gen_str, pos, result_found, no_more_solutions):
    for fir_set_l in range(len(gen_str[0])):
        if result_found == True: break
        fir_str = gen_str[0][fir_set_l]
        for sec_set_l in range(len(gen_str[1])):
            if result_found == True: break
            sec_str = gen_str[1][sec_set_l]
            #sprawdzamy 1) czy ostatnie litery sie zgadzaja AND 2) czy poziom pokrycia jest taki sam
            if fir_str[0][-1] == sec_str[0][-1] and pos[0][fir_set_l] == pos[1][sec_set_l]:
                temp_len = cre_nuc_len + (len(fir_str[0]) - pos[0][fir_set_l])
                if end_seq_len >= temp_len:
                    cre_nuc_len += (len(fir_str[0]) - pos[0][fir_set_l])
                    ols[0][fir_str[0]], ols[1][sec_str[0]] = True, True
                    set_first = (set_first[0] + fir_str[0][pos[0][fir_set_l]:], set_first[1] + sec_str[0][pos[1][sec_set_l]:])
                    paths[0].append(fir_str[1])
                    paths[1].append(sec_str[1])
                    paths[2].append(pos[0][fir_set_l])
                    paths[3].append(pos[1][sec_set_l])
                    end_re, cre_nuc_len, end_seq_len, set_first, paths, ols, sx, result_found, no_more_solutions \
                        = bin_chip_recursion(end_re, cre_nuc_len, end_seq_len, set_first, paths, ols, sx, result_found, no_more_solutions)
    return bin_chip_recursion(end_re, cre_nuc_len, end_seq_len, set_first, paths, ols, sx, result_found, True)

def fun_iter(size_sequence):
    start_time                                             = time.time()
    start, seq_len, probe1_pattern, probe2_pattern, s1, s2 = parse_file("sample_cases/positive_numbers_max_k_max_sqt_max_pos/" + str(size_sequence) + ".xml")

    Z_set_first, P_set_first                               = decode(start)
    ols1, ols2                                             = {x: False for x in s1}, {x: False for x in s2}
    ols1[Z_set_first], ols2[P_set_first]                   = True, True
    path_Z, path_P, depth_Z, depth_P                       = [list(ols1).index(Z_set_first)], [list(ols2).index(P_set_first)], [0], [0]
    created_nuc_length                                     = len(Z_set_first)
    #********* PACK VALUES *****************************************************
    set_first                                              = (Z_set_first, P_set_first)
    ols                                                    = (ols1, ols2)
    paths                                                  = (path_Z, path_P, depth_Z, depth_P)
    sx                                                     = (s1, s2)
    pack_arg                                               = (created_nuc_length, seq_len, set_first, paths, ols, sx)
    no_more_solutions, result_found                        = False, False

    result = bin_chip_recursion('', pack_arg[0], pack_arg[1], pack_arg[2], pack_arg[3], pack_arg[4], pack_arg[5], result_found, no_more_solutions)

    end_time = time.time() - start_time
    #print(time.time() - start_time)
    #print(result[0])
    return result[0], end_time

#--------------------------
#UNCOMMENT THIS IF YOU WANT RUN THIS ALGHORITHM

if __name__ == '__main__':
    result, end_time = fun_iter(400)
    print(result)
    print(end_time)
