#!python
#cython: language_level=3

def find_same_last_letter(letter, second_list, second_dict, name):
    for i in range(len(second_list)):
        if letter[-1] == second_list._get_value(i, name)[-1] and second_dict[second_list._get_value(i, name)] == False: return True, i
    return False, None

def duplicate(x):
    _size = len(x)
    repeated = []
    for i in range(_size):
        k = i + 1
        for j in range(k, _size):
            if x[i] == x[j] and x[i] not in repeated: repeated.append(x[i])
    return repeated

def endOverlap(a, b):   # <--https://stackoverflow.com/questions/1315559/how-good-is-startswith operacja zabierala nam najwiecej czasu, trzeba polepszyc evaluate
    max = 0
    for i in range(len(a)):
        if b.startswith(a[-i:]): max = i
    if max != 0: return max
    return max

def last_element_decode(seq_string):
    if seq_string[0] == "W" or seq_string[0] == "S": return seq_string.replace("A", "W").replace("T", "W").replace("G", "S").replace("C", "S")
    else: return seq_string.replace("A", "R").replace("T", "Y").replace("G", "R").replace("C", "Y")

def decode(string, both = False):
    first, second = '', ''
    letter = string[0]
    for i in range(len(string)):
        if string[i] == 'G':   first += 'S'; second += 'R'
        elif string[i] == 'A': first += 'W'; second += 'R'
        elif string[i] == 'T': first += 'W'; second += 'Y'
        elif string[i] == 'C': first += 'S'; second += 'Y'
        else:
            first += string[i]
            second += string[i]
    if both == True: return first, second
    elif letter == 'S' or letter == 'W': return first
    else: return second

def scale(element_of_list, number_of_permutations, list_k_max):
    result_index = element_of_list[0]
    result = list_k_max[result_index]
    for i in range(1, len(element_of_list)):
        cover = endOverlap(result, list_k_max[i])
        result = result[:len(result) - cover] + str(list_k_max[element_of_list[i]])
    return result

def merge(events, actors, index=None):
    if len(events) == 0: return []
    if index is None: index = set()
    merged = None
    tried_pairs = set()
    for event in events:
        for i, actor in enumerate(actors):
            pair = (event, actor)
            if pair not in index and pair not in tried_pairs:
                new_index = index.union([pair])
                rest = merge(events[1:], actors[:i] + actors[i + 1:], new_index)
                if rest is not None:
                    merged = [pair] + rest
                    break
                else: tried_pairs.add(pair)
        if merged is not None: break
    return merged