#!python
#cython: language_level=3

import xml.etree.ElementTree as Et
import pandas

def parse_file(file_path):
    try:
        tree = Et.parse(file_path)
    except IOError:
        print(f"Problem reading file {file_path}")
        raise
    root = tree.getroot()
    start = root.attrib["start"]
    seq_len = int(root.attrib["length"])
    probe1_pattern = root[0].attrib["pattern"]
    probe2_pattern = root[1].attrib["pattern"]
    s1 = []
    s2 = []
    for oligon in root[0]:
        s1.append(oligon.text)
    for oligon in root[1]:
        s2.append(oligon.text)

    ols1 = {x: False for x in s1}
    ols2 = {x: False for x in s2}

    s1df = pandas.DataFrame(s1, columns=['SW poll'])
    s2df = pandas.DataFrame(s2, columns=['RY poll'])

    return start, seq_len, probe1_pattern, probe2_pattern, ols1, ols2, s1df, s2df, s1, s2