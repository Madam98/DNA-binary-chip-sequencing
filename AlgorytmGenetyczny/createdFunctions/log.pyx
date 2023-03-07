#!python
#cython: language_level=3
import time
DEBUG = False               #<---- wyswietla podstawowe informacje dzialania programu
DEBUG_CROSSOVER = True    #<---- informacje zwiazane z krzyzowaniem populacji (NIE DZIAŁA!)
DEBUG_GENERATION = True    #<---- wyswietla informacje na temat nowo wygenerowanych generacji
st = time.time()

def log(s, DEBUG):
    if DEBUG:
        print (s)

def end_time():
    print("CZAS:                 " + str(time.time() - st))