# DNA-binary-chip-sequencingLaboratorium z bioinformatyki

Sekwencjonowanie DNA w oparciu o chipy binarne z błędami pozytywnymi Prowadzący laboratoria: **dr inż. Piotr Wawrzyniak**

Danylo Brushko Adam Miernicki

143137 135750

Numer grupy dziekańskiej: L13 Numer grupy dziekańskiej: L13     Termin zajęć: wtorek 8:00 Termin zajęć: wtorek 8:00         danylo.brushko@student.put.poznan.pl adam.miernicki@student.put.poznan.pl

Termin wymagany: 22-06-2022 Termin oddania: 15-09-2022

**Spis treści**

[**1 Opis problemu](#_page1_x50.40_y36.00) **1 [2 Definicja problemu](#_page1_x50.40_y132.56) 1**

[**3 Przedstawienie działania algorytmu dokładnego](#_page1_x50.40_y414.47) **1**

1. [Zbiór danych wstępnych .](#_page1_x50.40_y504.03) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 1
1. [Środowisko implementacji danego algorytmu . .](#_page1_x50.40_y646.50) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 1
1. [Działanie algorytmu ](#_page2_x50.40_y36.00). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 2
1. [Dekodowanie początkowej sekwencji ](#_page2_x50.40_y223.69). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 2
1. [Główna funkcja poszukiwań . ](#_page3_x50.40_y36.00). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 3
1. [Znajdywanie potencjalnych następników - search . ](#_page3_x50.40_y393.88). . . . . . . . . . . . . . . . . . . . . . . . . 3
4. [Szacowana złożoność obliczeniowa ](#_page4_x50.40_y36.00). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4
5. [Testy .](#_page4_x50.40_y116.65) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4
1. [Metodologia wykonywanych testów . ](#_page4_x50.40_y136.87). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4
1. [Wyniki ](#_page4_x50.40_y204.67). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4

[**4 Przedstawienie algorytmu genetycznego](#_page5_x50.40_y36.00) **5**

1. [Zbiór danych wstępnych .](#_page5_x50.40_y58.62) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 5
1. [Środowisko implementacji danego algorytmu . .](#_page5_x50.40_y200.24) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 5
1. [Opis działania algorytmu ](#_page5_x50.40_y281.83). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 5
1. [Wyznaczanie początkowych kandydatów . .](#_page7_x50.40_y36.00) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 7
1. [Permutacja powstałych kandydatów w celu zbudowania populacji . . .](#_page7_x50.40_y556.04) . . . . . . . . . . . . . . 7
1. [Sposób oceny osobników populacji oraz sposób selekcji . ](#_page8_x50.40_y355.01). . . . . . . . . . . . . . . . . . . . . . 8
1. [Sposób krzyżowania się osobników z populacji . ](#_page8_x50.40_y638.40). . . . . . . . . . . . . . . . . . . . . . . . . . . 8
1. [Omówienie błędów implementacyjnych oraz sugestia polepszenia algorytmu dla metody crossover 9](#_page9_x50.40_y334.14)
1. [Działanie mutacji . ](#_page9_x50.40_y581.66). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 9
1. [Omówienie błędów implementacyjnych oraz sugestia polepszenia algorytmu dla metody mutation 9](#_page9_x50.40_y697.68)
1. [Ulepszanie kandydatów po wygenerowaniu nowej populacji . ](#_page10_x50.40_y36.00). . . . . . . . . . . . . . . . . . . 10
1. [Warunek stopu .](#_page10_x50.40_y116.10) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 10
4. [Omówienie dodatkowych zależności występujących w kodzie . .](#_page10_x50.40_y388.43) . . . . . . . . . . . . . . . . . . . . . . 10
5. [Omówienie występujących błędów krytycznych w algorytmie . . ](#_page10_x50.40_y496.42). . . . . . . . . . . . . . . . . . . . . . 10
5. [Szacowana złożoność czasowa ](#_page11_x50.40_y241.17). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 11
5. [Testy .](#_page11_x50.40_y359.18) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 11
1. [Metodologia wykonywanych testów . ](#_page11_x50.40_y379.95). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 11
1. [Wyniki ](#_page12_x50.40_y36.00). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 12

[**5 Oznaczenie zmiennych w algorytmie oraz sposób uruchamiania programu ](#_page13_x50.40_y36.00)**13 [6 Podsumowanie wyników oraz omówienie pracy działania algorytmów ](#_page13_x50.40_y337.79)13**

0

1  **Opis problemu**

W niniejszym sprawozdaniu zostanie zaprezentowany algorytm rozwiąujący problem sekwencjonowania przez hybry- dyzację oparty o nieklasyczny chip DNA - chipie binarnym. W definicji problemu zostanie opisana zasada konstrukcji rozpatrywanego chipu oraz sposób działania algorytmu rozwiązujący problem SBH algorytmem dokładnym oraz algo- rytmem genetycznym. Przedstawiony zostanie także wynik eksperymentów obliczeniowych, w których testowane będą omawiane algorytmy.

2  **Definicja problemu**

Rozwiązując dany problem sekwencjonowania DNA w oparciu o chipy binarne z błędami pozytywnymi, warto na początku wyjaśnić kilka podstawowych pojęć składających się na definicjędanego zagadnienia. Proces sekwencjonowa-

nia składa się z dwóch faz: biochemicznej oraz obliczeniowej. Podczas procesu biochemicznego przeprowadzany jest eksperyment hybrydyzacyjny z pełną biblioteką oligonukletydów o pewnej ustalonej długości *n*, nazywanych *n*-merami. W wyniku przeprowadzonego eksperymentu otrzymujemy zbiór, który określany jest mianem spektrum zawierającym wszystkie podciągi o długości *n* możliwe do wyróżnienia w analizowanej sekwencji DNA. W procesie sekwencjonowa- nia wyróżniamy dwa możliwe do wystąpienia niedoskonałości, spowodowane procesem hybrydyzacji. W sytuacji, gdy

w naszym spektrum nie występuje n-merów, mimo występowania ich w badanej sekwencji, mówimy o występowaniu błędów negatywnych. Natomiast, w przypadku występowania błędów pozytywnych możemy w wynikowym spektrum wyznaczyć *n*-merowe sekwencje, które przed dokonaniem procesu sekwencjonowania nie występują.

Sam chip binarny DNA będzie wykorzystywał inny alfabet znaków do zapisu sond niż alfabet stosowany w podejściu klasycznym A, C, G, T. Jest to alfabet 8-literowy W, S, R, Y, A, C, G, T, który kolejno tworzy dwa 6-literowe alfabety W, S, A, C, G, T oraz R, Y, A, C, G, T. Chip składa się z dwóch sond i to w każdym z nich będzie stosowany jeden rodzaj alfabetu: *SZ* dla alfabetu W, S, A, C, G, T oraz *SP* dla alfabetu R, Y, A, C, G, T. Warto pamię- tać, że wszystkie litery, oprócz ostatniej, w sekwencji kodowane są za pomocą liter W, S, lub R, Y w zależności od sond. Każda kombinacja dwóch liter W, S lub R, Y odpowiada za kodowanie jednej litery A, C, G, T. Zastosowanie dwóch rodzajów sond umożliwia lepszą kontrolę błędów w procesie hybrydyzacji - proces przetwarzania i znajdowa- nia poszukiwanej sekwencji będzie wykorzystywał dwie równocześnie budowane ścieżki, które wzajemnie będą miały możliwość sprawdzania czy budowana sekwencja jest prawidłowa.

3  **Przedstawienie działania algorytmu dokładnego**

W celu rozwiązania przedstawionego problemu został napisany algorytm rozwiązujący dany problem sekwencjonowania w sposób dokładny. Oznacza to, że dany algorytm skupia się na znalezieniu potencjalnych wszystkich rozwiązań oraz że dla większych instancji danych wejściowych może okazać się algorytmem mocno kosztownym czasowo. Najważniejsze części funkcjonowania algorytmu zostaną przedstawione poniżej.

1. **Zbiór danych wstępnych**

Początkowe dane zostały pobrane ze strony [https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php ](https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php)Algorytm został przetestowany odpowiednio dla sond w których wynikowa długość sekwencji wynosiła odpowiednio n = 30,

40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500. Jeśli chodzi o dodatkowe konfigurowalne parametry znajdujące się na wspomnianej wcześniej stronie to zostały one ustawione odpowiednio na:

- k = 10
- sqpe ustawiony na maksymalną zawsze wartość równą zawsze n/4
- pose ustawiony na maksymalną zawsze wartość równą zawsze n/2
2. **Środowisko implementacji danego algorytmu**

Algorytm został stworzony w języku Python w wersji 3.9 z wykorzystaniem IDE Pycharm Community Edition.

3. **Działanie algorytmu**

W celu znalezienia pasujących sekwencji, warto przyjrzeć sie jak przygotowane są nasze dane wejściowe. Dla dwóch rodzajów sond *SZ* oraz *SP* zostały wygenerowane ciągi o długości *k* W, S oraz R, Y. W każdym ciągu ostatnia litera została zapisana za pomocą alfabetu A, G, C, T. Głównym zadaniem algorytmu będzie znajdywanie listy potencjalnych kandydatów posiadających największe pokrycie dla tworzonej sekwencji. W przypadku znalezienia takich kandydatów dla dwóch rodzajów sond, dokonujemy sprawdzenia, czy wybrani kandydaci dla obu sond posiadają takie same n- liczne pokrycie oraz czy są zakończone tym samym nukleotydem. W przypadku stwierdzenia poprawności algorytm akceptuje danych kandydatów z każdej z sond do budowanych ścieżek i zaczyna poszukiwania nastęnych możliwych kandydatów.

Jeżeli otrzymana scieżka będzie dluższa od oczekiwanej długości, to ostatni dodany następnik zostanie usunięty, a następnie będzie szukany następny kandydat, aż do znalezienia szukanej n-długiej sekwencji dla obu sond. Algorytm kończy się, gdy cała przestrzeń możliwych rozwiązań zostanie sprawdzona i program nie będzie w stanie znaleźć następnych możliwych kandydatów.

1. **Dekodowanie początkowej sekwencji**

Zgodnie z definicją działania chipów binarnych, początkowa sekwencja DNA (oznaczona w wygenerowanym pliku danych przez słowo *start* ) została zdekodowana do elementów zbiorów *Sz* oraz *Sp* kodujących odpowiednio sekwencję startową. Następnie dla dwóch wygenerowanych sekwencji startowych szukany był ich odpowiednik w danych zbiorach

w celu znalezienia pierwszych sekwencji rozpoczynajacych kodowanie.

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.001.png)

Rysunek 1: Przykład zdekodowania starowej sekwencji i znalezienia ich odpowiedników w sondzie *SZ* oraz *SP*

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.002.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.003.png)

2. **Główna funkcja poszukiwań**

Główną funkcją odpowiedzialną za znalezienie szukanych sekwencji jest **bin chip recursion**. Jest to funkcja rekuren- cyjna, której warunkiem stopu jest znalezienie szukanej sekwencji lub zakończenie programu przez brak możliwych do znalezienia i dopasowania następnych kandydatów.

Jeśli warunki stopu nie są spełnione, funkcja przystępuje do znalezienia następnych potencjalnych możliwych maksy- malnych dopasowanych sekwencji. Za znalezienie następnych możliwych kandydatów odpowiada funkcja **search**. Jeśli potencjalni kandydaci zostali znalezieni uruchamiana jest funkcja **addOutgoingVerticesToList**, która odpowiada

za dodawanie znalezionych kandydatów do listy rozwiązania i budowania tekstowo danej sekwencji.

**Algorytm 1:** Główna pętla poszukiwań![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.004.png)![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.005.png)

**1 bin chip recursion** *(first, second, created by algorithm sequence length)***: 2 if** *created by algorithm sequence length == expectedLength* **then![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.006.png)**

**3 return** create~~ by~~ algorithm~~ sequence~~ length

**4** *SZ* candidate, *SP* candidate, level ← *search*(S*Z parameters*, *SP parameters*);

**5 if** *return from add == True or created by algorithm sequence length > searched sequence length* **then 6** delete~~ last~~ sequence~~ for~~ both~~ polls

**7** addOutgoingVerticesToList

W głównej pętli poszukiwań następników dla naszej sekwencji warunkami stopu pętli jest:

- znalezienie sekwencji równej długości zmiennej *expectedLength* zadeklarowanej w pliku xml,
- brak możliwości znalezienia potencjalnych następników (funkcja *foundSolutions* ) poszukiwanej sekwencji dla sondy *SZ* oraz *Sp*

W przypadku spełnienia wszystkich warunków tworzona zostaje sekwencja poprzez dołączenie znalezionych sekwencji nukletydów. Głębokość nakładania się nukleotydów deteminuje zmienna *level* określana przez funkcję *findSolution*.

3. **Znajdywanie potencjalnych następników - search**

W celu znalezienia następnych potencjalnych kandydatów możliwych do przyłączenia dla danej sekwencji, funkcja search dekoduje ostatnie litery stworzonej dotychczas sekwencji, odpowiednio dla sond *SZ* oraz *SP* . Szukanie potenc- jalnych dopasowań rozpoczynamy od najwyższego możliwego dopasowania sekwencji. Oznacza to, że dla *k*-długiej sekwencji największe możliwe dopasowanie jest równe *k* - 1. W naszym pseudokodzie poziom dopasowania został oz- naczony przez zmienną *level*. Następnie utworzone słowniki *ols1* oraz *ols2* przekazują informację czy dane sekwencje ze zbiorów *s1* oraz *s2* były wykorzystane. Jeśli znalezieni kandydaci nie zostali wykorzystani, zostają oni dodani jako potencjalni możliwi kandydaci do dołączenia do budowanych sekwencji przez **bin chip recursion**.

**Algorytm 2:** Znajdywanie potencjalnych następników![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.007.png)![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.008.png)

**1 Search** *(firstString, secondString, ols*1*, ols*2*, s1, s2)***:**

**2** firstString, secondString ←transform~~ last~~ letters(*firstString, secondString*)![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.009.png)

**3** found~~ seq~~ and~~ index~~ list~~ 1, level~~ all~~ seq~~ list~~ 1, level1 ← found~~ element(first, s1, ols1)

**4** found~~ seq~~ and~~ index~~ list~~ 2, level~~ all~~ seq~~ list~~ 2, level2 ← found~~ element(second, s2, ols2) **5 if** *level1 == 0 or level2 == 0* **then**

**6 return** no~~ found~~ solution

**7 else**

**8 return** first, second, res1, res2, val1, val2

Algorytm kolejno poszukuje wszystkich możliwych następników dla obu rodzajów sond. W przypadku znalezienia podobnych sekwencji, z takim samym ostatnim nukleotydem, zwracana jest informacja o znalezieniu odpowiadającego kandydata.

Warto zauważyć, że algorytm został stworzony tak, że nie znajduje on tylko jednego maksymalnego kandydata. Algorytm przygotowuje listę wykorzystująć wszystkie możliwe sekwencje które pozostały, układając ją od najgorszych do najlepszych dopasowań. W efekcie czas trwania danego algorytmu jest dłuższy, w porównaniu do zastosowania wspomnianej wcześniej innej metody, za to jednak algorytm jest bardziej uniwersalny i mógłby z spowodzeniem po kilklu przekształceniach być zastosowany do np. znajdywania błędów negatywnych albo znajdywania pokryć w sekwencjach gdzie nie zawsze dane pokrycie musi być maksymalne.

4. **Szacowana złożoność obliczeniowa**

Algorytm w trakcie działania w najbardziej pesymistycznym przypadku do znalezienia każdego możliwego maksymal- nego dopasowania będzie musiał przejrzeć maksymalną ilość sekwencji która nie została wykorzystana. Oznacza to że złożoność czasowa będzie równa O(*nSW* ! + *mRY* !) (gdzie liczby *n* oraz *m* oznaczają liczność zbiorów sond *SZ* oraz *SP* ).

5. **Testy**
1. **Metodologia wykonywanych testów**

W testach wzieła udział wygenerowana wcześniej grupa sekwencji wspomniana w **3.1 Zbiór danych wstępnych**. W celu wyeliminowania rozbieżności wyników dla każdej sekwencji zostało wykonanych 30 testów i jako wynik został podany ich wynik arytmetyczny.

2. **Wyniki**

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.010.png)

Rysunek 2: Osiągnięte wyniki dla algorytmu dokładnego

Możemy zauważyć, że dla bardzo małych sekwencji algorytm działa bardzo szybko. Jednakże wraz ze wzrostem populacji zgodnie z szacowaną złożonością obliczeniową czas rośnie i dla n = 400 wynosi on już prawie 30 sekund. Możemy zauważyć, że dla coraz co większych sekwencji efektywność czasowa mocno spada, dlatego też ten algorytm ten nie powinien być stosowany dla dużych sekwencji.

4  **Przedstawienie algorytmu genetycznego**
1. **Zbiór danych wstępnych**

Początkowe dane zostały pobrane ze strony [https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php ](https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php)Algorytm został przetestowany odpowiednio dla sond w których wynikowa długość sekwencji wynosiła odpowiednio n = 30,

40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500. Jeśli chodzi o dodatkowe konfigurowalne parametry znajdujące się na wspomnianej wcześniej stronie to zostały one ustawione odpowiednio na:

- k = 10
- sqpe ustawiony na maksymalną zawsze wartość równą zawsze n/4
- pose ustawiony na maksymalną zawsze wartość równą zawsze n/2
2. **Środowisko implementacji danego algorytmu**

Algorytm został stworzony w języku Python w wersji 3.9. W celu minimalnego polepszenia rezultatów został wyko- rzystany dodatkowy interpreter Cython w celu translacji powstałego kodu w języku Python na język C. Czas pomiędzy wersją Python oraz Cython w wykonanych testach był marginalny co oznacza w rezultacie, że osiągniete rezultaty przez algorytm genetyczny będą porównywane z algorytmem dokładnym.

3. **Opis działania algorytmu**

Na samym początku opisu zaimplementowanego algorytmu trzeba uwzględnić fakt, że algorytm posiada pewien poważny błąd rzutujący na poprawność działania już na poziomie logicznym (autorzy nie przewidzieli możliwość występowania pewnej sytuacji; szczegółowy opis danego błędu znajduje się w rozdziale **”Omówienie występują- cych błędów w zaimplementowanym algorytmie”**.

W przesłanym programie za wykonywanie algorytmu genetycznego odpowiada klasa **GeneticAlgorithm**, posi- adająca jedną metodę odpowiedzialną za wykonanie algorytmu o nazwie *run()*. Warto zauważyć, że klasa nie ma określonej samej w sobie deklaracji funkcji selekcji, warunku stopu oraz generowania pierwszych kandydatów. Zamiast tego klasa przyjmuje zewnętrznie zadeklarowanej funkcje zdefiniowane przez twórcę kodu. Poprzez słowo ’kandydat’

w dalszej części opisu algorytmu będziemy rozumieć wygenerowane sekwencje składające się z maksymalnych pokryć, służące do budowania nastęnych populacji.

Działanie algorytmu możemy opisać za pomocą kilku występujących po sobie faz które się zapętlają. Są to:

1. **Faza generowania kandydatów** (występuje tylko raz na początku działania algorytmu)
1. **Budowanie populacji opartej o losową ilość permutacji na podstawie wygenerowanych kandydatów**
1. **Podanie wygenerowanej populacji selekcji w celu wytworzenia nowej populacji**
1. **Krzyżowanie wyselekcjonowanych osobników ze starej populacji w celu wytworzenia nowej**
1. **Potencjalne wystąpienie mutacji na wytworzonym nowym osobniku**
1. **Jeśli w nowej populacji wygenerowany osobnik posiada sekwencję, która składa się wszędzie z maksymalnych pokryć to występuje komunikat o znalezionej sekwencji i zakończenie programu**
1. **Jeśli w nowej populacji znaleziono maksymalne pokrycie lub też pokrycia to zastąpienie starej listy kandydatów, nowszą ze scalonymi kandydatami posiadającymi maksymalne pokrycie. Powrót do punktu drugiego**

W porównaniu do algorytmu dokładnego nie występuje tutaj faza pobierania i dekodowania startowej sekwencji (oznaczona w wygenerowanym pliku .xml przez słowo *start* ). W przesłanym pliku klasa *GenethicAlghorithm* co prawda przyjmuje pole o nazwie *self.start*, jednak nie jest ono wykorzystywane w finalnej wersji algorytmu.

**Algorytm 3:** GeneticAlgorithm![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.004.png)![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.005.png)

**1 GeneticAlgorithm:**

**2** SWSetPop, RYSetPop, SWGenCan, RYGenCan, PopLen ← first~~ generation~~ func ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.011.png)**3** counterSW = 0 counterRY = 0 **while** *True* **do**

**4** newPopulationSW, newPopulationRY ← select~~ model~~ SW~~ RY

**5**

**6** ###*generate for SW population*

**7 while** *newPopLen != PopLen* **do**

**8** child ← SWSetPop[randPopMember1].crossover(SWSetPop[randPopMember2]) **9 if** *child != None* **then**

**10 if** *genRandomNumber*   *mutationProbability* **then**

**11** child.mutation

**12 if** *not child.StringSeqSW in newPopulationSW.StringSeqSW* **then 13** newPopulationSW.append(child)

**14**

**15** ###*generate for RY population*

**16 while** *newPopLen != PopLen* **do**

**17** child ← RYSetPop[randPopMember1].crossover(RYSetPop[randPopMember2]) **18 if** *child != None* **then**

**19 if** *genRandomNumber*   *mutationProbability* **then**

**20** child.mutation

**21 if** *not child.StringSeqRY in newPopulationRY.StringSeqRY* **then 22** newPopulationRY.append(child)

**23**

**24** BestMatchSW ← max(newPopulationSW.child.fitness) **25** BestMatchRY ← max(newPopulationRY.child.fitness)

**26**

**27 if** *counterSW > 10* **then**

**28** SWSetPop = generate~~ new~~ populationSW

**29** counterSW = 0

**30 else**

**31**

**32** ###*if found max cover improve candidates and generated new population based on new created candidates* **33 if** *BestMatchSW.found max cover* **then**

**34** PopulationSW, SWGenCan← ImproveCanSW~~ AND~~ GenNewPopSW counterSW = 0

**35 else**

**36** SWSetPop ← newPopulationSW

**37 if** *counterRY > 10* **then**

**38** RYSetPop = generate~~ new~~ populationRY counterRY = 0

**39 else**

**40**

**41** ###*if found max cover improve candidates and generated new population based on new created candidates* **42 if** *BestMatchRY.found max cover* **then**

**43** PopulationRY, RYGenCan← ImproveCanRY~~ AND~~ GenNewPopRY

**44** counterRY = 0

**45 else**

**46** RYSetPop ← newPopulationRY

1. **Wyznaczanie początkowych kandydatów**

Pod pojęciem wyznaczania kandydatów rozumiemy wygenerowanie takich sekwencji z których następnie będziemy mogli ’budować’ nasze pokolenia. To czym musi się odznaczać kandydat to posiadanie maksymalnych pokryć w swoim stworzonym ciągu. W innym przypadku łączenie pokoleń posiadających łańcuchy nie mające maksymalnych pokryć nigdy nie doprowadziłoby nas do znalezienia szukanej końcowej sekwencji. Przykład wygenerowanych pierwszych kandydatów dla obydwu sond znajduje się poniżej.

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.012.png)

Wygenerowani kandydaci dla sekwencji n=30 dla sondy SW oraz RY

Możemy więc zauważyć, że pierwszy wygenerowany kandydat SSWSWSWSWWSWSSSSSSSSSSA składa się z części sond znajdujących się w pliku *.xml* po kolei z linii: 0, 2, 9 oraz 10. Po wygenerowaniu danych kandydatów funkcja sprawdzając **check candidates** sprawdza czy wszystkie wygenerowane sekwencje są długości: (*k*-1) + *n* (gdzie *k* oznacza długość jednej sekwencji w sondzie a *n* ilość wszystkich sekwencji). W przypadku stwierdzenia poprawności wygenerowanych sekwencji algorytm przystępuje do następnej fazy działania - tworzenia pierwszej populacji.

Samo wytworzenie początkowych kandydatów jest podobne na tym etapie do działania algorytmu dokładnego. Program iteracyjnie bierze sekwencje z pliku .xml i szuka sewkencji które będą miały maksymalne pokrycie z pobraną sekwencją. Czyli dla przykładu dla 10 elementowego zbioru dla pierwszej sewkencji zostanie zbadanych maksymal- nie 9 kandydatów, w przypadku znalezienia maksymalnego pokrycia sekwencje są ze sobą łączone i dla następnej sprawdzanej sekwencji istniejących kandydatów jest już tylko ośmiu. Podobny schemat działania występował w al- gorytmie dokładnym z tą różnicą, że program ten nie jest wywoływany rekurencyjnie, aż do znalezienia końcowej sekwencji tylko raz do wygenerowania pierwszych kandydatów.

Implementacja danej funkcjonalności znajduje się w pliku **firstPopulation.pyx** pod nazwą

**first sequence k minus 1 for both poll**.

2. **Permutacja powstałych kandydatów w celu zbudowania populacji**

Po wygenerowaniu pierwszych kandydatów następuje faza losowego permutowania powstałych kandydatów. Od razu trzeba podkreślić fakt, że nie są tworzone wszystkie możliwe permutacjea a jedynie jakaś ich część, gdyż wygenerowanie wszystkich możliwych permutacji jest niepotrzebne i bardzo czasochłone obliczeniowo. Dla sekwencji o długości 500 korzystając z zależności, że długość końcowej sekwencji jest równa: *n* + *l* (gdzie *n* oznacza długości sekwencji, a *l* ilość pokryć), to dla sekwencji o długości np. 15 potrzeba by było wygenerować 485! permutacji. Przy uwzględnieniu faktu,

że mamy do czynienia z błędami pozytywnymi i posiadamy nadmiarowe niepasujące sekwencje końcowa ilość permu-

tacji wyniosłaby (485 + *m*)! (gdzie *m* oznacza liczbę dodatkowych błędnych sekwencji). W celu optymalizacji pracy algorytmu ilość permutacji jest ograniczona do ilości docelowej populacji określonej przez **GenethicAlghorithm**. Je- dynym wyjątkiem od tej reguły jest fakt wystąpienia sytuacji gdy ilość permutacji z wygenerowanych początkowych kandydatów jest mniejsza niż wielkość populacji określona przez klasę **GenethicAlghoritm** (np. niemożliwe będzie wygenerowanie populacji równej 100 dla 3 kandydatów gdyż 3! jest równe 6). W przypadku wystąpienia tej sytu-

acji przedział populacji jest odpowiednio zmniejszany do połowy wielkości permutacji wygenerowanych kandydatów. Twórcy algorytmu staneli też przed wyborem jaki zakres permutacji powinien być generowany: czy zawsze ze stałą długością, czy też ze zmienną długością. Ostatecznie zdecydowaną się na drugą opcję gdyż stwierdzono, że sama długość powstałych wcześniej kandydatów też może być różna, a dodanie osobnej funkcji określającej długość stwor-

zonej sekwencji nie przyniosłaby żadnego pozytywnego wpływu jeśli chodzi o jakość i czas przetwarzania algorytmu. Przykładowe wygenerowanie populacji o wielkości 20 dla sekwencji *n*=30 znajduje się poniżej.

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.013.png)

Przykład wygenerowanej losowo populacji o długości 20 dla sekwencji n=30

Informacje o tworzonej sekwencji są przechowywane w zmiennej typu Dataframe, który swoją budową przypom- ina trochę strukturę bazodanową (w danym typie danych możemy przechowywać różne informacje w kolumnach, które posiadają swoje indetyfikatory. Dodatkowo wyszukiwanie informacji w danym typie jest dzięki indeksowaniu o wiele wydajniejsze). Wykorzystany Dataframe w algorytmie składa się z dwóch kolumn: kolumny [’SEQUENCE] - przechowywującej wartość typu *string* zbudowanej sekwencji oraz [’PATH’] przechowywującej informację z jakich kandydatów jest zbudowana. Warto zauważyć, że przy tworzeniu permutacji nie posiadamy żadnych powtórzeń, gdyż każdą sekwencję możemy użyć tylko raz w budowaniu danego łańcucha. Sam proces tworzenia permutacji opiera się na permutacji indetyfikatorów ścieżek wygenerowanych wcześniej kandydatów (kolumna ”*PATH* ”). Następnie po stworzeniu ścieżki składającej się z ”numerów” możemy zbudować sekwencję składającą się z liter. Zbudowana w ten sposób sekwencja przyda się nam w nastęnym etapie w którym będziemy mogli dokonać oceny zbudowanych sekwencji oraz selekcji wygenerowanych osobników.

3. **Sposób oceny osobników populacji oraz sposób selekcji**

Posiadając ścieżkę stworzonego osobnika możemy przejść do stworzenia go jako obiektu. Nowy osobnik w naszym algorytmie będzie zawsze nowym obiektem typu *Text*. Warto zauważyć, że definicja tej klasy dziedziczy z klasy abstrakcyjnej **Element** opis konstruktora, który mówi, że każdy stworzony obiekt powinien posiadać dopisaną do siebie wartość oceny *self.fitness*. Sama wartość oceny dla danego obiektu jest tworzona za pomocą metody **evalu- ate function** znajdującą się w klasie **Text**. Sam system oceny jest dosyć prosty i opiera się na dodawaniu do ogólnej oceny obiektu wartości wykrywanych pokryć dla danych łączeń kandydatów. Im większa wartość oceny tym lepsze czyli większe znalezione pokrycia między dwoma kandydatami. Warto wspomnieć też, że sekwencje które posiadają maksymalnie możliwe pokrycia są promowane podwójnie (czyli do funkcji oceny jeśli znajdziemy dla sekwencji o dłu- gości 15 pokrycie równe 14 to do końcowego wyniku dodawany jest wynik 14\*2 = 28). Taki sposób promocji ma na celu wyznaczanie kandydatów, którzy posiadają pokrycia maksymalne, gdyż bardzo łatwo mogłoby dojść do sytuacji gdzie znajdujemy sekwencję np. podsiadającą *n* pokryć równą *k*-2 (gdzie *k* oznacza długość sekwencji) oraz znaleźć drugą sekwencję, która posiada tylko jedno pokrycie *k*-1. W przytoczonym przykładzie to pierwsza sekwencja byłaby bardziej promowana zamiast drugiej która znalazła pokrycie maksymalne. Podwójna promocja pokryć maksymalnych

ma za zadanie zniwelowanie powstałych różnić podczas oceniania powstałego obiektu. Posiadając już w pełni wygen- erowaną populację (czyli o takiej wielkości jaka została zadeklarowana przez klasę **GenethicAlghorithm**) program przystępuje do selekcji powstałych osobników. Wybierany jest *x*% najlepszych osobników w oparciu o przypisaną im wartość *fitness*. Wartość *x* w programie opisana jest przez zmienną *selection percent* i jest możliwa do zmiany przez użytkownika z poziomu pliku main. Warto zauważyć, że mimo posiadania łańcucha sekwencji typu string wygen- erowanej w czasie losowania permutacji, metoda **evaluate function** ponownie dokonuje procesu budowania całej sekwencji od nowa. W celu przyśpieszenia danego programu jest możliwość eliminacji podwójnego etapu budowania ścieżek jednak dla jasności oraz czytelności kodu w oddawanej pracy została ona w takiej postaci.

4. **Sposób krzyżowania się osobników z populacji**

Po wybraniu najlepszych osobników z danej populacji algorytm przechodzi do budowania następnej populacji. W tym celu z klasy **GenethicAlghotihm** zostaje wywołana metoda **crossover**, która odpowiada za krzyżowanie danych populacji. Dla każdej z sond algorytm pobiera dwóch osobników (może pobrać nawet dwóch takich samych) oraz atrybut *self.cross number* w celu określenia w ilu miejscach ma nastąpić podział danej ścieżki. Podobnie jak to było w przypadku generowania permutacji program sprawdza i odpowiednio zmniejsza liczbę *self.cross number* w przypadku gdy liczba przecięć jest większa niż liczba możliwych podziałów w stworzonych ścieżkach. Sam podział ścieżek opiera

się na losowym podziale ścieżek *PATH* (posiadającą numerację odpowiadającym im kandyatom) dwóch osobników

w celu następnie losowego kolejno ich łączenia. Warto podkreślić fakt, że dane łączenie bierze pod uwagę kolejność łączeń. Dla dwóch sekwencji podzielonych w 3 miejscach wygenerujemy 4 części sekwencji. Każda z ich część będzie łączona losowo z każdego punktu podziału. Mając np. części 1 2 3 4 w pierwszej sekwencji oraz 5 6 7 8 w drugiej, w celu tworzenia sekwencji będą brane pod uwagę kolejno sekwencje: [1,5], [2,6], [3,7], [4,8]. Poniżej znajduje się przykład ilustrujący proces łączenia się dwóch ścieżek.

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.014.png)

Przykładowe łączenie dwóch ścieżek

Sekwencja [0, 8, 6, 3, 7, 5] oraz [3, 0, 2, 4] zostały podzielone w losowych miejscach. Sekwencja pierwsza została podzielona na części: [0], [8, 6, 3], [7], [5], a sekwencja druga na [3], [0], [2], [4]. Następnie losowo ze sobą zostały połączone wybrane fragmenty. Z sekwencji drugiej została pobrana pierwsza część, a z drugiej część druga, trzecia oraz czwarta. Możemy zauważyć jednak niestety, że w wyniku tych operacji może wystąpić zjawisko wystąpienia duplikacji (co jak było wspomniane wcześniej nie jest sytuacją możliwą w przypadku gdy chcemy znaleźć prawidłowe rozwiązanie). W celu eliminacji tego zjawiska, algorytm dokonuje detekcji danych powtarzających się punktów oraz generuje listę kandydatów którzy w danej ścieżce nie byli wykorzystani. Następnie algorytm dokonuje losowego wyboru elementu z listy, usuwa wykorzystany element z listy kandydatów nie wykorzystanych i podstawia pod powtarzający

się element. Proces ten powtarza się dopóki wszyscy zduplikowani kandydaci nie zostaną wyeliminowani. Ostatecznie metoda crossover może utworzyć nowy obiekt będący obiektem typu **Text**. będący potomkiem starego pokolenia, a następnie przystąpić do oceny powstałego obiektu.

5. **Omówienie błędów implementacyjnych oraz sugestia polepszenia algorytmu dla metody crossover**

Analizując zasadę działania **crossover** możemy zauważyć, że część pracy jest wykonywana dodatkowo i mogłaby zostać wyeliminowana. Chodzi dokładniej o typ danych które zostaje dzielony czyli ścieżkę typu PATH (posiadającą numerację kandydatów). Problem ten zaczyna być zauważalny dla bardzo dużych instacji i dużej ilości kandydatów. Jest to też jedna z przyczyn dlaczego algorytm był testowany dla stosunkowo małej ilości populacji (wynoszącej średnio do 100 osobników na populację). Po wygenerowaniu ścieżki typu PATH i wyeliminowaniu duplikatów, cała sekwencja budowana jest od początku co zajmuje niepotrzebny czas. Czyli dla 3 krzyżowań w przypadku dla 500 kandydatów i tak będziemy musieli znaleźć 499 wartości pokryć oraz wykonać 499 operacji łączeń kandydatów. Dla coraz większej ilości kandydatów problem ten staje się coraz bardziej widoczny. Proponowaną alternatywą jest dodanie dwóch dodatkowych kolumn do wartości dataframe. W oddawanej wersji typ ten dla generowanych naszych osobników posiada dwa pola: [’SEQUENCE’] oraz [’PATH’]. Dodanie nowych wartości [’NEW~~ STRING’] oraz [’COVER’] umożliwiłoby przyśpiesze- nie pracy danej metody. Pole [’NEW~~ STRING’] zbierałoby informacje oznaczającą jaka była wielkość łańcucha po każdym dołączeniu nowej sekwencji, a pole [’COVER’] jakie było między nimi pokrycie. Ostatecznie posiadając in- formację o: miejscu przecięcia, nowym dołączonym kandydacie, wielkości łańcucha gdzie znajduje się kandydat oraz pokryciu bylibyśmy w stanie poprawnie dokonać przecięcia oraz dołączenia nowej części sekwencji z innej populacji. Oznacza to, że liczba operacji potrzebnej do zbudowania nowego łańcucha (nie uwzględniając podstawieniach nowych wartości do [NEW~~ STRING] oraz [COVER]) zakładając, że nie wystąpi żadna duplikacja kandydatów zmalałaby do

3 operacji. Część pomysłów w celu przyśpieszenia danego algorytm jest wciąż w oddanym kodzie, jednak ze względu

na ograniczenia czasowe została ona zakomentowana i nie pełni ona teraz żadnej funkcji.

6. **Działanie mutacji**

Podczas tworzenia obiektu pobierając informację o częstotliwości występowania mutacji z klasy GenethicAlghorithm self.mutation~~ probability, może dojść do wystąpienia mutacji na świeżo stworzonym obiekcie. Metoda perform~~ mutation bierze![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.015.png) stworzoną ścieżkę nowego osobnika i na jej postawie tworzy ścieżkę kandydatów którzy jeszcze nie zostali wyko- rzystani. Następnie wybiera jeden element z danej listy i zastępuje jeden element ze stworzonej ścieżki. Twórcy pro-

gramu mogli wybrać dwie ścieżki działania mutacji: dokonania zamiany losowego elementu z wygenerowanej ścieżki i dokonanie oceny wykonanej mutacji lub też zamieniania po kolei jednego elementu z wygenerowanej ścieżki i wybrania mutacji która osiągneła najlepszy wynik dopasowania. Autorzy zdecydowali się na drugą wersję programu.

7. **Omówienie błędów implementacyjnych oraz sugestia polepszenia algorytmu dla metody muta- tion**

Sposób zaimplementowanej mutacji mimo bardzo dużego nakładu czasowego daje gwarancję, że dla losowo wybranej ścieżki (jeśli isniteje dla niego dopasowanie) to dla niego takie dopasowanie znajdzie. Zdecydowanie się na mutację każdego elementu po kolei niestety niesie ze sobą także wzrost czasu obliczeniowego potrzebnego na wykonanie funkcji

perform~~ mutation. Przykład rozwiązania tego problemu został juz przedstawiony w sekcji **4.3.5 Omówienie błędów implementacyjnych![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.016.png) oraz sugestia polepszenia algorytmu dla metody crossover**.

8. **Ulepszanie kandydatów po wygenerowaniu nowej populacji**

Każdy osobnik oprócz atrybutu fitness posiada także atrybut self.found~~ max~~ cover. Umieszczane tam są znalezione maksymalne dopasowania jakie udało się znaleźć podczas generowania nowej populacji. Jeśli w czasie selekcji, wyz- naczeni osobnicy posiadają znalezione maksymalne dopasowania to istniejący kandydaci zostają ”ulepszeni” to znaczy dwóch kandydatów zostaje ze sobą złączonych w jednego kandydata.

9. **Warunek stopu**

Wraz z ulepszaniem się kandydatów kiedy zostanie wygenerowany kandydat który swoją długością będzie równy szukanej sekwencji oznacza to, że szukana sekwencja została znaleziona. Gdy odpowiednio dla obu sond zostaną znalezione sekwencje, program przystępuje do dekodowania obu sekwencji do sekwencji końcowej i następuje koniec programu.

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.017.png)

Przykład wyjścia programu i zakodowanego wyniku

4. **Omówienie dodatkowych zależności występujących w kodzie**

W przesłanej pracy ze względu na charakter działania algorytmu genetycznego, który do poprawnej pracy nie potrze- buje komunikacji i wymiany wyników między dwoma sondami zdecydowano się podzielić pracę dla sondy SW oraz sondy RY na dwa potomne procesy. Dodatkowo, w celu polepszenia czasu działania algorytmu genetycznego w przy- padku gdy przez określoną ilość populacji algorytm nie będzie w stanie znaleźć ani maksymalnego dopasowania, ani polepszyć swojego wyniku, następuje wygenerowanie całkowicie nowej populacji. Testy wykazały, że implementacja tego mechanizmu znacznie przyśpiesza czas dojścia do wyniku.

5. **Omówienie występujących błędów krytycznych w algorytmie**

Analizując pseudokod zaimplementowanego algorytmu genetycznego oraz porównując go do np. zaprezentowanej przykładowej implementacji algorytmu na zajęciach możemy zauważyć, że została dodana do kodu nowa funkcjonal- ność. Mianowicie w przypadku znalezienia maksymalnego pokrycia następuje scalanie ze sobą kandydatów. Podejście to zostało zaproponowane przez autorów aby polepszyć szybkość i czas działania programu, gdzyż bez istnienia tej funkcji, algorytm nie był w stanie dojść do prawidłowego rozwiązania. Mogła odpowiadać za to źle zaimplementowana funkcja oceny, jednak autorom nie udało się ustalić gdzie leży przyczyna aż tak długiego czasu działania programu, i po wielokrotnym przerabianiu funkcji oceny zdecydowano się na zaimplementowanie operacji scalania kandydatów gdzie znaleziono maksymalne pokrycie. Rozwiązanie to działało, dla wygenerowanych wcześniej przykładach, jednak wraz z próbą zmniejszania długości wszystkich sekwencji algorytm zaczął czasami zwracać błędne wyniki. Autorzy doszli do wniosku, że z powodu zaimplementowania scalania kandydatów algorytm nie będzie w pełni działał poprawnie dla wygenerowanych sond w których znajduje się para sekwencji składająca się z identycznych liter i różniąca się tylko ostatnią literą. Jest to spowodowane faktem, że dla obu takich sekwencji z pary znalezione dopasowanie zawsze będzie poprawne, co dla algorytmu oznacza, że nie ma znaczenia którego (losowego) kandydata weźmie do łączenia dwóch kandydatów w celu stworzenia lepszej skróconej listy kandydatów. Autorzy kodu nie przewidzieli zaitnienia tej sytuacji i przepraszają za pojawienie się tego błędu.

W celu wyeliminowania danego negatywnego zjawiska można zastosować kilka metod. Pierwszą, a zarazem najprostszą wydaje się wyeliminowanie systemu polepszania kandydatów. Jednakże nawet jeśli dany problem wyeliminujemy to w takiej sytuacji dany problem wciąż wystęujpe i łatwo można wyobrazić sobie sytuację gdy w procesie budowania i promowania coraz to lepiej dopasowanych sekwencji zostanie dołączona sekwencja która nie powinna się tam znaleźć ale wraz z np. innymi błędami pozytywnymi tworzy bardzo długą sekwencję z maksymalnymi pokryciami. Wygen- erowana w ten sposób sekwencja będzie bardzo wysoko promowana i istnieje możliwość, że generowana sekwencja też z maksymalnym pokryciem ale z poprawną literą na końcu będzie odrzucana. Ostatecznie, autorzy stwierdzają, że mogłaby to być jedna z możliwości wyeliminowania danego błędu jednakże nie są oni do końca pewni czy ta metoda byłaby do końca poprawna.

Inną metodą proponowaną przez autorów byłoby sprawdzanie już na etapie znalezienia maksymalnego pokrycia przez jakiegoś kandydata, sprawdzenia czy w danej sondzie nie znajduje się inny kandydat, składający się z takich samych liter oprócz ostatniej. W przypadku znalezienia takich kandydatów, algorytm generowałby x ilość dołączeń do sek- wencji.

Na przykład:

Posiadając sekwencję **SWWWWSSSWSWSS** i znajdując jej dopasowanie **WWWWSSSWSWSST** sprawdzamy czy nie istnieją inni kandydaci zakończeni inną literą. W przypadku gdy znajdziemy kandydatów np. **WWWWSSS- WSWSSA** oraz **WWWWSSSWSWSSC** posiadamy 3 potencjalnych kandydatów którzy mają pokrycie maksy- malne. Na etapie więc polepszenia kandydatów wykorzystujemy pozostałych kandydatów (których algorytm w czasie swojej realizacji bezpośrednio nie znalazł ale są podobni do znalezionego kandydata) i generuje 3 ulepszonych kandy- datów: **SWWWWSSSWSWSST**, **SWWWWSSSWSWSSA**, **SWWWWSSSWSWSSC**.

Dalsza praca algorytmu powinna ułożyć się tak, że oczywiscie tylko ta sekwencja która jest poprawna doprowadzi program do prawidłowego wyniku i zakończy pomyślnie program. Pomysł ten był bliski implementacji przez autorów, jednakże znów przez ograniczenia czasowe nie został on zrealizowany.

6. **Szacowana złożoność czasowa**

Oszacowanie złożoności czasowej zaimplementowanego algorytmu autorzy szacują, że wynosi ona: **O(G SW(PM) + G RY(PM))**:

- G~~ SW oraz G~~ oznaczają odpowiednio liczbę generacji SW oraz RY
- P oznacza liczbę populacji
- M oznacza liczbę mutacji
7. **Testy**
1. **Metodologia wykonywanych testów**

W testach wzieła udział wygenerowana wcześniej grupa sekwencji wspomniana w **3.1 Zbiór danych wstępnych**. W celu wyeliminowania rozbieżności wyników dla każdej sekwencji zostało wykonanych 30 testów i jako wynik została padana uzyskana wartość arytmetyczna. Algorytm został przetestowany dla róznych parametrów początkowej popu- lacji oraz prawdopodobieństwa wystąpienia mutacji. W testach wartości te wyniosły odpowiednio:

- populacja = [10, 20, 40]
- mutacja = [0.00, 0.05, 0.1, 0.15]
2. **Wyniki**

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.018.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.019.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.020.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.021.png)

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.022.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.023.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.024.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.025.png)

![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.026.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.027.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.028.png) ![](Aspose.Words.96726926-4c80-4359-9e55-f17fbf1ad50d.029.png)

Analizując uzyskane wyniki możemy zauważyć, że wspomniane przez autorów problemy wydajnościowe jeśli chodzi o operacje **crossover** oraz zwłaszcza **mutation**. Wzrost występowania mutacji znacznie zwiększa ilość operacji która w rezultacie nie przyśpiesza dojścia do rozwiązania. Zmienienie zasady działania mutacji (na mutowanie jednego zamiast wszystkich i brania najlepszej mutacji) mogłoby przyśpieszyć operacje związane z mutacją. Jednak wciąż odsetek występowania mutacji powinien być nie wielki i dalsze testy pokazały, że mutacja na poziomie nie przekraczającej 3% jest dosyć jeszcze dobrą wartością dobraną empirycznie, która nie ma tak dużego wpływu na wynik. Testy wykazały także, że wzrost ilości populacji nie ma też jakiegoś dużego znaczenia jeśli chodzi o czas dojścia do danego wyniku. Jednak trzeba mieć na względzie fakt, że populacje które wzieły udział w teście były niewielkie i że dla bardzo dużych sekwencji z błędami pozytywnymi wzrost populacji mógłby mieć lepszy wpływ na polepszenie wyników obliczeń.

5  **Oznaczenie zmiennych w algorytmie oraz sposób uruchamiania pro- gramu**

Aby uruchomić odpowiednio algorytm dokładny należy uruchomić plik *main iter.py* natomiast aby uruchomić al- gorytm genetyczny należy uruchomić *setup.py*. W algorytmie dokładnym możemy zmienić plik wejściowy poprzez podstawienie odpowiedniego argumentu fla funkcji *fun iter*. W algorytmie genetycznym do zmiany w pliku *main.pyx* mamy takie zmienne jak:

- **population number** - liczba osobników na populację
- **mutation probability** - prawdopodobieństwo mutacji
- **size sequence** - numer wczytanego pliku z folderu **sample cases/positive numbers max k max sqt max pos/**
- **selection percent** - procent wybrania osobników do nowej populacji
- **TIMEOUT** - czas po jakim program jest przerywany
- **number of average** - liczba testów wykonywana w testach
- **run tests** - wybór czy chcemy uruchomić testy czy też raz wykonać program

Dodatkowo w pliku log.pyx znajdują się zmienne za włączenie debugowania działania programu.

- **DEBUG** - wyświetlanie podstawowych informacji na temat analizowanej sekwencji
- **DEBUG GENERATION** - wyświetlanie informacji o nowych generacjach
6  **Podsumowanie wyników oraz omówienie pracy działania algorytmów**

Analizując najlepsze osiągnięte wyniki pomiędzy dwoma algorytmami możemy zauważyć, że algorytm genetyczny osiąga podobne czasy (jednak zawsze wolniejsze) czasy w porównaniu do algorytmu dokładnego. Wyższy czas wynika z zastosowanych większych struktur danych oraz większej potrzebnej do wykonania wstępnych operacji które rzutują na dłuższy czas wykonania. Widać jednak, że dla coraz co większych sekwencji czas wykonywania algorytmu dokładnego rośnie o wiele szybciej w porównaniu do algorytmu genetycznego. Oznacza to, że dla większych sekwencji powinniśmy zastosować algorytm genetyczny.
13
