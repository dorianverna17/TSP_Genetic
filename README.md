# TSP_Genetic - Tema SM Anul 4 - CTI

Dranca Stefana-Ioana 341C1 <br/>
Verna Dorian-Alexandru 341C1 <br/>

# Travelling Salesman Problem

## Scurta introducere

Travelling Salesman Problem sau Problema Comis Voiajorului este o problema care se
enunta astfel:

Un vanzator ambulant vrea sa parcurga o serie de orase pentru a-si vinde bunurile.
El vrea sa gaseasca ruta cea mai scurta astfel incat sa poata vizita toate orasele
si sa se poate intoarce in orasul de unde a plecat. Gasiti insiruirea de orase in
ordinea in care acestea trebuie strabatute astfel incat costul sa fie minim

Exp.

Avem 4 orase, iar costurile pentru a strabate distanta dintre oricare doua dintre
orasele acestea sunt afisate prin matricea urmatoare:

0 10 15 20 <br/>
10 0 35 25 <br/>
15 35 0 30 <br/>
20 25 30 0 <br/>

Rezultatul pe care al trebui sa il obtinem ar fi urmatorul:<br/>
cost: 80<br/>
ruta: 1 -> 2 -> 4 -> 3 -> 1

## Dimensiunea problemei
- In cadrul problemei studiate, variabila care sta la baza calculului complexitatii, si in functie de care am realizat si testele este numarul de orase prin care negustorul doreste sa treaca.

## Implementarea secventiala

### Implementarea naiva
- Prima metoda de implementare secventiala este bazata pe backtracking (aceasta ar fi varianta naiva de implementare). Cu toate acestea, aceasta nu se dovedeste a fi o varianta buna de implementare, deoarece, pentru o dimensiune a problemei care depaseste 100 de orase, timpul de rulare creste considerabil

- Complexitate: O(N!)

### Implementare bazata pe algoritm genetic

- Un algoritm genetic este un algoritm prin intermediul caruia se doreste obtinerea unei solutii cat mai bune pentru o problema, pornind de la o solutie mai putin buna a acesteia (de obicei pornind de la o solutie aleasa random). Pe masura ce algoritmul genetic ruleaza, solutia curenta este imbunatatita.

- Urmatoarele notiuni sunt importante pentru un algoritm genetic:
	- Cromozom - o mica bucata din solutie care, in problema de fata, reprezinta un oras
	- Individ - reprezinta o solutie pentru problema noastra, este compus din mai multi cromozomi. In cadrul problemei noastre, mai multe orase (cromozomi) alaturate formeaza o ruta (individ)
	- Populatie - o serie de indivizi, intotdeauna se va lucra cu o lista de indivizi, lista din care cei mai buni vor fi pastrati, iar ceilalti vor fi inlocuiti cu variatii ale celor mai buni

- Varianta de algoritm genetic din aceasta tema realizeaza urmatorii pasi pentru a ajunge la o solutie buna:
	- genereaza o prima generatie de indivizi (random). Doar primul individ din aceasta generatie este realizat printr-o parcurgere a fiecarui oras si alegerea drumului minim din fiecare oras.
	- Intra in loop-ul care trece prin toate generatiile. Aici executa urmatoarele operatii:
		- calculeaza fitness-ul fiecarui individ
		- sorteaza indivizii crescator (fitness-ul cel mai scazut este cel mai bun) (qsort)
		- aplica algortmul de mutatie pe populatia de indivizi
		- schimba generatiile astfel incat populatia curenta sa fie cea obtinuta prin mutatie
	- Odata iesit din loop, primul individ din generatia finala trebuie sa fie destul de apropiat sau sa fie chiar solutia problemei.

- In cazul unui algoritm genetic de acest tip, se poate ajunge la solutie dupa un numar destul de mare de epoci alese
- In cazul nostru, am folosit 1000 de epoci (generatii)
- Detalii despre algoritmul de mutatie:
	- Primii 30% din indivizii din generatie (cei mai buni 30%) sunt pastrati in generatia urmatoare
	- Urmatorii 30% sunt obtinuti prin interschimbarea a doua gene a celor 30% indivizi alesi precedent
	- Urmatorii 30% sunt obtinuti prin interschimbarea a patru gene a celor 30% alesi prima oara
	- Restul indivizilor sunt obtinuti prin interschimbarea celor patru gene, similar cu pasul anterior, dar in alta ordine (tot a celor 30% indivizi alesi prima oara)
	- Pozitiile genelor unde are loc interschimbarea sunt determinate random si verificate mai apoi (pentru a ne asigura ca nu sunt egale)


### Implementare bazata pe programare dinamica

- Aceasta implementare nu este tratata aici, dar detalii despre aceasta pot fi gasite aici: https://www.geeksforgeeks.org/travelling-salesman-problem-using-dynamic-programming/

## Implementare paralela

- Implementarea paralela este bazata pe paralelizarea implementarii algoritmului genetic (implementare ce a fost prezentata anterior).

- Functiile din bucla de generatii care sunt realizate in paralel sunt urmatoarele:
	- compute_generation_fitness_openmp - fiecare thread va calcula fitness-ul pentru indivizii care se afla in start si end
	- mutate_generation_openmp - fiecare thread va face mutatia pe fiecare individ care se afla intre start si end
- Thread-ul 0 se va ocupa de sortarea generatiilor si de generarea pozitiilor random care sunt folosite la mutatie
- Tot ce este generat random, este generat de catre thread-ul 0, deaorece se doreste obtinerea acelorasi rezultate la rularile tuturor variantelor de implementare 

### OpenMP

- In cadrul implementarii OpenMP au fost paralelizate urmatoarele sectiuni de cod:
	- Sectiunea paralelizata este reprezentata de bucla care itereaza prin toate generatiile. Toata aceasta bucata de cod este inclusa intr-un bloc marcat cu directiva "#pragma omp parallel". In cadrul acestei sectiuni, este calculat id-ul thread-ului, indexul de start si indexul de sfarsit pentru indivizii din populatie de care acest thread se va ocupa in procesare.

### Pthreads

- In cadrul implementarii Pthreads au fost paralelizate urmatoarele sectiuni de cod:
	- Sectiunea paralela este aceeasi ca cea de la OpenMP, din acest punct de vedere, cele doua implementari nu prezinta nicio diferenta, ambele paralelizeaza aceleasi secvente de cod + functii.

### MPI
### Hibrid (MPI + OpenMP)

## Concluzii

- Pe parcursul implementarii temei, am observat ca o parte din sectiunile de cod prezente nu se merita a fi paralelizate, de exemplu, am realizat o implementare a sortarii mergesort paralelizata, dar timpii erau mai mari atunci cand se folosea. Din rezultatele obtinute la profiling am aflat ca se petrecea considerabil mai mult timp in functia de comparare a indivizilor, decat atunci cand se folosea functia qsort built in.
- Am mai incercat un approach in care pastram primii 30% indivizi, apoi facem mutatie pe urmatorii 60%, apoi pe restul in regeneram in totalitate random. Tot pe profiling am descoperita ca se pierdea prea mult timp la apelurile functiei rand(), asa ca am decis sa fac o mutatie similara si pe ultimii indivizi astfel incat sa obtinem rezultate cat mai bune.
- Profiling am facut cu IntelVTune si cu valgrind, mai multe detalii sunt in directorul profiling

- Timpii sunt mai mari pe implementarile paralele pentru input-uri mici. Aparent, avem un overhead destul de mare legat de pornirea thread-urilor si de sincronizarea acestora. Se obtine speedup doar atunci cand dimensiunea problemei se mareste. Se poate observa in partea de profiling cum timpul de rulare scade in favoarea implementarilor paralele.

## Detalii rulare, organizare arhiva

- Proiectul se poate rula in urmatoarele moduri:
	- Utilizand comanda "./main <varianta_de_rezolvare> <fisier_input>", unde:
		- Varianta de rezolvare este una dintre urmatoarele:
			- sequential_naive
			- sequential_genetic
			- parallel_openmp
			- parallel_pthreads
			- parallel_mpi
			- parallel_mpi_omp
		- Fisierul de input reprezinta calea catre fisierul de input
	- Utilizand scriptul run.sh ("./run.sh"). Acesta executa toate comanda anterioara pentru toate input-urile date si folosind toate variantele de rezolvare
	- Utilizand regulile de rulare din fisierul Makefile

- In cadrul arhivei avem urmatoarea structura:
	- Directorul input - contine testele pe care le-am generat pentru a testa functionalitatea algoritmului si a paralelizarilor
	- generate_input.py - fisierul cu ajutorul carora se genereaza testele (python3 generate_input.py <file_to_generate> <no_cities>)
	- main.c - fisierul sursa unde se afla functia main, entry_point-ul programului nostru
	- Makefile - fisierul are o regula de build, una de clean si mai multe reguli de rulare in cazul in care se doreste a fi folosite
	- run.sh - fisier ce ruleaza toate implementarile pe toate input-urile existente
	- README.md
	- Directorul utils - contine doua fisiere helper:
		- genetic_utils.h - contine mai multe functii care sunt folosite de toate variantele de implementare (secvential + paralel)
		- graph.h - contine definirea structurii care retine configuratia oraselor (o retinem sub forma de matrice), precum si cateva functii pentru citire si printare
	- Directorul sequential - contine urmatoarele implementari secventiale:
		- implementarea naiva - in directorul naive avem TSP.h
		- implementarea cu algoritm genetic - in directorul genetic avem TSP.h
	- Directorul parallel - contine urmatoarele implementari paralele:
		- openmp - in directorul openmp avem TSP.h
		- pthreads - in directorul pthreads avem TSP.h
		- mpi - in directorul mpi avem TSP.h
		- hibrid (mpi + openmp) - in directorul mpi_omp avem TSP.h
	- Directorul profiling - contine detalii, screenshot-uri legate de partea de profiling

## Feedback tema

Tema a fost ok, in general dificultatea a depins de problema aleasa, iar problema aleasa a fost potrivita pentru a lucra doua persoane la ea. A fost un proiect bun si de folos pentru a consolida cunostiintele legate de paralelizarea anumitor implementari si evaluarea performantelor acestora.
