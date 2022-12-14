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
	In cadrul problemei studiate, variabila care sta la baza calculului complexitatii, si in functie de care am realizat si testele este numarul de orase prin care negustorul doreste sa treaca.

## Implementarea secventiala

### Implementarea naiva
	- Prima metoda de implementare secventiala este bazata pe backtracking (aceasta ar fi varianta naiva de implementare). Cu toate acestea, aceasta nu se dovedeste a fi o varianta buna de implementare, deoarece, pentru o dimensiune a problemei care depaseste 100 de orase, timpul de rulare creste considerabil

	- Complexitate: O(N!)

### Implementare bazata pe algoritm genetic

	- Un algoritm genetic este un algoritm prin intermediul caruia se doreste obtinerea unei solutii cat mai bune pentru o problema, pornind de la o solutie mai putin buna a acesteia (de obicei pornind de la o solutie aleasa random). Pe masura ce algoritmul genetic ruleaza, solutia curenta este imbunatatita.

	- Urmatoarele notiuni sunt importante pentru un algoritm genetic:
		-> Cromozom - o mica bucata din solutie care, in problema de fata, reprezinta un oras
		-> Individ - reprezinta o solutie pentru problema noastra, este compus din mai multi cromozomi. In cadrul problemei noastre, mai multe orase (cromozomi) alaturate formeaza o ruta (individ)
		-> Populatie - o serie de indivizi, intotdeauna se va lucra cu o lista de indivizi, lista din care cei mai buni vor fi pastrati, iar ceilalti vor fi inlocuiti cu variatii ale celor mai buni

	- Varianta de algoritm genetic din aceasta tema realizeaza urmatorii pasi pentru a ajunge la o solutie buna:
		-> genereaza o prima generatie de indivizi (random). Doar primul individ din aceasta generatie este realizat printr-o parcurgere a fiecarui oras si alegerea drumului minim din fiecare oras.
		-> 

### Implementare bazata pe programare dinamica

## Implementare paralela

	- Implementarea paralela este bazata pe paralelizarea implementarii algoritmului genetic (implementare ce a fost prezentata anterior).

### OpenMP
### Pthreads
### MPI
### Hibrid (MPI + OpenMP)

## Concluzii
## Detalii rulare, organizare arhiva
## Feedback tema


