import random
import sys

input_file = sys.argv[1]
cities_no = int(sys.argv[2])

random.seed(1)

roads = []
for i in range(0, cities_no):
    roads.append([])
    for j in range(0, cities_no):
        roads[i].append(0)

for i in range(0, cities_no):
    for j in range(i + 1, cities_no):
        roads[i][j] = random.randint(1, 30)

for i in range(0, cities_no):
    for j in range(0, i):
        roads[i][j] = roads[j][i]

f = open(input_file, "w")
f.write(str(cities_no) + "\n")
for i in range(0, cities_no):
    for j in range(0, cities_no):
        f.write(str(roads[i][j]) + " ")
    f.write("\n")
