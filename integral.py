import matplotlib.pyplot as plt
x = []
y = []
labels = []
labels.append('Область определения')
labels.append('Фигура')
with open("output.txt") as fin:
    while True:
        line = fin.readline()
        if not line:
            break
        x.append(list(map(float, line.strip().split())))
        line = fin.readline()
        y.append(list(map(float, line.strip().split())))

plt.plot(x[0], y[0], "-r", label="Ограничивающая область")
plt.plot(x[1], y[1]);
plt.legend()

plt.show()


