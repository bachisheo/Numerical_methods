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

for curPlot in range (1, x.__len__(), 2):
    # plt.plot(x[curPlot], y[curPlot], label=labels[curPlot])
    plt.plot(x[curPlot], y[curPlot], ".r")
    plt.plot(x[curPlot + 1], y[curPlot + 1])
plt.legend()

plt.show()


