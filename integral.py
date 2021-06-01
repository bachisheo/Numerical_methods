import matplotlib.pyplot as plt
x = []
y = []
with open("output.txt") as fin:
    while True:
        line = fin.readline()
        if not line:
            break
        x.append(list(map(float, line.strip().split())))
        line = fin.readline()
        y.append(list(map(float, line.strip().split())))

plt.plot(x[3], y[3],'.', label="in")
plt.plot(x[4], y[4],'.', label="out")
plt.plot(x[0], y[0],'.', label="Исходные точки")
plt.plot(x[1], y[1], "-r", label="Ограничивающая область")
plt.plot(x[2], y[2], label="Интерполирование")


plt.legend()

plt.show()


