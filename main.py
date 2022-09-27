import scipy.stats as sps
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
size = 200
v = 94
p = 0.25 + 0.005*v
lamb = 0.8 + 0.02*v

def binomial_law(p):
    sample = sps.binom(n=5+v%16, p=p).rvs(size=size)
    print('Выборка: \n', sample[:size])
    sample.sort()
    print('Выборка сортированная: \n', sample[:size])

    xi = []
    ni = []
    wi = []
    s = []

    xi.append(0)
    ni.append(1)
    wi.append(0)
    s.append(0)

    for i in range(size):
        if not i == size - 1:
            if not sample[i] == sample[i + 1]:
                xi[len(xi) - 1] = sample[i]
                ni[len(ni) - 1] = list(sample).count(sample[i])
                wi[len(wi) - 1] = ni[len(ni) - 1] / size
                s[len(s) - 1] = round(sum(wi), 5)
                ni.append(1)
                wi.append(0)
                s.append(0)
                xi.append(0)
        else:
                xi[len(xi) - 1] = sample[i]
                ni[len(ni) - 1] = list(sample).count(sample[i])
                wi[len(wi) - 1] = ni[len(ni) - 1] / size
                s[len(s) - 1] = round(sum(wi), 5)

    vib_sr = 0
    vib_dis = 0
    vib_kvad = 0
    print("-----------------------------")
    print("Эксперементальные значения: ")
        #Выборочное среднее
    for i in range(len(xi)):
        vib_sr += xi[i] * wi[i]
    print("Выборочное среднее: ", round(vib_sr,5))
        #Выборочная дисперсия
    for i in range(len(xi)):
        vib_dis += pow((xi[i] - vib_sr), 2)*wi[i]
    print("Выборочная дисперсия: ", round(vib_dis,5))
        #Выборочное срднее квадратическое откланение
    vib_kvad = pow(vib_dis, 0.5)
    print("Выборочное срднее квадратическое откланение: ", round(vib_kvad,5))

    ni_count = 0
    ni_buf = []
    for i in range(len(ni)):
        ni_buf.append(False)

    vib_mod = 0
    index = 0
        #Выборочная мода
    for i in range(len(ni_buf)):
        if max(ni) == ni[i]:
            ni_buf[i] = True
            ni_count += 1
            index = i
        else:
            ni_buf[i] = False

    if ni_count == 1:
        vib_mod = xi[index]
    else:
        for i in range(index - ni_count + 1, index):
            if ni_buf[i] == ni_buf[i+1]:
                vib_mod = (xi[index - ni_count + 1] + xi[index])/2
            else:
                vib_mod = -1
                break
    if not vib_mod == -1:
        print("Выборка по моде: ", round(vib_mod, 5))
    else:
        print("Выборка по моде: Не существует")

    vib_med = 0

    for i in range(len(s)):
        if s[i] > 0.5:
            vib_med = xi[i]
            break
        if s[i] == 0.5:
            vib_med = (xi[i]+xi[i+1])/2
            break
    print("Выборочная медиана: ", round(vib_med, 5))

    u3 = 0
    for i in range(len(xi)):
        u3 += pow(xi[i], 3) * wi[i]
    u4 = 0
    for i in range(len(xi)):
        u4 += pow(xi[i], 4) * wi[i]

    vib_asi = (2*pow(vib_sr, 3) - 3*(pow(vib_sr, 2) + vib_dis)*vib_sr + u3)/pow(vib_kvad,3)
    print("Выборочный коэффициент асимметрии: ", round(vib_asi, 5))
    vib_eks = (u4 - 4 * u3*vib_sr + 6 *(pow(vib_sr, 2) + vib_dis) * pow(vib_sr, 2) - 3 * pow(vib_sr,4))/pow(vib_kvad,4) - 3
    print("Выборочный коэффициент эксцесса: ", round(vib_eks, 5))
    print("-----------------------------")
    print("Теоретические значения: ")
    teor_sr = round((5+v%16) * p, 5)
    print("Математическое ожидание: ", teor_sr)
    teor_dis = round(((5+v%16) * p * (1-p)), 5)
    print("Дисперсия: ", teor_dis)
    teor_kvad = round(pow((5 + v % 16) * p * (1 - p), 0.5), 5)
    print("Среднее квадратическое отклонение: ", teor_kvad)
    teor_mod = 0
    if ((5 + v % 16 + 1) * p)%1 == 0:
        teor_mod = round((5 + v % 16) * p - 0.5, 5)
        print("Мода: ", teor_mod)
    else:
        teor_mod = round(((5 + v % 16 + 1) * p)//1, 5)
        print("Мода: ", teor_mod)
    teor_med = round((5 + v % 16)*p)
    print("Медиана: ", teor_med)
    teor_asi = round(((1-p) - p)/pow((5 + v % 16) * p * (1 - p), 0.5), 5)
    print("Коэффициент асимметии: ", teor_asi)
    teor_eks = round((1-6*p*(1-p))/((5 + v % 16) * p * (1 - p)), 5)
    print("Коэффициент эксцесса: ", teor_eks)
    print("-----------------------------")
    print("Абсолютное отклонение: ")
    print("Выборочное среднее: ", round(abs(round(vib_sr, 5) - teor_sr), 5))
    print("Выборочная дисперсия: ", round(abs(round(vib_dis, 5) - teor_dis), 5))
    print("Выборочное среднее квадратичное отклонение: ", round(abs(round(vib_kvad, 5) - teor_kvad), 5))
    print("Выборочная мода: ", round(abs(round(vib_mod, 5) - teor_mod), 5))
    print("Выборочная медиана: ", round(abs(round(vib_med, 5) - teor_med), 5))
    print("Выборочный коэффициент асимметрии: ", round(abs(round(vib_asi, 5) - teor_asi), 5))
    print("Выборочный коэффициент эксцесса: ", round(abs(round(vib_eks, 5) - teor_eks), 5))
    print("-----------------------------")
    print("Относительное отклонение: ")
    if not teor_sr == 0:
        print("Выборочное среднее: ", (round(abs(round(vib_sr, 5) - teor_sr), 5))/(teor_sr))
    if teor_sr == 0:
        print("-")
    if not teor_dis == 0:
        print("Выборочная дисперсия: ", (round(abs(round(vib_dis, 5) - teor_dis), 5))/(teor_dis))
    if teor_dis == 0:
        print("-")
    if not teor_kvad == 0:
        print("Выборочное среднее квадратичное отклонение: ", (round(abs(round(vib_kvad, 5) - teor_kvad), 5))/(teor_kvad))
    if teor_kvad == 0:
        print("-")
    if not teor_mod == 0:
        print("Выборочная мода: ", (round(abs(round(vib_mod, 5) - teor_mod), 5))/(teor_mod))
    if teor_mod == 0:
        print("-")
    if not teor_med == 0:
        print("Выборочная медиана: ", (round(abs(round(vib_med, 5) - teor_med), 5))/(teor_med))
    if teor_med == 0:
        print("-")
    if not teor_asi == 0:
        print("Выборочный коэффициент асимметрии: ", (round(abs(round(vib_asi, 5) - teor_asi), 5))/(teor_asi))
    if teor_asi == 0:
        print("-")
    if not teor_eks == 0:
        print("Выборочный коэффициент эксцесса: ", (round(abs(round(vib_eks, 5) - teor_eks), 5))/(teor_eks))
    if teor_eks == 0:
        print("-")

    buf = xi[0]

    for i in range(xi[0]):
        xi.insert(0, buf - i - 1)
        wi.insert(0, 0)
        s.insert(0, 0)
    P = []
    for i in range(len(xi)):
        P.append(math.comb(len(xi), i) * pow(p, i) * pow((1 - p), len(xi) - i))
    maxMod = []
    print('x_i :', xi)
    print('n_i :', ni)
    print('w_i :', wi)
    print('S :', s)
    for i in range(len(P)):
       P[i] = round(P[i],5)
    print('p_i :', P)
    for i in range(len(wi)):
        maxMod.append(abs(wi[i] - P[i]))
        print('|w_i - p_i| = ', round(maxMod[i], 5))
    print('Максимальная разность :', round(max(maxMod), 5))

    plt.axis([0, xi[len(xi)-1], 0, 1])
    plt.plot(xi, wi)
    plt.plot(xi, P, color="red")

    subplot = plt.subplot()
    subplot.tick_params(which='both', width=1)
    subplot.tick_params(which='major', length=7)

    subplot.minorticks_on()
    subplot.xaxis.set_major_locator(tick.MultipleLocator(1))
    subplot.xaxis.set_minor_locator(tick.MultipleLocator(0.5))
    subplot.yaxis.set_major_locator(tick.MultipleLocator(0.1))
    subplot.yaxis.set_minor_locator(tick.MultipleLocator(0.05))
    plt.show()

    plt.grid(True)
    plt.axis([0, xi[len(xi)-1]+1, 0, 1])
    for i in range(len(xi)):
        plt.plot([i, i+1], [s[i], s[i]], color="blue", linewidth=2.0)

    subplot = plt.subplot()
    subplot.tick_params(which='both', width=1)
    subplot.tick_params(which='major', length=7)

    subplot.minorticks_on()
    subplot.xaxis.set_major_locator(tick.MultipleLocator(1))
    subplot.xaxis.set_minor_locator(tick.MultipleLocator(0.5))
    subplot.yaxis.set_major_locator(tick.MultipleLocator(0.1))
    subplot.yaxis.set_minor_locator(tick.MultipleLocator(0.05))
    plt.show()

def geometric_law(p):
    sample = sps.geom.rvs(p=p, size=size)
    for i in range(len(sample)):
        sample[i] += -1
    print('Выборка: \n', sample[:size])
    sample.sort()
    print('Выборка сортированная: \n', sample[:size])
    xi = []
    ni = []
    wi = []
    s = []

    xi.append(0)
    ni.append(1)
    wi.append(0)
    s.append(0)

    for i in range(size):
        if not i == size - 1:
            if not sample[i] == sample[i + 1]:
                xi[len(xi) - 1] = sample[i]
                ni[len(ni) - 1] = list(sample).count(sample[i])
                wi[len(wi) - 1] = ni[len(ni) - 1] / size
                s[len(s) - 1] = round(sum(wi), 5)
                ni.append(1)
                wi.append(0)
                s.append(0)
                xi.append(0)
        else:
            xi[len(xi) - 1] = sample[i]
            ni[len(ni) - 1] = list(sample).count(sample[i])
            wi[len(wi) - 1] = ni[len(ni) - 1] / size
            s[len(s) - 1] = round(sum(wi), 5)

    vib_sr = 0
    vib_dis = 0
    vib_kvad = 0
    print("Эксперементальные значения: ")
    # Выборочное среднее
    for i in range(len(xi)):
        vib_sr += xi[i] * wi[i]
    print("Выборочное среднее: ", round(vib_sr, 5))
    # Выборочная дисперсия
    for i in range(len(xi)):
        vib_dis += pow((xi[i] - vib_sr), 2) * wi[i]
    print("Выборочная дисперсия: ", round(vib_dis, 5))
    # Выборочное срднее квадратическое откланение
    vib_kvad = pow(vib_dis, 0.5)
    print("Выборочное срднее квадратическое откланение: ", round(vib_kvad, 5))

    ni_count = 0
    ni_buf = []
    for i in range(len(ni)):
        ni_buf.append(False)

    vib_mod = 0
    index = 0
    # Выборочная мода
    for i in range(len(ni_buf)):
        if max(ni) == ni[i]:
            ni_buf[i] = True
            ni_count += 1
            index = i
        else:
            ni_buf[i] = False

    if ni_count == 1:
        vib_mod = xi[index]
    else:
        for i in range(index - ni_count + 1, index):
            if ni_buf[i] == ni_buf[i + 1]:
                vib_mod = (xi[index - ni_count + 1] + xi[index]) / 2
            else:
                vib_mod = -1
                break
    if not vib_mod == -1:
        print("Выборка по моде: ", round(vib_mod, 5))
    else:
        print("Выборка по моде: Не существует")

    vib_med = 0

    for i in range(len(s)):
        if s[i] > 0.5:
            vib_med = xi[i]
            break
        if s[i] == 0.5:
            vib_med = (xi[i] + xi[i + 1]) / 2
            break
    print("Выборочная медиана: ", round(vib_med, 5))

    u3 = 0
    for i in range(len(xi)):
        u3 += pow(xi[i], 3) * wi[i]
    u4 = 0
    for i in range(len(xi)):
        u4 += pow(xi[i], 4) * wi[i]

    vib_asi = (2 * pow(vib_sr, 3) - 3 * (pow(vib_sr, 2) + vib_dis) * vib_sr + u3) / pow(vib_kvad, 3)
    print("Выборочный коэффициент асимметрии: ", round(vib_asi, 5))
    vib_eks = (u4 - 4 * u3*vib_sr + 6 *(pow(vib_sr, 2) + vib_dis) * pow(vib_sr, 2) - 3 * pow(vib_sr,4))/pow(vib_kvad,4) - 3
    print("Выборочный коэффициент эксцесса: ", round(vib_eks, 5))
    print("-----------------------------")
    print("Теоретические значения: ")
    teor_sr = round((1 - p) / p, 5)
    print("Математическое ожидание: ", teor_sr)
    teor_dis = round((1 - p) / pow(p, 2), 5)
    print("Дисперсия: ", teor_dis)
    teor_kvad = round(pow((1 - p), 0.5) / p, 5)
    print("Среднее квадратическое отклонение: ", teor_kvad)
    teor_mod = 0
    print("Мода: ", teor_mod)
    #if (math.log(2) / math.log(1 - p)) % 1 == 0:
    #    teor_med = round(-(math.log(2) / math.log(1 - p)) - 0.5, 5)
    #    print("Медиана: ", teor_med)
    #else:
    #    teor_med = round((-(math.log(2) / math.log(1 - p))) % 1, 5)
    #    print("Медиана: ", teor_med)
    if (math.log(2) / math.log(1 - p)) % 1 != 0:
        medTeor = (-1 * (math.log(2) / math.log(1 - p))) // 1
    else:
        medTeor = (-1 * (math.log(2) / math.log(1 - p))) - 0.5
    print(medTeor)
    teor_asi = round((2 - p) / pow(1 - p, 0.5), 5)
    print("Коэффициент асимметрии: ", teor_asi)
    teor_eks = round(6 + pow(p, 2) / (1 - p), 5)
    print("Коффициент эксцесса: ", teor_eks)
    print("-----------------------------")
    print("Абсолютное отклонение: ")
    print("Выборочное среднее: ", round(abs(round(vib_sr, 5) - teor_sr), 5))
    print("Выборочная дисперсия: ", round(abs(round(vib_dis, 5) - teor_dis), 5))
    print("Выборочное среднее квадратичное отклонение: ", round(abs(round(vib_kvad, 5) - teor_kvad), 5))
    print("Выборочная мода: ", round(abs(round(vib_mod, 5) - teor_mod), 5))
    print("Выборочная медиана: ", round(abs(round(vib_med, 5) - medTeor), 5))
    print("Выборочный коэффициент асимметрии: ", round(abs(round(vib_asi, 5) - teor_asi), 5))
    print("Выборочный коэффициент эксцесса: ", round(abs(round(vib_eks, 5) - teor_eks), 5))
    print("-----------------------------")
    print("Относительное отклонение: ")
    if not teor_sr == 0:
        print("Выборочное среднее: ", (round(abs(round(vib_sr, 5) - teor_sr), 5)) / (teor_sr))
    if teor_sr == 0:
        print("-")
    if not teor_dis == 0:
        print("Выборочная дисперсия: ", (round(abs(round(vib_dis, 5) - teor_dis), 5)) / (teor_dis))
    if teor_dis == 0:
        print("-")
    if not teor_kvad == 0:
        print("Выборочное среднее квадратичное отклонение: ",
              (round(abs(round(vib_kvad, 5) - teor_kvad), 5)) / (teor_kvad))
    if teor_kvad == 0:
        print("-")
    if not teor_mod == 0:
        print("Выборочная мода: ", (round(abs(round(vib_mod, 5) - teor_mod), 5)) / (teor_mod))
    if teor_mod == 0:
        print("-")
    if not medTeor == 0:
        print("Выборочная медиана: ", (round(abs(round(vib_med, 5) - medTeor), 5)) / (medTeor))
    if medTeor == 0:
        print("-")
    if not teor_asi == 0:
        print("Выборочный коэффициент асимметрии: ", (round(abs(round(vib_asi, 5) - teor_asi), 5)) / (teor_asi))
    if teor_asi == 0:
        print("-")
    if not teor_eks == 0:
        print("Выборочный коэффициент эксцесса: ", (round(abs(round(vib_eks, 5) - teor_eks), 5)) / (teor_eks))
    if teor_eks == 0:
        print("-")

    buf = xi[0]

    for i in range(xi[0]):
        xi.insert(0, buf - i - 1)
        wi.insert(0, 0)
        s.insert(0, 0)
    P = []
    for i in range(len(xi)):
        P.append((pow((1-p), i)) * p)
    maxMod = []
    print('x_i :', xi)
    print('n_i :', ni)
    print('w_i :', wi)
    print('S :', s)
    for i in range(len(P)):
        P[i] = round(P[i], 5)
    print('p_i :', P)
    for i in range(len(wi)):
        maxMod.append(abs(wi[i]-P[i]))
        print('|w_i - p_i| = ', round(maxMod[i],5))
    print('Максимальная разность :', round(max(maxMod),5))

    plt.axis([0, xi[len(xi) - 1], 0, 1])
    plt.plot(xi, wi)
    plt.plot(xi, P, color="red")

    subplot = plt.subplot()
    subplot.tick_params(which='both', width=1)
    subplot.tick_params(which='major', length=7)

    subplot.minorticks_on()
    subplot.xaxis.set_major_locator(tick.MultipleLocator(1))
    subplot.xaxis.set_minor_locator(tick.MultipleLocator(0.5))
    subplot.yaxis.set_major_locator(tick.MultipleLocator(0.1))
    subplot.yaxis.set_minor_locator(tick.MultipleLocator(0.05))
    plt.show()

    plt.grid(True)
    plt.axis([0, xi[len(xi) - 1] + 1, 0, 1])
    for i in range(len(xi)):
        plt.plot([i, i + 1], [s[i], s[i]], color="blue", linewidth=2.0)

    subplot = plt.subplot()
    subplot.tick_params(which='both', width=1)
    subplot.tick_params(which='major', length=7)

    subplot.minorticks_on()
    subplot.xaxis.set_major_locator(tick.MultipleLocator(1))
    subplot.xaxis.set_minor_locator(tick.MultipleLocator(0.5))
    subplot.yaxis.set_major_locator(tick.MultipleLocator(0.1))
    subplot.yaxis.set_minor_locator(tick.MultipleLocator(0.05))
    plt.show()

def poisson_law(p):
    sample = sps.poisson.rvs(mu=lamb, size=size)
    print('Выборка: \n', sample[:size])
    sample.sort()
    print('Выборка сортированная: \n', sample[:size])

    xi = []
    ni = []
    wi = []
    s = []

    xi.append(0)
    ni.append(1)
    wi.append(0)
    s.append(0)

    for i in range(size):
        if not i == size - 1:
            if not sample[i] == sample[i + 1]:
                xi[len(xi) - 1] = sample[i]
                ni[len(ni) - 1] = list(sample).count(sample[i])
                wi[len(wi) - 1] = ni[len(ni) - 1] / size
                s[len(s) - 1] = round(sum(wi), 5)
                ni.append(1)
                wi.append(0)
                s.append(0)
                xi.append(0)
        else:
            xi[len(xi) - 1] = sample[i]
            ni[len(ni) - 1] = list(sample).count(sample[i])
            wi[len(wi) - 1] = ni[len(ni) - 1] / size
            s[len(s) - 1] = round(sum(wi), 5)

    vib_sr = 0
    vib_dis = 0
    vib_kvad = 0
    print("Эксперементальные значения: ")
    # Выборочное среднее
    for i in range(len(xi)):
        vib_sr += xi[i] * wi[i]
    print("Выборочное среднее: ", round(vib_sr, 5))
    # Выборочная дисперсия
    for i in range(len(xi)):
        vib_dis += pow((xi[i] - vib_sr), 2) * wi[i]
    print("Выборочная дисперсия: ", round(vib_dis, 5))
    # Выборочное срднее квадратическое откланение
    vib_kvad = pow(vib_dis, 0.5)
    print("Выборочное срднее квадратическое откланение: ", round(vib_kvad, 5))
    ni_count = 0
    ni_buf = []
    for i in range(len(ni)):
        ni_buf.append(False)

    vib_mod = 0
    index = 0
    # Выборочная мода
    for i in range(len(ni_buf)):
        if max(ni) == ni[i]:
            ni_buf[i] = True
            ni_count += 1
            index = i
        else:
            ni_buf[i] = False

    if ni_count == 1:
        vib_mod = xi[index]
    else:
        for i in range(index - ni_count + 1, index):
            if ni_buf[i] == ni_buf[i + 1]:
                vib_mod = round(((xi[index - ni_count + 1] + xi[index]) / 2), 5)
            else:
                vib_mod = -1
                break
    if not vib_mod == -1:
        print("Выборка по моде: ", round(vib_mod, 5))
    else:
        print("Выборка по моде: Не существует")

    vib_med = 0

    for i in range(len(s)):
        if s[i] > 0.5:
            vib_med = xi[i]
            break
        if s[i] == 0.5:
            vib_med = (xi[i] + xi[i + 1]) / 2
            break
    print("Выборочная медиана: ", round(vib_med, 5))

    u3 = 0
    for i in range(len(xi)):
        u3 += pow(xi[i], 3) * wi[i]
    u4 = 0
    for i in range(len(xi)):
        u4 += pow(xi[i], 4) * wi[i]

    vib_asi = (2 * pow(vib_sr, 3) - 3 * (pow(vib_sr, 2) + vib_dis) * vib_sr + u3) / pow(vib_kvad, 3)
    print("Выборочный коэффициент асимметрии: ", round(vib_asi, 5))
    vib_eks = (u4 - 4 * u3*vib_sr + 6 *(pow(vib_sr, 2) + vib_dis) * pow(vib_sr, 2) - 3 * pow(vib_sr,4))/pow(vib_kvad,4) - 3
    print("Выборочный коэффициент эксцесса: ", round(vib_eks, 5))
    print("-----------------------------")
    print("Теоретические значения: ")
    print("Математическое ожидание: ", lamb)
    print("Дисперсия: ", lamb)
    print("Среднее квадратическое отклонение: ", round(pow(lamb, 0.5), 5))
    print("Мода: ", round(lamb//1, 5))
    teor_med = round((lamb+1/3-0.02/lamb)//1, 5)
    print("Медиана: ", teor_med)
    print("Коэффициент асимметрии: ", round(pow(lamb,-0.5), 5))
    print("Коффициент эксцесса: ", round(pow(lamb, -1), 5))
    print("-----------------------------")
    print("Абсолютное отклонение: ")
    print("Выборочное среднее: ", round(abs(round(vib_sr, 5) - lamb), 5))
    print("Выборочная дисперсия: ", round(abs(round(vib_dis, 5) - lamb), 5))
    print("Выборочное среднее квадратичное отклонение: ", round(abs(round(vib_kvad, 5) - round(pow(lamb, 0.5), 5)), 5))
    print("Выборочная мода: ", round(abs(round(vib_mod, 5) - round(lamb//1, 5)), 5))
    print("Выборочная медиана: ", round(abs(round(vib_med, 5) - teor_med), 5))
    print("Выборочный коэффициент асимметрии: ", round(abs(round(vib_asi, 5) - round(pow(lamb,-0.5), 5)), 5))
    print("Выборочный коэффициент эксцесса: ", round(abs(round(vib_eks, 5) - round(pow(lamb, -1), 5)), 5))
    print("-----------------------------")
    print("Относительное отклонение: ")
    print("Выборочное среднее: ", round((round(abs(round(vib_sr, 5) - lamb), 5)/lamb), 5))
    print("Выборочная дисперсия: ", round((round(abs(round(vib_dis, 5) - lamb), 5)/lamb), 5))
    print("Выборочное среднее квадратичное отклонение: ", round((round(abs(round(vib_kvad, 5) - round(pow(lamb, 0.5), 5)), 5)/round(pow(lamb, 0.5), 5)), 5))
    print("Выборочная мода: ", round((round(abs(round(vib_mod, 5) - round(lamb//1, 5)), 5)/round(lamb//1, 5)), 5))
    if not teor_med == 0:
        print("Выборочная медиана: ", round((round(abs(round(vib_med, 5) - teor_med), 5)/teor_med), 5))
    if teor_med == 0:
        print("-")
    print("Выборочный коэффициент асимметрии: ", round((round(abs(round(vib_asi, 5) - round(pow(lamb,-0.5), 5)), 5)/round(pow(lamb,-0.5), 5)), 5))
    print("Выборочный коэффициент эксцесса: ", round((round(abs(round(vib_eks, 5) - round(pow(lamb, -1), 5)), 5)/round(pow(lamb, -1), 5)), 5))

    buf = xi[0]

    for i in range(xi[0]):
        xi.insert(0, buf - i - 1)
        wi.insert(0, 0)
        s.insert(0, 0)
    P = []
    for i in range(len(xi)):
        P.append((pow(lamb, i) * math.exp(-lamb))/math.factorial(i))

    maxMod = []
    print('x_i :', xi)
    print('n_i :', ni)
    print('w_i :', wi)
    print('S :', s)
    for i in range(len(P)):
        P[i] = round(P[i], 5)
    print('p_i :', P)
    for i in range(len(wi)):
        maxMod.append(abs(wi[i] - P[i]))
        print('|w_i - p_i| = ', round(maxMod[i],5))
    print('Максимальная разность :', round(max(maxMod),5))

    plt.axis([0, xi[len(xi) - 1], 0, 1])
    plt.plot(xi, wi)
    plt.plot(xi, P, color="red")

    subplot = plt.subplot()
    subplot.tick_params(which='both', width=1)
    subplot.tick_params(which='major', length=7)

    subplot.minorticks_on()
    subplot.xaxis.set_major_locator(tick.MultipleLocator(1))
    subplot.xaxis.set_minor_locator(tick.MultipleLocator(0.5))
    subplot.yaxis.set_major_locator(tick.MultipleLocator(0.1))
    subplot.yaxis.set_minor_locator(tick.MultipleLocator(0.05))
    plt.show()

    plt.grid(True)
    plt.axis([0, xi[len(xi) - 1] + 1, 0, 1])
    for i in range(len(xi)):
        plt.plot([i, i + 1], [s[i], s[i]], color="blue", linewidth=2.0)

    subplot = plt.subplot()
    subplot.tick_params(which='both', width=1)
    subplot.tick_params(which='major', length=7)

    subplot.minorticks_on()
    subplot.xaxis.set_major_locator(tick.MultipleLocator(1))
    subplot.xaxis.set_minor_locator(tick.MultipleLocator(0.5))
    subplot.yaxis.set_major_locator(tick.MultipleLocator(0.1))
    subplot.yaxis.set_minor_locator(tick.MultipleLocator(0.05))
    plt.show()

binomial_law(p)
geometric_law(p)
poisson_law(p)

