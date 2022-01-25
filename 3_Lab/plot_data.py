import matplotlib.pyplot as plt

explicit_file = open('output/Explicit.txt')
explicit_data = [[float(x) for x in line.split()] for line in explicit_file]
implicit_file = open('output/Implicit.txt')
implicit_data = [[float(x) for x in line.split()] for line in implicit_file]

exp_x_data = [arr_el[0] for arr_el in explicit_data]
exp_u1_data = [arr_el[1] for arr_el in explicit_data]
exp_u2_data = [arr_el[2] for arr_el in explicit_data]

imp_x_data = [arr_el[0] for arr_el in implicit_data]
imp_u1_data = [arr_el[1] for arr_el in implicit_data]
imp_u2_data = [arr_el[2] for arr_el in implicit_data]

plt.plot(exp_x_data, exp_u1_data, label='explicit u1')  # label=r"$u_1^{explicit} $")
# plt.hold(True)
plt.plot(exp_x_data, exp_u2_data, label='explicit u2')
# plt.hold(True)
plt.plot(imp_x_data, imp_u1_data, label='implicit u1')
# plt.hold(True)
plt.plot(imp_x_data, imp_u2_data, label='implicit u2')

plt.legend()
plt.title("Solution for diff. equations obtained by exp. and impl. Euler method")
plt.show()

# plt.pause(10)

# pauseE = input('pause ')
# plt.close()
