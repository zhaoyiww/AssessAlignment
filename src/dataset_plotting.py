from timeit import default_timer as time
import numpy as np
import matplotlib.pylab as plt
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rc('font', family="Times New Roman", size=20)

time_0 = time()
print('Start running ...')

mAP_5 = np.array([[0.88, 0.89, 0.97, 0.93], [0.93, 0.93, 0.99, 0.99]])
mAP_5_9 = np.array([[0.64, 0.64, 0.76, 0.70], [0.77, 0.79, 0.85, 0.87]])

plt.figure(figsize=(8, 6))
fig, ax1 = plt.subplots()
t = [1, 2]
colors = ['r', 'b']
symbols = ['^', '*', 'o', '+']
for i in range(len(t)):
    if i == 0:
        plt.scatter(t[i] - 0.2, mAP_5[i, 0], c='r', marker=symbols[0], label='mAP@0.5 F R-CNN', s=150)
        plt.scatter(t[i] - 0.1, mAP_5[i, 1], c='r', marker=symbols[1], label='mAP@0.5 M R-CNN', s=150)
        plt.scatter(t[i] + 0.1, mAP_5[i, 2], c='r', marker=symbols[2], label='mAP@0.5 YOLOX', s=150)
        plt.scatter(t[i] + 0.2, mAP_5[i, 3], c='r', marker=symbols[3], label='mAP@0.5 YOLOv7', s=150)

        plt.scatter(t[i] - 0.2, mAP_5_9[i, 0], c='b', marker=symbols[0], label='mAP@0.5:0.95 F R-CNN', s=150)
        plt.scatter(t[i] - 0.1, mAP_5_9[i, 1], c='b', marker=symbols[1], label='mAP@0.5:0.95 M R-CNN', s=150)
        plt.scatter(t[i] + 0.1, mAP_5_9[i, 2], c='b', marker=symbols[2], label='mAP@0.5:0.95 YOLOX', s=150)
        plt.scatter(t[i] + 0.2, mAP_5_9[i, 3], c='b', marker=symbols[3], label='mAP@0.5:0.95 YOLOv7', s=150)
    else:
        a1 = plt.scatter(t[i] - 0.2, mAP_5[i, 0], c='r', marker=symbols[0], s=150)
        a2 = plt.scatter(t[i] - 0.1, mAP_5[i, 1], c='r', marker=symbols[1], s=150)
        a3 = plt.scatter(t[i] + 0.1, mAP_5[i, 2], c='r', marker=symbols[2], s=150)
        a4 = plt.scatter(t[i] + 0.2, mAP_5[i, 3], c='r', marker=symbols[3], s=150)

        b1 = plt.scatter(t[i] - 0.2, mAP_5_9[i, 0], c='b', marker=symbols[0], s=150)
        b2 = plt.scatter(t[i] - 0.1, mAP_5_9[i, 1], c='b', marker=symbols[1], s=150)
        b3 = plt.scatter(t[i] + 0.1, mAP_5_9[i, 2], c='b', marker=symbols[2], s=150)
        b4 = plt.scatter(t[i] + 0.2, mAP_5_9[i, 3], c='b', marker=symbols[3], s=150)

ax1.set_xticks([1, 2])
# ax1.set_xticklabels(('Target_R', 'Target_S', 'Target_C'))
ax1.set_xticklabels(('Target_R', 'Target_C'))
ax1.set_yticks([0.45, 0.5, 0.625, 0.75, 0.875, 1.0, 1.05])
ax1.set_yticklabels(('', '0.50', '', '0.75', '', '1.0', ''))
plt.grid()

from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
plt.legend([(a1, a2, a3, a4), (b1, b2, b3, b4)], [('mAP@.5 F R-CNN, M R-CNN, YOLOX, YOLOv7'),
                                                  ('mAP@[.5, .95] F R-CNN, M R-CNN, YOLOX, YOLOv7')],
           handler_map={tuple: HandlerTuple(ndivide=None)},
           loc=8, bbox_to_anchor=(0.5, -0.35), borderaxespad=0.)

fig.savefig('datasets.pdf', dpi=600, bbox_inches='tight')
plt.close()

print(f'Total running time: {time() - time_0}')
print('All done!')
