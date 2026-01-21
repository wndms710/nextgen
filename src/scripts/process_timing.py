import sys
import os
import glob
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from collections import defaultdict

colors=['red', 'blue', 'green', 'purple', 'black']
colors_legend=['rhs','rk', 'ghost', 'sgs_mat', 'report', '']

def fmt(x):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}.10^{{{}}}$'.format(a, b)

def plot_percentages(p, g_ncores, test_tag):
    fig, ax = plt.subplots()
    n_groups=len(g_ncores)
    index = np.arange(n_groups)
    bar_width=0.05
    for i in range(0,12):
        if i>=1 and i<=9:
            color_idx=(i+2)%3
            legend_idx=(i+2)%3
            if i>3:
                legend_idx=5
        elif i==10:
            color_idx=3
            legend_idx=color_idx
        elif i==11:
            color_idx=4
            legend_idx=color_idx
        else:
            color_idx=4
            legend_idx=5

        plt.bar(index+bar_width*float(i), [p[key][i] for key in g_ncores], bar_width, alpha=0.5, color=colors[color_idx], label=colors_legend[legend_idx])

    plt.ylabel('% of Time')
    plt.xlabel('#Cores')
    plt.xticks(index + bar_width, g_ncores)
    plt.legend(loc=2,prop={'size':6})
    plt.tight_layout()    
    plt.savefig('percentages_%s.pdf' %(test_tag))

g_speed_avg=[]
g_speed_std=[]
g_ncores=[]
g_cv_core=[]
percentages=defaultdict(list)

test_tag=sys.argv[1]

for outfile in glob.glob('%s.o*' %(test_tag)):
    in_file = open(outfile, 'r')
    ncores=-1
    cv_core=-1
    norm_speed=[]
    time_step=[]
    for line in in_file.readlines():
        if 'cores,' in line:
            columns=line.split()
            ncores=int(columns[3])
        if 'normalized' in line:
            columns=line.split()
            norm_speed.append(float(columns[9]))
            cv_core=float(columns[16])
        if 'since' in line:
            columns=line.split()
            time_step.append(float(columns[5]))
        if 'percent of time' in line:
            columns=line.split()
            percentages[ncores].append(float(columns[6]))

    if cv_core>0:
        if len(norm_speed)>0:
            avg_speed=np.average(norm_speed)
            std_speed=np.std(norm_speed)

        print(ncores, cv_core, avg_speed, std_speed)

        if cv_core in g_cv_core:
            idx=g_cv_core.index(cv_core)
            g_ncores[idx]=ncores
            #take the best
            g_speed_avg[idx]=min(g_speed_avg[idx],avg_speed)
            g_speed_std[idx]=min(g_speed_std[idx],std_speed)
        else:
            g_ncores.append(ncores)
            g_cv_core.append(cv_core)
            g_speed_avg.append(avg_speed)
            g_speed_std.append(std_speed)

    in_file.close()


ideal_speedup=[]
for i in range(0,len(g_ncores)):
    ideal_speedup.append(g_ncores[i]/np.min(g_ncores))

g_cv_core_sort = [x for (y,x) in sorted(zip(ideal_speedup,g_cv_core))]
g_speed_avg_sort = [y for (y,x) in sorted(zip(g_speed_avg,g_cv_core))]
g_speed_std_sort = [y for (y,x) in sorted(zip(g_speed_std,g_cv_core))]
g_ncores_sort = [x for (y,x) in sorted(zip(ideal_speedup,g_cv_core))]

print(g_cv_core_sort, g_speed_avg_sort, g_speed_std_sort)

plt.figure()
plt.plot(g_cv_core_sort, g_speed_avg_sort, ls='-', marker='s')
plt.ylim([1,6])
plt.xlim([1000,1000000])
plt.ylabel('Normalized Speed', fontsize=18)
plt.xlabel('#CVs/core', fontsize=18)
plt.gca().invert_xaxis()
plt.xscale('log')
plt.tick_params(axis='both', which='major', labelsize=16)
plt.savefig('normalized_speed_%s.pdf' %(test_tag))

cur_speedup=[]
for i in range(0,len(g_ncores)):
    cur_speedup.append((g_ncores[i]/np.min(g_ncores))*(np.min(g_speed_avg)/g_speed_avg[i]))

g_ncores_sort = [x for (y,x) in sorted(zip(ideal_speedup,g_ncores))]
ideal_speedup_sort = [y for (y,x) in sorted(zip(ideal_speedup,g_ncores))]
cur_speedup_sort   = [y for (y,x) in sorted(zip(cur_speedup,g_ncores))]
g_cv_core_sort = [x for (y,x) in sorted(zip(ideal_speedup,g_cv_core))]

print(g_ncores_sort, ideal_speedup_sort, cur_speedup_sort)

plt.figure()
plt.plot(g_cv_core_sort, cur_speedup_sort, marker='s')
plt.plot(g_cv_core_sort, ideal_speedup_sort, ls='--')
plt.ylabel('Speedup', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel('#CVs/core', fontsize=18)
plt.xscale('log')
plt.yscale('log')
plt.xlim([1000,1000000])
plt.gca().invert_xaxis()

plt.savefig('speedup_%s.pdf' %(test_tag))

plot_percentages(percentages, g_ncores_sort, test_tag)
