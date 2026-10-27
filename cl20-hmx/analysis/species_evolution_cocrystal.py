# counting species
#import matplotlib.pyplot as plt
#import statsmodels.api as sm

# species and isomers ====要改===
H2O = ['[H][O][H]']
CO2 = ['[O]=[C][O]', '[O]=[C]=[O]']
N2 = ['[N]#[N]']
NO2 = ['[O]=[N][O]', '[O][N][O]']
NO = ['[N]=[O]']
#ONNO = ['[O][N][N][O]', '[O]=[N][N][O]', '[O]=[N][N]=[O]']
HONO = ['[H][O][N][O]', '[H][O][N]=[O]']
HNCO = ['[H][N]=[C]=[O]', '[H][N][C]=[O]', '[H][N][C][O]']
N2O = ['[N][N][O]']
# HO = ['[H][O]']
HMX = ['[H][C]1([H])[N]([N](=[O])=[O])[C]([H])([H])[N]([N](=[O])=[O])[C]([H])([H])[N]([N](=[O])=[O])[C]([H])([H])[N]1[N](=[O])=[O]']
CL20 = ['[H][C]12[N]([N](=[O])=[O])[C]3([H])[N]([N](=[O])=[O])[C]1([H])[N]([N](=[O])=[O])[C]1([H])[N]([N](=[O])=[O])[C]3([H])[N]([N](=[O])=[O])[C]1([H])[N]2[N](=[O])=[O]']
# List of species and labels===要改===
# initial decomposition
#TOTAL = [NO2, CH3, CH3NO, CH2NO2, CH4, NO, OH]
# LABEL = ['NM', 'NO2', 'CH3', 'CH3NOOH', 'CH2NO2']

# CO2 generation
# TOTAL = [NM, CH3OH, CH3O, HCHO, CO, CO2]
# LABEL = ['NM', 'CH3OH', 'CH3O', 'HCHO', 'CO', 'CO2']

# H2O generation
# TOTAL = [NM, HONO, CH3NOOH, H2O, OH]
# LABEL = ['NM', 'HONO', 'CH3NOOH', 'H2O', 'OH']

# N2 generation
# TOTAL = [NM, NO, ONNO, N2, N2O]
# LABEL = ['NM', 'NO', 'ONNO', 'N2', 'N2O']

# CH4 generation
# TOTAL = [NM, CH3OH, CH3O, CH2O, CH4, C2H6, CH2]
# LABEL = ['NM', 'CH3OH', 'CH3O', 'CH2O', 'CH4', 'C2H6', 'CH2']

# CN generation
# TOTAL = [NM, CN, HCN, HCNH, NCNCN]
# LABEL = ['NM', 'CN', 'HCN', 'HCNH', 'NCNCN']

# HNCO generation
# TOTAL = [NM, HNCO, NCO, CH2N]
# LABEL = ['NM', 'HNCO', 'NCO', 'CH2N']

# product
TOTAL = [CL20, HMX, NO2, NO, HONO, HNCO, N2O, N2, H2O, CO2]
#LABEL = ['NM', 'H2O', 'CO2', 'N2', 'HNCO', 'CO', 'NH3', 'NH4', 'H2', 'HCN', 'CN']
#TOTAL = [NCN]
#TOTAL = [HNC, NCN, HCNH, NCO, NO2, CH3, CH3NOOH, CH2NO2, CH3NO]
#LABEL = ['CH3', 'NO2', 'CH2NO2', 'CH3NO', 'CH4', 'NO', 'OH', 'HCHO', 'HCNO']


# basic settings
stepsize = 0.0001  # unit: ps  ===可改===
f = open("dump.lammpstrj.species",'r')  # read file ===可改===

# initiation
time = []
TOTAl_num = []
for specie in TOTAL:
    TOTAl_num.append([]) # each specie has a list to store the number at each time step
byt = f.readlines()

# counting species step by step
for line in byt:
    curline = line.strip().split()
    time.append(round(stepsize*float(curline[1].replace(':','')),4))   # save time as x-axis

    for specie in TOTAL:
        j = TOTAL.index(specie)
        TOTAl_num[j].append(int(0))
        for b in specie:
            for a in curline:
                i = curline.index(a)
                if a == b:
                    TOTAl_num[j][-1] = TOTAl_num[j][-1] + int(curline[i+1])

# Save data
f1 = open('data.txt','w')
for i in range(len(time)):
    f1.write(str('%.4f'%(time[i]))+'\t')
    for j in range(len(TOTAL)):
        f1.write(str(TOTAl_num[j][i])+'\t')
    f1.write('\n')
f1.close()

# Plot specie
# lowess=sm.nonparametric.lowess  # fitting setteings
# for i in range(len(TOTAl_num)):
    # plt.plot(time, TOTAl_num[i], lw=0.5, label=LABEL[i])
    # z = lowess(TOTAl_num[i],time,frac=0.1)
    # #plt.plot(z[:,0],z[:,1],lw=1, color='r')
    # plt.plot(z[:, 0], z[:, 1], lw=2, label=LABEL[i])
# plt.xlabel('time(ps)')
# plt.ylabel('number of species')
# plt.legend(bbox_to_anchor=(1.05, 0), loc=3)
# # plt.ylim([-0.05, 0.5])
# plt.savefig('species evolution.png',dpi=600, bbox_inches='tight')
# plt.show()
f.close()




