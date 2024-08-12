import numpy as np
import matplotlib.pyplot as plt

with open(r"Temperatures_bf.txt") as inp:
    T_bf = np.array(inp.read().split()).astype(float)

with open(r"Temperatures_af.txt") as inp:
    T_af = np.array(inp.read().split()).astype(float)

with open(r"Energies_bf.txt") as inp:
    E_bf = np.array(inp.read().split()).astype(float)

with open(r"Energies_af.txt") as inp:
    E_af = np.array(inp.read().split()).astype(float)

with open(r"Settings.txt") as inp:
    Settings = np.array(inp.read().split()).astype(float)


N = int(Settings[0])
iterations = int(Settings[1])

MC_STEPS = int(T_bf.size/(N*iterations))
L = int(T_bf.size/N)



probabilities_bf = np.zeros((L, N-1))
probabilities_af = np.zeros((L, N-1))



E_bf_tmp = E_bf.reshape((L, N)).copy()
E_af_tmp = E_af.reshape((L, N)).copy()

T_bf_tmp = T_bf.reshape((L, N)).copy()
T_af_tmp = T_af.reshape((L, N)).copy()

for i in range(L):
    for j in range(N-1):
        probabilities_bf[i][j] = '{:10.3f}'.format(np.exp((E_bf_tmp[i][j+1]-E_bf_tmp[i][j])*(1/T_bf_tmp[i][j+1]-1/T_bf_tmp[i][j])))
        probabilities_af[i][j] = '{:10.3f}'.format(np.exp((E_af_tmp[i][j+1]-E_af_tmp[i][j])*(1/T_af_tmp[i][j+1]-1/T_af_tmp[i][j])))

probabilities_bf = probabilities_bf.flatten()
probabilities_af = probabilities_af.flatten()

probabilities_bf *= 100
probabilities_af *= 100



for i in range(L):
    Figure, axis = plt.subplots(1, 2, figsize=(15, 6))
    axis[0].set_title(f"MC step: {i//iterations+1}/{MC_STEPS}, iteration: {i%iterations+1}/{iterations}.  "
              f"\n probabilities before: {probabilities_bf[(N-1)*i:(N-1)*(i+1)]} %, p = {'{:10.3f}'.format(probabilities_bf[(N-1)*i:(N-1)*(i+1)].mean())} %")
    axis[0].set_xlabel('Temperatures')
    axis[0].set_ylabel('Energies')
    [axis[0].text(T_bf[i*N+j], E_bf[i*N+j], j+1) for j in range(N)]
    axis[0].grid(True)
    axis[0].scatter(T_bf[N*i:N*(i+1)], E_bf[N*i:N*(i+1)], marker='o', label = 'before')



    axis[1].set_title(f"MC step: {i//iterations+1}/{MC_STEPS}, iteration: {i%iterations+1}/{iterations}.  "
              f"\n probabilities after: {probabilities_af[(N-1)*i:(N-1)*(i+1)]} %, p = {'{:10.3f}'.format(probabilities_af[(N-1)*i:(N-1)*(i+1)].mean())} %")
    axis[1].set_xlabel('Temperatures')
    axis[1].set_ylabel('Energies')
    [axis[1].text(T_af[i*N+j], E_af[i*N+j], j+1) for j in range(N)]
    axis[1].grid(True)
    axis[1].scatter(T_af[N*i:N*(i+1)], E_af[N*i:N*(i+1)], marker='o', label = 'after')

    plt.savefig(f"stepsOfNewApproach/{i+1}_parallelTempering_steps.jpg")