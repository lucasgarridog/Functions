filename = r"C:\Users\Lucas Garrido\Documents\PhD\Documentos\2023 Brasil\Datos_X64Zn\7lizn_e18p14.dat"

import numpy as np

lucas = np.loadtxt(filename)
ind = np.argsort(lucas[:,0])
lucas = lucas[ind]
theta = []
sigma = []
error = []
for f in lucas:
    theta.append(np.round(f[0],1))
    sigma.append(np.round(f[1],3))
    error.append(np.round(f[2],3))

caca = int(len(theta)/2) # si el n de lineas del txt es impar pon +1 aqui
print("\\begin{table}[htbp]\n\\centering\n\\begin{tabular}{cccc}\n$\\theta_{CM}$ ($^\\circ$) & $\\sigma_{el}/\\sigma_R$ & $\\theta_{CM}$ ($^\\circ$) & $\\sigma_{el}/\\sigma_R$ \\bigstrut[b]\\\\\n\\hline")
for i in range(caca):
    if i == 0:
        print(theta[i], "&", "%.3f " % sigma[i], "$ \pm $", "%.3f " % error[i], "&", theta[i + caca], "&",
              "%.3f " % sigma[i + caca], "$ \pm $", "%.3f " % error[i + caca], "\\bigstrut[t] \\\\")
    else:
        print(theta[i], "&", "%.3f " % sigma[i], "$ \pm $", "%.3f " % error[i], "&", theta[i+caca], "&", "%.3f " % sigma[i+caca], "$ \pm $", "%.3f " % error[i+caca], "\\\\")
print("\\hline\n\\end{tabular}\n\\caption{Datos experimentales correspondientes a la reacción \\ce{^8B} + \\ce{^64Zn} con energía incidente $E_{CM}=$ 34.22 MeV.}\n\\end{table}")