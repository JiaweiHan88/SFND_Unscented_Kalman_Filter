import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data_r = pd.read_csv('../build/NIS_RADAR.txt', sep=",")
data_r.head()
plt.figure(figsize=(10,6))
plt.plot(data_r, label='NIS RADAR Data')
plt.plot([0,len(data_r)],[7.815,7.815],'r--',lw=2, label='Chi square = 7.815 for 3 DOF')
plt.xlabel('x')
plt.ylabel('y')
plt.title('NIS Radar')
plt.legend()
plt.show()