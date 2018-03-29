#coding=utf-8
#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf  import PdfPages

pp=PdfPages('cumulative.pdf')

data=np.genfromtxt ("cumulative.txt")

index=data[:,1:350]
value=data[:,2:350]
#print len(index), len(value), index, value

x_vec,y_vec=[],[]
for i,j in zip (index, value):
    x_vec.append(i[0])
    y_vec.append(i[1])



fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot(x_vec, y_vec, 'r.')
pp.savefig()
pp.close()
