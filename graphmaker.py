import pandas as pd
from matplotlib import pyplot as plt
import numpy
import pylab
# reads the csv file
data = pd.read_csv('ex_csv.csv')
# selects the rows you want (right now it is selecting all rows with promoter size = 500)
data = data.iloc[110:163]
# makes a column go in ascending order rather than random
data = data.sort_values([' number_of_sites'], ascending=True)
# pulls out only the data points you are interested in
print(data[[' motif', ' number_of_sites', ' e_val', ' distance', ' success_rate',
         ' promoter_size']])

success = data[' success_rate']
numsite = data[' number_of_sites']

# plots your variable on x-y-axis and "o" will be the trendline
plt.plot(numsite, success, "o")

# calc the trendline (copied code from a stack overflow answer)
z = numpy.polyfit(numsite, success, 1)
p = numpy.poly1d(z)
pylab.plot(numsite,p(numsite),"r--")
# add labels
plt.xlabel('number of sites')
plt.ylabel('success rate')
plt.show()

