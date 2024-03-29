
Consider a series-connected multijunction device. The voltage across cell $i$ is

$V_i(I)=\frac{n_ikT}{q}  \cdot \ln(\frac{I_{p_i} - I}{I_{0_1}} + 1)$

The total power is

$P = \sum_i I \cdot V_i(I)$

If we have two cells, write

$I_{p,1}=(1-f)\cdot I_p$

$I_{p,2}=f\cdot I_p$

$P = \frac{n_1kT}{q}  \cdot I \cdot \ln(\frac{(1-f)\cdot I_p - I}{I_{0,1}} + 1) + \frac{n_2kT}{q}  \cdot I \cdot \ln(\frac{f\cdot I_p - I}{I_{0,2}} + 1)$

Find $f$ for maximum power:

$\frac{\partial P}{\partial f} = 0$

$0 = -\frac{n_1kT}{q}  \cdot I \cdot (\frac{I_p}{(1-f)\cdot I_p - I+I_{0,1}})\
+\frac{n_2kT}{q}  \cdot I \cdot (\frac{I_p}{f\cdot I_p - I+I_{0,2}})$

Find $f$ and $1-f$:

$f=(\frac{n_2}{n_1+n_2}) + (\frac{n_2-n_1}{n_1+n_2})\cdot\frac{I}{I_p}\
+(\frac{n_2}{n_1+n_2})\cdot \frac{I_{0,1}}{I_p}\
-(\frac{n_1}{n_1+n_2})\cdot \frac{I_{0,2}}{I_p}$

$1-f=(\frac{n_1}{n_1+n_2}) - (\frac{n_2-n_1}{n_1+n_2})\cdot\frac{I}{I_p}\
-(\frac{n_2}{n_1+n_2})\cdot \frac{I_{0,1}}{I_p}\
+(\frac{n_1}{n_1+n_2})\cdot \frac{I_{0,2}}{I_p}$

The photocurrents are

$I_{p,1}=(\frac{n_2}{n_1+n_2}) \cdot I_p + (\frac{n_1-n_2}{n_1+n_2})\cdot I +\
(\frac{n_2}{n_1+n_2})\cdot I_{0,1}\
-(\frac{n_1}{n_1+n_2})\cdot I_{0,2}$

$I_{p,2}=(\frac{n_1}{n_1+n_2}) \cdot I_p - (\frac{n_1-n_2}{n_1+n_2})\cdot I -\
(\frac{n_2}{n_1+n_2})\cdot I_{0,1}\
+(\frac{n_1}{n_1+n_2})\cdot I_{0,2}$

