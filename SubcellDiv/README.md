
Consider a series-connected multijunction device. The voltage across cell $i$ is

&ensp; $V_i(I)=\frac{n_ikT}{q}  \cdot \ln(\frac{I_{p_i} - I}{I_{0_i}} + 1)$

The total power is

&ensp; $P = \sum_i I \cdot V_i(I)$

If we have two cells, write

&ensp; $I_{p,1}=(1-f)\cdot I_p$

&ensp; $I_{p,2}=f\cdot I_p$

&ensp; $P = \frac{n_1kT}{q}  \cdot I \cdot \ln(\frac{(1-f)\cdot I_p - I}{I_{0,1}} + 1) + \frac{n_2kT}{q}  \cdot I \cdot \ln(\frac{f\cdot I_p - I}{I_{0,2}} + 1)$

Find $f$ for maximum power:

&ensp; $\frac{\partial P}{\partial f} = 0$

&ensp; $0 = -\frac{n_1kT}{q}  \cdot I \cdot (\frac{I_p}{(1-f)\cdot I_p - I+I_{0,1}})
+\frac{n_2kT}{q}  \cdot I \cdot (\frac{I_p}{f\cdot I_p - I+I_{0,2}})$

Find $1-f$ and $f$:

&ensp; $1-f=(\frac{n_1}{n_1+n_2}) - (\frac{n_2-n_1}{n_1+n_2})\cdot\frac{I}{I_p}
-(\frac{n_2}{n_1+n_2})\cdot \frac{I_{0,1}}{I_p}
+(\frac{n_1}{n_1+n_2})\cdot \frac{I_{0,2}}{I_p}$

&ensp; $f=(\frac{n_2}{n_1+n_2}) + (\frac{n_2-n_1}{n_1+n_2})\cdot\frac{I}{I_p}
+(\frac{n_2}{n_1+n_2})\cdot \frac{I_{0,1}}{I_p}
-(\frac{n_1}{n_1+n_2})\cdot \frac{I_{0,2}}{I_p}$

The photocurrents are

&ensp; $I_{p,1}=(\frac{n_2}{n_1+n_2}) \cdot I_p + (\frac{n_1-n_2}{n_1+n_2})\cdot I +
(\frac{n_2}{n_1+n_2})\cdot I_{0,1}
-(\frac{n_1}{n_1+n_2})\cdot I_{0,2}$

&ensp; $I_{p,2}=(\frac{n_1}{n_1+n_2}) \cdot I_p - (\frac{n_1-n_2}{n_1+n_2})\cdot I -
(\frac{n_2}{n_1+n_2})\cdot I_{0,1}
+(\frac{n_1}{n_1+n_2})\cdot I_{0,2}$

