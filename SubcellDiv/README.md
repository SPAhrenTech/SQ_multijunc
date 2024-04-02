
The optimization for maximum power of the partionining among constituent subcells of the total photocurrent available to a series-connected tandem can be performed simultaneously with the determination of the operating point for maximum power. In fact, the optimal photocurrent partitioning can be evaluated for any specified operating current, allowing determination of the maximum power point by variation of the device current as a single parameter. Some insight is gained by displaying the I-V curves for the individual subcells with that of the tandem device.
In the optimal configuration, more photocurrent should be allocated to subcells with larger ideality factor, due to the larger increase in voltage across those subcells. Differences in subcell dark currents have little effect on the optimal partitioning of photocurrent for the tandem device.

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

