#pragma rtGlobals=1		// Use modern global access method.

Constant hc=1239.842 //eV.nm
Constant c=2.997925E+8 //m/s
Constant kB=8.617342E-05 //eV/K
Constant qc=1.602176E-19 //J/eV

Function F_theta(theta)
Variable theta

	return pi*sin(theta/180*pi)^2//
End

Function flux_E(E,T,F,mu)
Variable E,T,F,mu

	return 2*F*c/hc^3*E^2/(exp((E-mu)/(kB*T))-1)*(1e9)^3//#/eV/m^2/s
End

Function flux_lambda(lambda,T,F,mu)
Variable lambda,T,F,mu

	return hc/lambda^2*flux_E(hc/lambda,T,F,mu)//#/nm/m^2/s
End

Function irr_E(E,T,F,mu)
Variable E,T,F,mu

	return E*flux_E(E,T,F,mu)*qc//W/eV/m^2
End

Function irr_lambda(lambda,T,F,mu)
Variable lambda,T,F,mu
	
	Variable E=hc/lambda
	return E*flux_lambda(lambda,T,F,mu)*qc//W/nm/m^2
End

Function j_E(E,T,F,mu)
Variable E,T,F,mu

	return qc*flux_E(E,T,F,mu)//A/eV/m^2
End

Function j_lambda(lambda,T,F,mu)
Variable lambda,T,F,mu

	return qc*flux_lambda(lambda,T,F,mu)//A/nm/m^2
End

Function dndE(E,T)
variable E,T
	
	return 8*pi*E^2/(hc)^3/(exp(E/(kB*T))-1)*(1e9^3)//#/eV/m^3
End

Function dudE(E,T)
variable E,T

	return qc*E*dndE(E,T)//J/eV/m^3
End

Function find_I2_inf()
variable T

	NVAR vI2_inf=root:I2_inf
	vI2_inf = 2*zeta(3,1)	
	return vI2_inf
End

Function find_I3_inf()
variable T

	NVAR vI3_inf=root:I3_inf
	vI3_inf = 6*zeta(4,1)
	return vI3_inf 
End

Function find_I2_E(E,T)
Variable E,T

	variable I2=0
	variable I2_n = 0
	variable x=E/kB/T
	variable n = 1
	do
		I2_n=2/n^3-(x^2/n+2*x/n^2+2/n^3)*exp(-n*x)
		I2+=I2_n	
		n+=1
	while (I2_n > 1e-18)	
	return I2
End

Function find_I3_E(E,T)
Variable E,T

	variable I3=0
	variable I3_n = 0
	variable x=E/kB/T
	variable n = 1
	do
		I3_n=6/n^4-(x^3/n+3*x^2/n^2+6*x/n^3+6/n^4)*exp(-n*x)
		I3+=I3_n	
		n+=1
	while (I3_n > 1e-18)	
	return I3
End

// photon density within [0,E]
Function n_E(E,T) // #/m^3
Variable E,T

	return 8*pi*(kB*T/hc)^3*find_I2_E(E,T)*1e27	
End

// energy density within [0,E]
Function u_E(E,T) // J/m^3
Variable E,T

	return 8*pi*kB*T*(kB*T/hc)^3*find_I3_E(E,T)*qc*1e27	
End

// photon density within [E1,E2]
Function n_dE(E1,E2,T) // #/m^3
Variable E1,E2,T

	return 8*pi*(kB*T/hc)^3*(find_I2_E(E2,T)-find_I2_E(E1,T))*1e27	
End

// energy density within [E2,E2]
Function u_dE(E1,E2,T) // J/m^3
Variable E1,E2,T

	return 8*pi*kB*T*(kB*T/hc)^3*(find_I3_E(E2,T)-find_I3_E(E1,T))*qc*1e27	
End

// photon flux within [E1,E2]
Function nr_dE(E1,E2,T) // #/s.m^2
Variable E1,E2,T

	return c/4*n_dE(E1,E2,T)	
End

// spectral irradiance within [E1,E2]
Function Lr_dE(E1,E2,T) // W/m^2
Variable E1,E2,T

	return c/4*u_dE(E1,E2,T)	
End

// single-junction J-V absorbing within [E1,E2]
Function J_V(V,E1,E2,Xc,Ts,Ta,n)
Variable V,E1,E2,Xc,Ts,Ta,n

	NVAR vTheta_s=root:theta_s //sun semi-angle
	variable fs=sin(vTheta_s)^2
	variable fc=1

	variable Jph=qc*Xc*fs*nr_dE(E1,E2,Ts)
	variable J0=qc*fc*nr_dE(E1,100,Ta)
	
	return Jph-J0*(exp(V/(n*kB*Ta))-1)
End

// max single-junction V absorbing within [E1,E2]
Function find_Vmax(E1,E2,Xc,Ts,Ta,n)
Variable E1,E2,Xc,Ts,Ta,n

	NVAR vTheta_s=root:theta_s //sun semi-angle
	variable fs=sin(vTheta_s)^2
	variable fc=1

	variable V=E1	
	variable dV=1e-4
	
	if(E2<=E1)
		return 0
	endif
	variable Pmax=V*J_V(V,E1,E2,Xc,Ts,Ta,n)
	
	variable iter=0
	variable Vplus,Vminus
	variable Pplus,Pminus
	do
		Vplus=V+dV
		Pplus=Vplus*J_V(Vplus,E1,E2,Xc,Ts,Ta,n)

		Vminus=V-dV
		Pminus=Vminus*J_V(Vminus,E1,E2,Xc,Ts,Ta,n)
	
		variable improved=0
		if(Pplus>Pmax)
			V=Vplus
			Pmax=Pplus
			improved=1
		endif
		if(Pminus>Pmax)
			V=Vminus
			Pmax=Pminus
			improved=1
		endif
		
		if(improved)
			dV=dV*1.5
		else
			dV=dV/2
		endif
		
		iter+=1
	while((dV>1e-18 )||(iter<100))
	return V
End

Function find_Pmax(E1,E2,Xc,Ts,Ta,n)
Variable E1,E2,Xc,Ts,Ta,n

	variable Vm=find_Vmax(E1,E2,Xc,Ts,Ta,n)
	Variable Jm=J_V(Vm,E1,E2,Xc,Ts,Ta,n)
	return Vm*Jm
End

// double-junction J-V absorbing absorbing with [E1,E2] and [E2,inf]
Function J_V_2junc(V,E1,E2,Xc,Ts,Ta,n1,n2)
Variable V,E1,E2,Xc,Ts,Ta,n1,n2
	
	variable J=0
	return J
	
	variable V1=E1
	variable dV1=0.01
	variable alpha=0.5
	variable J1,J2
	variable iter=0
	do
		J1=J_V(V1,E1,E2,Xc,Ts,Ta,n1)
		J2=J_V(V-V1,E2,100,Xc,Ts,Ta,n2)
				
		variable f0=J1-J2
		
		variable V1_plus=V1+dV1/2
		J1=J_V(V1_plus,E1,E2,Xc,Ts,Ta,n1)
		J2=J_V(V-V1_plus,E2,100,Xc,Ts,Ta,n2)
		variable f_plus=J1-J2
		
		variable V1_minus=V1-dV1/2
		J1=J_V(V1_minus,E1,E2,Xc,Ts,Ta,n1)
		J2=J_V(V-V1_minus,E2,100,Xc,Ts,Ta,n2)
		variable f_minus=J1-J2
		
		variable dfdV1=(f_plus-f_minus)/dV1
		variable V1inc=-alpha*f0/dfdV1
		
		if(abs(V1inc)>0.1)
			V1inc=sign(V1inc)*0.1
		endif
		V1+=V1inc
		dV1=V1inc/1000
		iter+=1
					
	while((abs(f0)>1e-18)||(iter<500))
			
	return J_V(V1,E1,E2,Xc,Ts,Ta,n1)
End

//
Function find_Vmax_2junc(E1,E2,Xc,Ts,Ta,n1,n2)
Variable E1,E2,Xc,Ts,Ta,n1,n2

	variable V=E1+(E2-E1)/2	
	variable dV=1e-3
	
	variable Pmax=V*J_V_2junc(V,E1,E2,Xc,Ts,Ta,n1,n2)
		
	variable iter=0
	do
		variable Vplus=V+dV
		variable Pplus=Vplus*J_V_2junc(Vplus,E1,E2,Xc,Ts,Ta,n1,n2)

		variable Vminus=V-dV
		variable Pminus=Vminus*J_V_2junc(Vminus,E1,E2,Xc,Ts,Ta,n1,n2)
	
		variable improved=0
		if(Pplus>Pmax)
			V=Vplus
			Pmax=Pplus
			improved=1
		endif
		if(Pminus>Pmax)
			V=Vminus
			Pmax=Pminus
			improved=1
		endif
		
		if(improved)
			dV=dV*1.5
		else
			dV=dV/2
		endif
		
		iter+=1
	while((dV>1e-8 )&&(iter<400))
	return V
End

Function find_Pmax_2junc(E1,E2,Xc,Ts,Ta,n1,n2)
Variable E1,E2,Xc,Ts,Ta,n1,n2

	variable Vm=find_Vmax_2junc(E1,E2,Xc,Ts,Ta,n1,n2)
	Variable Jm=J_V_2junc(Vm,E1,E2,Xc,Ts,Ta,n1,n2)
	return Vm*Jm
End