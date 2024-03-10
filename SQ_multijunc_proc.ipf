#pragma rtGlobals=1		// Use modern global access method.

Constant hc=1239.842 //eV.nm
Constant c=2.997925E+8 //m/s
Constant kB=8.617342E-05 //eV/K
Constant qc=1.602176E-19 //J/eV

Function F_theta(theta)
Variable theta

	return pi*sin(theta/180*pi)^2//
End

// max concentraion
Function Xmax()

	NVAR vTheta_s=root:theta_s //sun semi-angle
	variable fs=sin(vTheta_s)^2

	return 1/fs
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
		I2_n=(x^2/n+2*x/n^2+2/n^3)*exp(-n*x)
		I2-=I2_n	
		n+=1
	while (I2_n > 1e-18)	
	return I2 + 2*zeta(3,1)
End

Function find_I3_E(E,T)
Variable E,T

	variable I3=0
	variable I3_n = 0
	variable x=E/kB/T
	variable n = 1
	do
		I3_n=(x^3/n+3*x^2/n^2+6*x/n^3+6/n^4)*exp(-n*x)
		I3-=I3_n	
		n+=1
	while (I3_n > 1e-18)	
	return I3+6*zeta(4,1)
End

Function find_I2_dE(E1,E2,T)
Variable E1,E2,T

	variable I2=0
	variable I2_n = 0
	variable x1=E1/kB/T
	variable x2=E2/kB/T
	variable n = 1
	do
		I2_n=(x2^2/n+2*x2/n^2+2/n^3)*exp(-n*x2)-(x1^2/n+2*x1/n^2+2/n^3)*exp(-n*x1)
		I2-=I2_n	
		n+=1
	while (I2_n > 1e-18)	
	return I2
End

Function find_I3_dE(E1,E2,T)
Variable E1,E2,T

	variable I3=0
	variable I3_n = 0
	variable x1=E1/kB/T
	variable x2=E2/kB/T
	variable n = 1
	do
		I3_n=(x2^3/n+3*x2^2/n^2+6*x2/n^3+6/n^4)*exp(-n*x2)-(x1^3/n+3*x1^2/n^2+6*x1/n^3+6/n^4)*exp(-n*x1)
		I3-=I3_n	
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

// photon flux within [0,E]
Function nr_E(E,T) // #/s.m^2
Variable E,T

	return c/4*n_E(E,T)	
End

// spectral irradiance within [0,E]
Function Lr_E(E,T) // W/m^2
Variable E,T

	return c/4*u_E(E,T)	
End

// photon density within [E1,E2]
Function n_dE(E1,E2,T) // #/m^3
Variable E1,E2,T

	return (E2>E1)?8*pi*(kB*T/hc)^3*find_I2_dE(E1,E2,T)*1e27	:0
End

// energy density within [E2,E2]
Function u_dE(E1,E2,T) // J/m^3
Variable E1,E2,T

	return (E2>E1)?8*pi*kB*T*(kB*T/hc)^3*find_I3_dE(E1,E2,T)*qc*1e27	 :0
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

// total power density (W/m^2)
Function Ptot(Xc,Ts)
Variable Xc,Ts

	NVAR vTheta_s=root:theta_s //sun semi-angle
	variable fs=sin(vTheta_s)^2

	return Xc*fs*Lr_E(30,Ts)
End

//--------------single-junction-----------------
// single-junction J(V) absorbing within [E1,E2]
Function J_1junc(V,E1,E2,Xc,Ts,Ta,n)
Variable V,E1,E2,Xc,Ts,Ta,n

	NVAR vTheta_s=root:theta_s //sun semi-angle
	variable fs=sin(vTheta_s)^2
	variable fc=1

	variable Jph=qc*Xc*fs*nr_dE(E1,E2,Ts)
	variable J0=qc*fc*nr_dE(E1,100,Ta)
	
	return Jph-J0*(exp(V/(n*kB*Ta))-1)
End

// single-junction P(V) absorbing within [E1,E2]
Function P_1junc(V,E1,E2,Xc,Ts,Ta,n)
Variable V,E1,E2
Variable Xc,Ts,Ta,n
	
	if(E2<=E1)
		return 0
	endif
	variable P=V*J_1junc(V,E1,E2,Xc,Ts,Ta,n)	
	return P
End

// optimization: single-junction power
Function wP_1junc(wPars,V)
	Wave wPars
	Variable V

	variable Xc=wPars[0]
	variable Ts=wPars[1]
	variable Ta=wPars[2]	
	variable n=wPars[3]
	variable E1=wPars[4]
	variable E2=wPars[5]

	variable P=P_1junc(V,E1,E2,Xc,Ts,Ta,n)
	return P
End

// max-power single-junction V absorbing within [Eg,Ec]
function findV_maxP_1junc(Eg,Ec,Xc,Ts,Ta,n)
Variable Eg,Ec,Xc,Ts,Ta,n

	string sP="root:"

	Make/D/O/N=6 $sP+"tempPars"/wave=wTempPars
	wTempPars[0]=Xc
	wTempPars[1]=Ts
	wTempPars[2]=Ta
	wTempPars[3]=n //n1
	wTempPars[4]=Eg
	wTempPars[5]=Ec

	Make/D/O/N=1 $sP+"tempData"/wave=wTempData
	wTempData[0]=Eg/2

	Optimize/A/X=wTempData/S=1 wP_1junc,wTempPars
	KillWaves $sP+"tempData"

	variable Vm=V_maxloc
	return Vm
End

function doFindMaxP_1junc()
	string sP="root:"
	
	wave w_1juncPars=$sP+"pars_1junc"

	variable Xc=w_1juncPars[0]
	variable Ts=w_1juncPars[1]
	variable Ta=w_1juncPars[2]	
	variable n=w_1juncPars[3]
	variable Eg=w_1juncPars[4]
	variable Ec=w_1juncPars[5]
		
	variable Vm=findV_maxP_1junc(Eg,Ec,Xc,Ts,Ta,n)	
	variable Jm=J_1junc(Vm,Eg,Ec,Xc,Ts,Ta,n)
	variable Pm=Vm*Jm	
	
	w_1juncPars[6]=Vm
	w_1juncPars[7]=Jm
	w_1juncPars[8]=Pm
End

Macro max_power_single_junction()
	doFindMaxP_1junc()
End

//--------------find J(V): two-junction, two-contact-----------------
// optimization: two-junction, two-contact current
Function wdJ_2junc_2cont(wPars,dV)
	Wave wPars
	Variable dV

	variable Xc=wPars[0]
	variable Ts=wPars[1]
	variable Ta=wPars[2]	
	variable n1=wPars[3]
	variable n2=wPars[4]
	variable Eg1=wPars[5]
	variable Eg2=wPars[6]
	variable Ec=wPars[7]
	variable V=wPars[8]

	variable V1=(V-dV)/2,V2=(V+dV)/2
	variable J1=J_1junc(V1,Eg1,Eg2,Xc,Ts,Ta,n1)
	variable J2=J_1junc(V2,Eg2,Ec,Xc,Ts,Ta,n2)
	return J1-J2
End

//two-junction, two-contact current absorbing within [E1,E2] and [E2,inf]
function find_dV_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2,V)
Variable Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2,V

	string sP="root:"
	Make/D/O/N=9 $sP+"tempPars"/wave=wTempPars
	wTempPars[0]=Xc
	wTempPars[1]=Ts
	wTempPars[2]=Ta
	wTempPars[3]=n1 //n1
	wTempPars[4]=n2 //n2
	wTempPars[5]=Eg1
	wTempPars[6]=Eg2
	wTempPars[7]=Ec
	wTempPars[8]=V

	Make/D/O/N=1 $sP+"tempData"/wave=wTempData
	wTempData[0]=0

	FindRoots/X=wTempData wdJ_2junc_2cont,wTempPars
	KillWaves $sP+"tempData"

	variable dV=V_root
	return dV
End

function doFind_dV_2junc_2cont()

	string sP="root:"	
	wave wPars=$sP+"pars_2junc_2cont"

	variable Xc=wPars[0]
	variable Ts=wPars[1]
	variable Ta=wPars[2]	
	variable n1=wPars[3]
	variable n2=wPars[4]
	variable Eg1=wPars[5]
	variable Eg2=wPars[6]
	variable Ec=wPars[7]
	variable V=wPars[8]

	variable dV=find_dV_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2,V)
	variable V1=(V-dV)/2,V2=(V+dV)/2	
	variable J=J_1junc(V1,Eg1,Eg2,Xc,Ts,Ta,n1)	
	
	wPars[9]=J
	wPars[10]=V*J
	wPars[11]=V1
	wPars[12]=V1*J
	wPars[13]=V2
	wPars[14]=V2*J
End

Macro current_2junction_2contact()
	doFind_dV_2junc_2cont()
End

//--------------find max power point: two-junction, two-contact-----------------
// optimization: two-junction, two-contact current
Function wP_2junc_2cont(wPars,V)
	Wave wPars
	Variable V

	variable Xc=wPars[0]
	variable Ts=wPars[1]
	variable Ta=wPars[2]	
	variable n1=wPars[3]
	variable n2=wPars[4]
	variable Eg1=wPars[5]
	variable Eg2=wPars[6]
	variable Ec=wPars[7]

	variable dV=find_dV_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2,V)
	variable V1=(V-dV)/2,V2=(V+dV)/2
	variable J=J_1junc(V1,Eg1,Eg2,Xc,Ts,Ta,n1)	
	variable P=V*J	
	return P
End

//two-junction, two-contact current absorbing within [E1,E2] and [E2,inf]
function find_Vmax_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2)
Variable Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2

	string sP="root:"
	Make/D/O/N=8 $sP+"tempPars"/wave=wTempPars
	wTempPars[0]=Xc
	wTempPars[1]=Ts
	wTempPars[2]=Ta
	wTempPars[3]=n1 //n1
	wTempPars[4]=n2 //n2
	wTempPars[5]=Eg1
	wTempPars[6]=Eg2
	wTempPars[7]=Ec

	Make/D/O/N=1 $sP+"tempData2"/wave=wTempData
	wTempData[0]=(Eg1+Eg2)/2

	Optimize/A/X=wTempData/S=1 wP_2junc_2cont,wTempPars
	KillWaves $sP+"tempData2"

	variable Vm=V_maxloc
	return Vm
End

function doFind_Vmax_2junc_2cont()

	string sP="root:"	
	wave wPars=$sP+"pars_2junc_2cont"

	variable Xc=wPars[0]
	variable Ts=wPars[1]
	variable Ta=wPars[2]	
	variable n1=wPars[3]
	variable n2=wPars[4]
	variable Eg1=wPars[5]
	variable Eg2=wPars[6]
	variable Ec=wPars[7]

	variable V=find_Vmax_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2)
	variable dV=find_dV_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2,V)
	variable V1=(V-dV)/2,V2=(V+dV)/2	
	variable J=J_1junc(V1,Eg1,Eg2,Xc,Ts,Ta,n1)	
	
	wPars[8]=V
	wPars[9]=J
	wPars[10]=V*J
	wPars[11]=V1
	wPars[12]=V1*J
	wPars[13]=V2
	wPars[14]=V2*J
End

Macro max_power_2junction_2contact()
	doFind_Vmax_2junc_2cont()
End

//--------------find optimal Eg1,Eg2: two-junction, two-contact-----------------
// optimization: two-junction, two-contact Eg1, Eg2
Function wPmax_Eg_2junc_2cont(wPars,Eg1,Eg2)
	Wave wPars
	Variable Eg1,Eg2

	variable Xc=wPars[0]
	variable Ts=wPars[1]
	variable Ta=wPars[2]	
	variable n1=wPars[3]
	variable n2=wPars[4]
	variable Ec=wPars[7]

	variable V=find_Vmax_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2)
	variable dV=find_dV_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2,V)
	variable V1=(V-dV)/2,V2=(V+dV)/2
	variable J=J_1junc(V1,Eg1,Eg2,Xc,Ts,Ta,n1)	
	variable P=V*J	
	return P
End

//two-junction, two-contact current absorbing within [E1,E2] and [E2,inf]
function find_Pmax_Eg_2junc_2cont(Ec,Xc,Ts,Ta,n1,n2)
Variable Ec,Xc,Ts,Ta,n1,n2

	string sP="root:"
	Make/D/O/N=8 $sP+"tempPars"/wave=wTempPars
	wTempPars[0]=Xc
	wTempPars[1]=Ts
	wTempPars[2]=Ta
	wTempPars[3]=n1 //n1
	wTempPars[4]=n2 //n2
	wTempPars[7]=Ec

	Make/D/O/N=2 $sP+"tempData3"/wave=wTempData
	wTempData[0]=0.9
	wTempData[1]=1.7

	Optimize/A/X=wTempData/S=1 wPmax_Eg_2junc_2cont,wTempPars
	variable Eg1=wTempData[0]
	variable Eg2=wTempData[1]
	KillWaves wTempData

	Make/D/O/N=2 $sP+"result2"/wave=wResult
	
	wResult[0]=Eg1
	wResult[1]=Eg2
End

function doFind_Pmax_Eg_2junc_2cont()

	string sP="root:"	
	wave wPars=$sP+"pars_2junc_2cont"

	variable Xc=wPars[0]
	variable Ts=wPars[1]
	variable Ta=wPars[2]	
	variable n1=wPars[3]
	variable n2=wPars[4]
	variable Ec=wPars[7]

	find_Pmax_Eg_2junc_2cont(Ec,Xc,Ts,Ta,n1,n2)
	wave wResult=$sP+"result2"
	variable Eg1=wResult[0]
	variable Eg2=wResult[1]

	variable V=find_Vmax_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2)
	variable dV=find_dV_2junc_2cont(Eg1,Eg2,Ec,Xc,Ts,Ta,n1,n2,V)
	variable V1=(V-dV)/2,V2=(V+dV)/2	
	variable J=J_1junc(V1,Eg1,Eg2,Xc,Ts,Ta,n1)	
	
	wPars[5]=Eg1
	wPars[6]=Eg2
	wPars[8]=V
	wPars[9]=J
	wPars[10]=V*J
	wPars[11]=V1
	wPars[12]=V1*J
	wPars[13]=V2
	wPars[14]=V2*J
End

Macro opt_bandgaps_2junction_2contact()
	doFind_Pmax_Eg_2junc_2cont()
End