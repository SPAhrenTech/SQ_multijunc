#pragma rtGlobals=1		// Use modern global access method.

Constant hc=1239.842//eV·nm
Constant c=2.997925E+8//m/s
Constant kB=8.617342E-05//eV/K
Constant qc=1.602176E-19//J/eV
	
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

Function n_E(E1,E2,T)
variable E1,E2,T

	variable dE=1e-2
	variable E=E1+dE/2
	variable tot=0
	do
		if(E>E2)
			break
		endif
		tot+=dndE(E,T)*dE
		E+=dE
	while(1)
	return tot	
End

Function u_E(E1,E2,T)
Variable E1,E2,T

	variable dE=1e-2
	variable E=E1+dE/2
	variable tot=0
	do
		if(E>E2)
			break
		endif
		tot+=dudE(E,T)*dE
		E+=dE
	while(1)
	return tot	
End

Function Nr_E(E1,E2,T)
Variable E1,E2,T

	return c/4*n_E(E1,E2,T)	
End

Function Lr_E(E1,E2,T)
Variable E1,E2,T

	return c/4*u_E(E1,E2,T)	
End

Function fnE_T(wnEdiff_T,E1,E2)
Wave wnEdiff_T

	variable E1,E2
	return wnEdiff_T(E1)(E2)	
End

Function fuE_T(wuEdiff_T,E1,E2)
Wave wuEdiff_T

	variable E1,E2
	return wuEdiff_T(E1)(E2)	
End

Function fNrE_T(wnEdiff_T,E1,E2)
wave wnEdiff_T

	variable E1,E2
	return c/4*fnE_T(wnEdiff_T,E1,E2)	
End

Function fLrE_Ts(wuEdiff_T,E1,E2)
Wave wuEdiff_T
Variable E1,E2

	return c/4*fuE_T(wuEdiff_T,E1,E2)	
End

Struct f_info
	variable Ftot
	variable nSteps
End

Function find_nE_inf_T(T)
variable T

	variable Ftot = 0.
	variable Eprev = 0
	variable fprev = 0
	variable nSteps=0
	variable dEref = kB*T/1e6
	
	variable E = dEref
	variable f=dndE(E,T)
	variable fAvg=(f+fprev)/2
	
	variable dfdEref=(f-fprev)/(E-Eprev)
	variable dE=dEref
	Ftot+=fAvg*dEref
	fprev=f
	Eprev=E
	E+=dE
	do	
		f=dndE(E,T)
		fAvg=(f+fprev)/2
		Ftot+=fAvg*dE

		variable dfdE=(f-fprev)/dE
		variable dfdEx=sqrt(dfdE^2+(1.e-9)^2)	
		dE=(dfdEref/dfdEx)*dEref

		fprev=f
		Eprev=E
		E+=dE

		nSteps+=1
	while((dfdE>0)||(abs(dfdE)>dfdEref/1e4))
	print nSteps		
	return Ftot
End

Function find_nE_T(n_T,E_T,T,Npoints)
wave n_T,E_T
variable T,Npoints

	Redimension /N=Npoints n_T,E_T
	variable nxp=100
	
	//integrate
	variable Ftot = 0.
	variable Eprev = 0
	variable fprev = 0
	variable nSteps=0
	variable dEref = kB*T/1e6
	
	variable E = dEref
	variable f=dndE(E,T)
	variable fAvg=(f+fprev)/2
	
	variable dfdEref=(f-fprev)/(E-Eprev)
	variable dE=dEref
	Ftot+=fAvg*dEref
	fprev=f
	Eprev=E
	E+=dE
	do	
		f=dndE(E,T)
		fAvg=(f+fprev)/2
		Ftot+=fAvg*dE

		variable dfdE=(f-fprev)/dE
		variable dfdEx=sqrt(dfdE^2+(1.e-9)^2)	
		dE=(dfdEref/dfdEx)*dEref

		fprev=f
		Eprev=E
		E+=dE

		nSteps+=1
	while((dfdE>0)||(abs(dfdE)>dfdEref/1e4))

	variable tot=0
	variable i,j
	for(i=0;i<nx;i+=1)	
		nE_T[i]=tot

		for(j=0;j<nxp;j+=1)		
			variable E=E0+i*dE+j*dE/nxp
			variable f_E=dndE(E,T)
			tot+=f_E*dE/nxp			
		endfor
	endfor
End
			
//
Function find_uE_T(uE_T,T)
wave uE_T
variable T

	variable nx=DimSize(uE_T,0)
	variable nxp=100
	variable dE=dimDelta(uE_T,0)
	variable E0=dimOffset(uE_T,0)
	
	variable tot=0
	variable i,j
	for(i=0;i<nx;i+=1)	
		uE_T[i]=tot

		for(j=0;j<nxp;j+=1)		
			variable E=E0+i*dE+j*dE/nxp
			variable f_E=dudE(E,T)
			tot+=f_E*dE/nxp			
		endfor
	endfor
End
			
//
Function find_Ptot()

	variable Ts=5960//K
	variable theta_s=atan(1.39E6/1.50E8/2)//sun semi-angle
	variable fs=sin(theta_s)^2

	return fs*Lr_E(0,10,Ts)
End

//
Function J_V(V,Eg1,Eg2,Xc)
Variable V,Eg1,Eg2,Xc


	wave wnEdiff_Ts=root:nEdiff_Ts
	wave wuEdiff_Ts=root:uEdiff_Ts
	
	wave wnEdiff_Ta=root:nEdiff_Ta
	wave wuEdiff_Ta=root:uEdiff_Ta

	variable theta_s=atan(1.39E+06/1.50E+08/2)
	variable fs=sin(theta_s)^2
	variable fc=1
	variable n=1

	variable Jph=qc*Xc*fs*fNrE_T(wnEdiff_Ts,Eg1,Eg2)
	variable J0=qc*fc*fNrE_T(wnEdiff_Ta,Eg1,10)
	
	return Jph-J0*(exp(V/(n*kB*Ta))-1)
End

//
Function find_Vmax(Eg1,Eg2,Xc)
Variable Eg1,Eg2,Xc

	variable theta_s=atan(1.39E+06/1.50E+08/2)//sun semi-angle
	variable fs=sin(theta_s)^2
	variable fc=1

	wave wnEdiff_Ts=root:nEdiff_Ts
	wave wuEdiff_Ts=root:uEdiff_Ts
	
	wave wnEdiff_Ta=root:nEdiff_Ta
	wave wuEdiff_Ta=root:uEdiff_Ta

	variable J_ph=qc*Xc*fs*fNrE_T(wnEdiff_Ts,Eg1,Eg2)
	variable J_0=qc*fc*fNrE_T(wnEdiff_Ta,Eg1,10)
	
	variable V_oc=kB*Ta*ln(J_ph/J_0+1)
	
	variable dV=1e-4
	variable V=V_oc/2	
	
	if(Eg2<=Eg1)
		return 0
	endif
	variable Pmax=V*J_V(V,Eg1,Eg2,Xc)
	
	variable iter=0
	variable Vplus,Vminus
	variable Pplus,Pminus
	do
		Vplus=V+dV
		Pplus=Vplus*J_V(Vplus,Eg1,Eg2,Xc)

		Vminus=V-dV
		Pminus=Vminus*J_V(Vminus,Eg1,Eg2,Xc)
	
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
	while((dV>1e-6 )&&(iter<100))
	return V
End

Function J_V_2junc(V,Eg1,Eg2,Xc)
Variable V,Eg1,Eg2,Xc

	variable theta_s=atan(1.39E+06/1.50E+08/2)
	variable fs=sin(theta_s)^2
	variable fc=1
	variable n=1

	variable kT=kB*Ta
	
	wave wnEdiff_Ts=root:nEdiff_Ts
	wave wuEdiff_Ts=root:uEdiff_Ts
	
	wave wnEdiff_Ta=root:nEdiff_Ta
	wave wuEdiff_Ta=root:uEdiff_Ta

	variable Jph_a=qc*Xc*fs*fNrE_T(wnEdiff_Ts,Eg1,Eg2)
	variable J0_a=qc*fc*fNrE_T(wnEdiff_Ta,Eg1,10)

	variable Jph_b=qc*Xc*fs*fNrE_T(wnEdiff_Ts,Eg2,10)
	variable J0_b=qc*fc*fNrE_T(wnEdiff_Ta,Eg2,10)

	variable Jnet_a=Jph_a+J0_a
	variable Jnet_b=Jph_b+J0_b
	
	variable cA=1
	variable cB=-(Jnet_a+Jnet_b)
	variable cC=Jnet_a*Jnet_b-J0_a*J0_b*(exp(V/kT)-1)
	
	variable arg=cB^2-4*cA*cC
	if(arg<0)
		return -1
	endif
	
	variable Jplus=(-cB+sqrt(arg))/2/cA
	variable Jminus=(-cB-sqrt(arg))/2/cA
	
	variable J=Jminus
	return J
	
	variable Va=kT*ln((Jph_a-J)/J0_a+1)
	variable dV=0.01
	variable alpha=0.5
	variable Ja,Jb
	variable iter=0
	do
		Ja=J_V(Va,Eg1,Eg2,Xc)
		Jb=J_V(V-Va,Eg2,10,Xc)
				
		variable f0=Ja-Jb
		
		variable Va_plus=Va+dV/2
		Ja=J_V(Va_plus,Eg1,Eg2,Xc)
		Jb=J_V(V-Va_plus,Eg2,10,Xc)
		variable f_plus=Ja-Jb
		
		variable Va_minus=Va-dV/2
		Ja=J_V(Va_minus,Eg1,Eg2,Xc)
		Jb=J_V(V-Va_minus,Eg2,10,Xc)
		variable f_minus=Ja-Jb
		
		variable dfdV=(f_plus-f_minus)/dV
		variable Vinc=-alpha*f0/dfdV
		
		if(abs(Vinc)>0.1)
			Vinc=sign(Vinc)*0.1
		endif
		Va+=Vinc
		dV=Vinc/1000
		iter+=1
					
	while((abs(f0)>1e-6)&&(iter<500))
			
	return J_V(Va,Eg2,10,Xc)
End

//
Function find_Vmax_2junc(Eg1,Eg2,Xc)
Variable Eg1,Eg2,Xc

	variable theta_s=atan(1.39E+06/1.50E+08/2)//sun semi-angle
	variable fs=sin(theta_s)^2
	variable fc=1

	wave wnEdiff_Ts=root:nEdiff_Ts
	wave wuE_Ts=root:uEdiff_Ts
	
	wave wnEdiff_Ta=root:nEdiff_Ta
	wave wuEdiff_Ta=root:uEdiff_Ta
	
	//if(Eg1>=Eg2)
		//return 0
	//endif
	
	variable J1_ph=qc*Xc*fs*fNrE_T(wnEdiff_Ts,Eg2,10)
	variable J1_0=qc*fc*fNrE_T(wnEdiff_Ta,Eg2,10)
	
	variable J2_ph=qc*Xc*fs*fNrE_T(wnEdiff_Ts,Eg1,Eg2)
	variable J2_0=qc*fc*fNrE_T(wnEdiff_Ta,Eg1,10)

	variable V1_oc=kB*Ta*ln(J1_ph/J1_0+1)
	variable V2_oc=kB*Ta*ln(J2_ph/J2_0+1)	
	variable Voc=V1_oc+V1_oc
	
	variable V=Voc/2	
	variable dV=1e-3
	
	variable Pmax=V*J_V_2junc(V,Eg1,Eg2,Xc)
		
	variable iter=0
	do
		variable Vplus=V+dV
		variable Pplus=Vplus*J_V_2junc(Vplus,Eg1,Eg2,Xc)

		variable Vminus=V-dV
		variable Pminus=Vminus*J_V_2junc(Vminus,Eg1,Eg2,Xc)
	
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

Function find_Pmax(Eg1,Eg2,Xc)
Variable Eg1,Eg2,Xc

	variable Vm=find_Vmax(Eg1,Eg2,Xc)
	Variable Jm=J_V(Vm,Eg1,Eg2,Xc)
	return Vm*Jm
End
	
Function V_J(J,Eg1,Eg2,Xc)
Variable J,Eg1,Eg2,Xc

	variable theta_s=atan(1.39E+06/1.50E+08/2)
	variable fs=sin(theta_s)^2
	variable fc=1
	variable n=1

	variable kT=kB*Ta
	
	wave wnEdiff_Ts=root:nEdiff_Ts
	wave wuEdiff_Ts=root:uEdiff_Ts
	
	wave wnEdiff_Ta=root:nEdiff_Ta
	wave wuEdiff_Ta=root:uEdiff_Ta

	variable Jph=qc*Xc*fs*fNrE_T(wnEdiff_Ts,Eg1,Eg2)
	variable J0=qc*fc*fNrE_T(wnEdiff_Ta,Eg1,10)
	
	variable V=kT*ln((Jph-J)/J0+1)
	return V
End

Function find_Pmax_2junc(Eg1,Eg2,Xc)
Variable Eg1,Eg2,Xc

	variable Vm=find_Vmax_2junc(Eg1,Eg2,Xc)
	Variable Jm=J_V_2junc(Vm,Eg1,Eg2,Xc)
	return Vm*Jm
End
	
Function opt_Pmax_2junc(w,Eg1,Eg2)
Variable Eg1,Eg2
Wave w
	return find_Pmax_2junc(Eg1,Eg2,w[0])
End