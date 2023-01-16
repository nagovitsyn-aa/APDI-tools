(* ::Package:: *)

BeginPackage["AlphaKappa`",{"BasicFunc`"}]

\[Alpha]\[Kappa];


Begin["`Private`"];

\[Alpha]\[Kappa][xdecay_, \[Omega]0_, \[Omega]uh_]= 
With[{\[Omega]s=\[Omega]0-\[Omega]uh,\[Omega]dce=\[Omega]ce[xdecay], \[Omega]dpe=\[Omega]pe[xdecay],k0=kX[\[Omega]0,xdecay],quh=qminus[\[Omega]uh,xdecay]},
With[{Nm=(c quh)/\[Omega]uh,Ns=(c kX[\[Omega]s,xdecay])/\[Omega]s,N0=(c k0)/\[Omega]0,\[Gamma]sx=1-\[Omega]dpe^2/(\[Omega]s^2-\[Omega]dce^2-\[Omega]dpe^2),\[Gamma]sy=1-\[Omega]dce^2/\[Omega]s^2 \[Omega]dpe^2/(\[Omega]s^2-\[Omega]dce^2-\[Omega]dpe^2)},
With[{\[CapitalOmega]0=\[Omega]s(\[Omega]0^2/\[Omega]dce^2-1)\[Gamma]sy+\[Gamma]sx(2\[Omega]0-\[Omega]s), \[CapitalOmega]s=\[Omega]0(\[Omega]s^2/\[Omega]dce^2-1)+\[Gamma]sx(2\[Omega]s-\[Omega]0)},
	{\[Alpha],\[Kappa]}=-(quh/c^2) el/(me c) (\[Omega]dce^2 \[Omega]dpe^2)/(\[Omega]uh^2-\[Omega]dce^2)\[Omega]s/(\[Omega]s^2-\[Omega]dce^2)1/(\[Omega]0^2-\[Omega]dce^2){\[Omega]s \[Omega]uh(Nm(\[Gamma]sx+(\[Omega]s \[Omega]0)/\[Omega]dce^2-\[Omega]0/\[Omega]s)+N0(2+\[Omega]0/\[Omega]s)\[Gamma]sx-Ns \[Gamma]sx),\[Omega]0(\[CapitalOmega]0 N0+\[CapitalOmega]s Ns)}
]]]

End[];
EndPackage[];


