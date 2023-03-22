(* ::Package:: *)

(* ::Text:: *)
(*Nm=(c qm)/\[Omega]m;Ns=(c ks)/\[Omega]s;N0=(c k0)/\[Omega]0;*)
(*\[Gamma]sx=1-\[Omega]pe^2/(\[Omega]s^2-\[Omega]ce^2-\[Omega]pe^2);*)
(*\[Gamma]sy=1-\[Omega]ce^2/\[Omega]s^2 \[Omega]pe^2/(\[Omega]s^2-\[Omega]ce^2-\[Omega]pe^2);*)
(*\[CapitalOmega]0=\[Omega]s(\[Omega]0^2/\[Omega]ce^2-1)\[Gamma]sy+\[Gamma]sx(2\[Omega]0-\[Omega]s);*)
(*\[CapitalOmega]s=\[Omega]0(\[Omega]s^2/\[Omega]ce^2-1)+\[Gamma]sx(2\[Omega]s-\[Omega]0);*)
(*\[Kappa]=-(qm/c^2) e/(m c) (\[Omega]ce^3 \[Omega]pe^2)/(\[Omega]m^2-\[Omega]ce^2)\[Omega]s/(\[Omega]s^2-\[Omega]ce^2)\[Omega]0/(\[Omega]0^2-\[Omega]ce^2)(\[CapitalOmega]0 N0+\[CapitalOmega]s Ns)*)
(*\[Alpha]=-(qm/c^2) e/(m c) (\[Omega]ce^2 \[Omega]pe^2)/(\[Omega]m^2-\[Omega]ce^2)\[Omega]s^2/(\[Omega]s^2-\[Omega]ce^2)\[Omega]m/(\[Omega]0^2-\[Omega]ce^2)(Nm(\[Gamma]sx+(\[Omega]s \[Omega]0)/\[Omega]ce^2-\[Omega]0/\[Omega]s)+N0(2+\[Omega]0/\[Omega]s)\[Gamma]sx-Ns \[Gamma]sx)*)


BeginPackage["AlphaKappa`",{"BasicFunc`"}]
\[Alpha]\[Kappa];
\[Alpha]\[Kappa]2;


Begin["`Private`"];

\[Alpha]\[Kappa][xdecay_, \[Omega]0_, \[Omega]uh_]= 
With[{\[Omega]s=\[Omega]0-\[Omega]uh,\[Omega]dce=\[Omega]ce[xdecay], \[Omega]dpe=\[Omega]pe[xdecay],k0=kX[\[Omega]0,xdecay],quh=qminus[\[Omega]uh,xdecay]},
With[{Nm=(c quh)/\[Omega]uh,Ns=(c kX[\[Omega]s,xdecay])/\[Omega]s,N0=(c k0)/\[Omega]0,\[Gamma]sx=1-\[Omega]dpe^2/(\[Omega]s^2-\[Omega]dce^2-\[Omega]dpe^2),\[Gamma]sy=1-\[Omega]dce^2/\[Omega]s^2 \[Omega]dpe^2/(\[Omega]s^2-\[Omega]dce^2-\[Omega]dpe^2)},
With[{\[CapitalOmega]0=\[Omega]s(\[Omega]0^2/\[Omega]dce^2-1)\[Gamma]sy+\[Gamma]sx(2\[Omega]0-\[Omega]s), \[CapitalOmega]s=\[Omega]0(\[Omega]s^2/\[Omega]dce^2-1)+\[Gamma]sx(2\[Omega]s-\[Omega]0)},
	{\[Alpha],\[Kappa]}=-(quh/c^2) el/(me c) (\[Omega]dce^2 \[Omega]dpe^2)/(\[Omega]uh^2-\[Omega]dce^2)\[Omega]s/(\[Omega]s^2-\[Omega]dce^2)1/(\[Omega]0^2-\[Omega]dce^2){\[Omega]s \[Omega]uh(Nm(\[Gamma]sx+(\[Omega]s \[Omega]0)/\[Omega]dce^2-\[Omega]0/\[Omega]s)+N0(2+\[Omega]0/\[Omega]s)\[Gamma]sx-Ns \[Gamma]sx),\[Omega]0(\[CapitalOmega]0 N0+\[CapitalOmega]s Ns)}
]]];


\[Alpha]\[Kappa]2[\[Omega]m_]= 
With[{\[Omega]s=\[Omega]0-\[Omega]mn,\[Omega]dce=\[Omega]ce[xdecay], \[Omega]dpe=\[Omega]pe[xdecay],k0=kX[\[Omega]0,xdecay],quh=qminus[\[Omega]uh,xdecay]},
With[{Nm=(c quh)/\[Omega]uh,Ns=(c kX[\[Omega]s,xdecay])/\[Omega]s,N0=(c k0)/\[Omega]0,\[Gamma]sx=1-\[Omega]dpe^2/(\[Omega]s^2-\[Omega]dce^2-\[Omega]dpe^2),\[Gamma]sy=1-\[Omega]dce^2/\[Omega]s^2 \[Omega]dpe^2/(\[Omega]s^2-\[Omega]dce^2-\[Omega]dpe^2)},
With[{\[CapitalOmega]0=\[Omega]s(\[Omega]0^2/\[Omega]dce^2-1)\[Gamma]sy+\[Gamma]sx(2\[Omega]0-\[Omega]s), \[CapitalOmega]s=\[Omega]0(\[Omega]s^2/\[Omega]dce^2-1)+\[Gamma]sx(2\[Omega]s-\[Omega]0)},
	{\[Alpha],\[Kappa]}=-(quh/c^2) el/(me c) (\[Omega]dce^2 \[Omega]dpe^2)/(\[Omega]uh^2-\[Omega]dce^2)\[Omega]s/(\[Omega]s^2-\[Omega]dce^2)1/(\[Omega]0^2-\[Omega]dce^2){\[Omega]s \[Omega]uh(Nm(\[Gamma]sx+(\[Omega]s \[Omega]0)/\[Omega]dce^2-\[Omega]0/\[Omega]s)+N0(2+\[Omega]0/\[Omega]s)\[Gamma]sx-Ns \[Gamma]sx),\[Omega]0(\[CapitalOmega]0 N0+\[CapitalOmega]s Ns)}
]]];

End[];
EndPackage[];
