(* ::Package:: *)

(* ============================================================================
   RGImprovedRateYukawa.wl -- RG-improved decay rate with a Yukawa fermion
   ----------------------------------------------------------------------------
   The fermion companion to RGImprovedRate.wl.  A single real scalar in the
   thin-wall model

       V = (lambda/8) (phi^2 - v^2)^2 + lambda Delta v^3 (phi - v)

   is coupled to a Dirac fermion by a Yukawa, m_psi = y phi.  The one-loop rate
   now has TWO fluctuation sectors -- the scalar determinant (+) and the fermion
   determinant (-1)^F (-) -- and TWO running couplings, lambda and y.  Two
   qualitatively new features relative to the pure scalar:

     * beta_lambda acquires the destabilising Yukawa term  -16 y^4  (the
       thin-wall analogue of the top-driven Standard-Model instability);
     * the field anomalous dimension is nonzero already at one loop,
       gamma_phi = 2 h y^2 (h = 1/16 pi^2) -- the SECOND running invariant that
       was absent for the pure scalar.

   The two-loop beta-functions come from RGBeta (Thomsen, arXiv:2101.08265) for
   the real-scalar + Dirac-fermion Yukawa theory; in the V superset lambda phi^4/8
   normalisation (RGBeta quartic lam = lambda/8, so lambda = 8 lam) and with
   y the physical Yukawa (m_psi = y phi, matching gamma_phi = 2 h y^2),

     beta_lambda = h (9 l^2 + 8 l y^2 - 16 y^4)
                 + h^2 (-51 l^3 - 36 l^2 y^2 + 28 l y^4 + 128 y^6),
     beta_y      = h (5 y^3) + h^2 (3/4 l^2 y - 6 l y^3 - 57/4 y^5),
     beta_m2/m2  = h (3 l + 4 y^2) + h^2 (-15/2 l^2 - 12 l y^2 - 2 y^4),
     gamma_phi   = h (2 y^2) + h^2 (3/4 l^2 - 5 y^4).

   Pipeline per (kappa0, mu), kappa = y/sqrt(lambda) the coupling ratio the
   fermion determinant depends on:
     1. run {lambda, y, m2} to mu (1- or 2-loop); mhat(mu) = sqrt(lambda) v,
        v^2 = -4 m2/lambda; kappa(mu) = y(mu)/sqrt(lambda(mu)).
     2. S0(mu) = S0hat(Delta)/lambda(mu) (the fermion has no vev; the bounce is
        the scalar bounce, driven by lambda alone).
     3. ln det_tot(mu) = Sigma^phi(mu) - Sigma^psi(mu):
          Sigma^phi(mu) = Sigma^phi(mhat) + I3^phi Log[mu/mhat]     (scalar, kappa-free),
          Sigma^psi(mu) = Sigma^psi_ren(kappa(mu)) + (beta^psi + I_g^psi) Log[mu/mhat]
                          (fermion: RECOMPUTED at kappa(mu), then MSbar-shifted;
                           beta^psi = I2^psi + I3^psi is the POTENTIAL running the
                           polygon captures; I_g^psi = -y^2 Int rho^3 phibar'^2 is the
                           GRADIENT moment it discards, supplied here from the bounce).
     4. ln(Gamma/V) = PolyDetRate[ln det_tot, S0] + 4 Log[mhat].

   THE NEW FEATURE (gamma_phi) AND ITS FIX (I_g^psi).  Unlike the pure scalar, the
   Yukawa activates the field anomalous dimension gamma_phi = 2 h y^2 at ONE loop.
   The polygonal fermion determinant is first-order (delta-free matching) and
   DISCARDS the within-shell gradient (d m_psi)^2 = y^2 (d phi)^2 that carries the
   wave-function running, so the NAIVE polygonal one-loop band is large (~3.6%), NOT
   the tiny Baratella value of the pure scalar (~0.01%): it is dominated by the
   uncaptured gradient.  The fix is the fermion analogue of the gauge gradient moment
   I_g (Rnu_gauge.tex): the missing determinant running is the CONTINUUM moment

       I_g^psi = -y^2 Int rho^3 phibar'^2 drho ,

   fixed by the same bridge as the gauge, gamma_phi Int (d phibar)^2 = -1/4 I_g^psi
   with gamma_phi = 2 y^2/16 pi^2.  Adding I_g^psi to the fermion determinant's
   running (NO ad-hoc action term -- exactly as the gauge puts I_g in the
   determinant) COLLAPSES the band from ~3.6% to ~0.7%.  The residual is the honest
   N^2LO + finite-N (the empirical minimum sits near -0.8 y^2 M_g rather than the
   analytic -1, the ~20% offset being the finite-N error of the polygonal potential
   running beta^psi = I_2^psi + I_3^psi).  This completes the Z_phi bookkeeping the
   companion note flagged; the fermion determinant now carries both its potential
   (beta^psi) and gradient (I_g^psi) running, like the gauge.

   Output: RGImprovedRateYukawa_band.{pdf,png}
   Run:    wolframscript -file RGImprovedRateYukawa.wl
           (needs PolyDet + FindBounce; RGBeta optional -- see RGImprovedRate.wl
            for the in-script extraction; here the validated values are used.)
   ============================================================================ *)

If[! MemberQ[$Packages, "PolyDet`"],
   Quiet@Check[Needs["PolyDet`"],
     PacletDirectoryLoad[FileNameJoin[Drop[FileNameSplit[$InputFileName], -3]]];
     Needs["PolyDet`"]]];
Needs["FindBounce`"];
(* the N=3 closed-form thin-wall fermion inputs live next to the figure scripts *)
fermInc = FileNameJoin[Append[Drop[FileNameSplit[$InputFileName], -4], "notebooks"], "ferm_tw_inputs.wls"];
If[FileExistsQ[fermInc], Get[fermInc],
   (* inline fallback: M_s = kappa phi_s on the exact case-b N=3 polygon *)
   fermTWInputs[Dl0_, kap_] := Module[{Dl = SetPrecision[Dl0, 50], allr, tT, tB, tF, Vf, a1, a2, R1, D2, D3, R0, R2},
     allr = x /. NSolve[1/2 x^3 - 1/2 x + Dl == 0, x, WorkingPrecision -> 50];
     {tT, tB, tF} = Sort[Re /@ allr];
     Vf[p_] := 1/8 (p^2 - 1)^2 + Dl (p - 1);
     a1 = (Vf[tB] - Vf[tT])/(tB - tT)/8; a2 = (Vf[tF] - Vf[tB])/(tF - tB)/8;
     R1 = (tF - tT)/2/(Sqrt[a1 (tB - tT)] - Sqrt[-a2 (tF - tB)]);
     D2 = Sqrt[(tB - tT)/a1]; D3 = Sqrt[(tF - tB)/(-a2)];
     R0 = Sqrt[R1 (R1 - D2)]; R2 = Sqrt[R1 (R1 + D3)];
     {kap {tT, (tT + tB)/2, (tB + tF)/2, tF}, {R0, R1, R2}}]];

h = 1/(16 Pi^2) // N;
(* full two-loop beta-functions (RGBeta-validated), pure functions of {l, y} *)
bLam[l_, y_, o_] := h (9 l^2 + 8 l y^2 - 16 y^4) +
   If[o == 2, h^2 (-51 l^3 - 36 l^2 y^2 + 28 l y^4 + 128 y^6), 0];
bYuk[l_, y_, o_] := h (5 y^3) + If[o == 2, h^2 (3/4 l^2 y - 6 l y^3 - 57/4 y^5), 0];
bMm[l_, y_, o_]  := h (3 l + 4 y^2) + If[o == 2, h^2 (-15/2 l^2 - 12 l y^2 - 2 y^4), 0];

(* couplings run to u = ln(mu/mhat0); mu = mhat0 = 1 at u = 0, lambda0 = v0 = 1,
   y0 = kappa0 (since sqrt(lambda0) = 1).  m2 = -lambda v^2/4 = -1/4 initially. *)
flow[kap0_, o_] := NDSolveValue[{
     la'[u] == bLam[la[u], yy[u], o], yy'[u] == bYuk[la[u], yy[u], o],
     mm'[u] == bMm[la[u], yy[u], o] mm[u],
     la[0] == 1, yy[0] == kap0, mm[0] == -1/4}, {la, yy, mm}, {u, -1, 1}];

(* ---- thin-wall scalar bounce + reference determinant ---------------------- *)
muB = Exp[EulerGamma/2]/(2 Sqrt[Pi]);
Vpp[p_] := (3 p^2 - 1)/2;  VV[p_, Dl_] := 1/8 (p^2 - 1)^2 + Dl (p - 1);
prepare[Dl_] := Module[{sp, tT, tF, bf, m2s, rr},
   sp = Sort[Re[x /. NSolve[x^3 - x + 2 Rationalize[Dl, 0] == 0, x]]];
   {tT, tF} = {sp[[1]], sp[[3]]};
   bf = FindBounce[VV[p, Dl], p, {N@tF, N@tT}, FieldPoints -> 31, Dimension -> 4];
   If[Head[bf] =!= BounceFunction, Return[$Failed, Module]];
   {m2s, rr} = PolyDetInputs[Vpp, bf];
   (* fermion gradient moment via the paclet: PolyDetIgPsi[bf, #&] = -Int rho^3 phibar'^2,
      so I_g^psi(y) = y^2 * IgUnit -- factor y^2 out to rescale to the running y(mu) *)
   <|"bf" -> bf, "S0hat" -> bf["Action"], "Dl" -> Dl, "IgUnit" -> PolyDetIgPsi[bf, # &],
     "lndetPhi0" -> PolyDetRen[m2s, rr, "Scale" -> muB, "Target" -> 10^-3],
     "I3phi" -> PolyDetIp[m2s, rr, 3][[3]]|>];

(* fermion renormalised determinant + running coefficient at coupling ratio kappa,
   on the SAME N=31 bounce (sectorMasses gives the linear mass M_s = kappa phi_s) *)
fermDet[bf_, kap_] := Module[{Mf, rf},
   {Mf, rf} = PolyDet`Private`sectorMasses[kap # &, bf];
   <|"sigPsiRen" -> PolyDetRenPsi[Mf, rf, "Scale" -> muB, "Target" -> 10^-3],
     "betaPsi" -> Total[PolyDetIpPsi[Mf, rf, 3][[2 ;; 3]]]|>];

(* ln(Gamma/V) at scale mu (ratio to mhat0).  igOn -> True adds the fermion
   determinant gradient moment I_g^psi = -y^2 Mg to the fermion running (the piece
   the first-order polygonal fermion chain discards) -- the gauge-analogous fix. *)
rate[dat_, laF_, yF_, mmF_, mu_, igOn_ : False] := Module[
   {u = Log[mu], laU, yU, m2U, mh, kapU, sPhi, fd, betaFull, sPsi, lndetTot},
   laU = laF[u]; yU = yF[u]; m2U = mmF[u];
   mh = Sqrt[laU] Sqrt[-4 m2U/laU];  kapU = yU/Sqrt[laU];
   sPhi = dat["lndetPhi0"] + dat["I3phi"] Log[mu/mh];
   fd = fermDet[dat["bf"], kapU];
   betaFull = fd["betaPsi"] + If[igOn, yU^2 dat["IgUnit"], 0];   (* potential beta^psi + gradient I_g^psi *)
   sPsi = fd["sigPsiRen"] + betaFull Log[mu/mh];
   lndetTot = sPhi - sPsi;                              (* fermion enters with -1 *)
   PolyDetRate[lndetTot, dat["S0hat"]/laU] + 4 Log[mh]];

(* ---- scan: 1-loop band without and WITH the Z_phi counterterm ------------- *)
kap0 = 6/10;                                            (* moderate Yukawa, lambda stays > 0 *)
grid = Range[0.02, 0.09, 0.007];
mus = {1/2, 1/Sqrt[2], 1, Sqrt[2], 2};
Print["Yukawa RG-improved rate, kappa0 = ", N@kap0, " ; scanning ", Length[grid], " Delta ..."];
scan = Table[Module[{dat = prepare[Dl], f1},
    If[dat === $Failed, Nothing,
      f1 = flow[kap0, 1];
      Module[{cen, rN, rZ},
        cen = rate[dat, f1[[1]], f1[[2]], f1[[3]], 1];
        rN = rate[dat, f1[[1]], f1[[2]], f1[[3]], #, False] & /@ mus;   (* naive: polygon only *)
        rG = rate[dat, f1[[1]], f1[[2]], f1[[3]], #, True] & /@ mus;    (* + fermion I_g^psi *)
        <|"Dl" -> Dl, "cen" -> cen, "bN" -> {Min[rN], Max[rN]}, "bG" -> {Min[rG], Max[rG]}|>]]],
   {Dl, grid}];
Print["Delta | ln(G/V) | naive band % | +I_g^psi band %"];
Scan[Print[#["Dl"], "  ", NumberForm[#["cen"], 8],
     "  ", NumberForm[100 (#["bN"][[2]] - #["bN"][[1]])/Abs[#["cen"]], 3],
     "  ", NumberForm[100 (#["bG"][[2]] - #["bG"][[1]])/Abs[#["cen"]], 3]] &, scan];

(* ---- figure: the naive band (gamma_phi uncaptured) vs the Z_phi-corrected band *)
colN = RGBColor[0.84, 0.37, 0.10]; colZ = RGBColor[0.13, 0.36, 0.62];
sty = Sequence[Frame -> True, FrameStyle -> Black, Axes -> False,
   GridLines -> Automatic, GridLinesStyle -> GrayLevel[0.85],
   LabelStyle -> Directive[FontFamily -> "Latin Modern Roman", 13, Black],
   ImageSize -> 500, PlotRangePadding -> Scaled[0.03]];
pct[e_, k_] := 100 (e[k] - e["cen"])/Abs[e["cen"]];
bNlo = {#["Dl"], pct[#, "bN"][[1]]} & /@ scan;  bNhi = {#["Dl"], pct[#, "bN"][[2]]} & /@ scan;
bGlo = {#["Dl"], pct[#, "bG"][[1]]} & /@ scan;  bGhi = {#["Dl"], pct[#, "bG"][[2]]} & /@ scan;
ymax = 1.15 Max[Abs[Join[bNlo, bNhi][[All, 2]]]];
figBand = Show[
   ListLinePlot[{bNlo, bNhi}, PlotStyle -> {{colN, AbsoluteThickness[1.3]}, {colN, AbsoluteThickness[1.3]}},
      Filling -> {1 -> {2}}, FillingStyle -> Directive[colN, Opacity[0.28]], InterpolationOrder -> 2],
   ListLinePlot[{bGlo, bGhi}, PlotStyle -> {{colZ, AbsoluteThickness[1.3]}, {colZ, AbsoluteThickness[1.3]}},
      Filling -> {1 -> {2}}, FillingStyle -> Directive[colZ, Opacity[0.5]], InterpolationOrder -> 2],
   Graphics[{GrayLevel[0.4], AbsoluteThickness[1], Dashing[{0.02, 0.015}],
      Line[{{First[grid], 0}, {Last[grid], 0}}]}],
   sty, PlotRange -> {{First[grid], Last[grid]}, {-ymax, ymax}},
   FrameLabel -> {"\[CapitalDelta]", "\[Delta] ln(\[CapitalGamma]/V) / |ln(\[CapitalGamma]/V)|   [%]"},
   PlotLabel -> None,
   Epilog -> {
     Text[Style["polygon only", 12, colN], Scaled[{0.70, 0.90}]],
     Text[Style["(\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Phi]\)]\) uncaptured)", 10, colN], Scaled[{0.70, 0.83}]],
     Text[Style["+ \!\(\*SubsuperscriptBox[\(I\), \(g\), \(\[Psi]\)]\) gradient moment", 12, colZ], Scaled[{0.66, 0.585}]],
     Text[Style["(\!\(\*SubsuperscriptBox[\(I\), \(g\), \(\[Psi]\)]\) = -\!\(\*SuperscriptBox[\(y\), \(2\)]\)\[Integral]\!\(\*SuperscriptBox[\(\[Rho]\), \(3\)]\)\!\(\*SuperscriptBox[OverscriptBox[\(\[Phi]\), \(.\)], \(2\)]\))", 10, colZ], Scaled[{0.66, 0.515}]],
     Text[Style["scalar + Yukawa fermion, \!\(\*SubscriptBox[\(\[Kappa]\), \(0\)]\) = 0.6", 10.5, GrayLevel[0.3]], Scaled[{0.36, 0.12}]]}];
outDir = If[$InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
Export[FileNameJoin[{outDir, "RGImprovedRateYukawa_band.pdf"}], figBand];
Export[FileNameJoin[{outDir, "RGImprovedRateYukawa_band.png"}], figBand, ImageResolution -> 150];
Print["wrote RGImprovedRateYukawa_band.{pdf,png} in ", outDir];
