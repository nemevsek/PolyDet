(* ::Package:: *)

(* ============================================================================
   RGImprovedRate.wl -- RG-improved false-vacuum decay rate with PolyDet + RGBeta
   ----------------------------------------------------------------------------
   A standalone PolyDet example: the one-loop decay rate of the single-scalar
   thin-wall model,

       V = (lambda/8) (phi^2 - v^2)^2 + lambda Delta v^3 (phi - v) ,

   RG-improved with the TWO-LOOP beta functions of the real-singlet quartic
   theory, obtained from the RGBeta package (Thomsen, arXiv:2101.08265) when
   available and from the (RGBeta-verified) hardcoded coefficients otherwise.

   The pipeline per (Delta, mu):
     1. couplings run to mu:  lambda(mu), m2(mu)  [m2 = -lambda v^2/4, the phi^2
        coefficient], with 1-loop or 2-loop truncation; v(mu)^2 = -4 m2/lambda;
        mhat(mu) = Sqrt[lambda] v (the tree FV mass).  Delta is held fixed: for
        this family t/(lambda v^3) is one-loop RG-invariant (beta_lambda/lambda
        cancels (3/2) dln v^2/dln mu exactly), and no cubic is generated (odd
        1PI functions vanish; the tilt only enters 1PR trees).
     2. S0(mu) = S0hat(Delta)/lambda(mu)   [exact family scaling in D = 4],
        S0hat from FindBounce at lambda = v = 1.
     3. ln det_ren(mu) = ln det_ren(mhat) + I_3 Log[mu/mhat(mu)]   [the exact
        MSbar running d Sigma^r/d ln mu = I_3 of the renormalised determinant],
        the reference determinant from the standard FindBounce -> PolyDetInputs
        -> PolyDetRen route (machine inputs suffice: the offset assembly).
     4. ln(Gamma/V / mhat0^4) = PolyDetRate[lndet, S0] + 4 Log[mhat(mu)/mhat0].

   POWER COUNTING.  Rescaling phi = v phi~ puts the coupling out front,
   S = (1/lambda) Shat[phi~], so 1/lambda plays the role of 1/hbar and the
   loop expansion is an expansion in lambda:
       ln(Gamma/V) = -Shat/lambda    (LO,   O(1/lambda): the bounce)
                     - (1/2) ln det' (NLO,  O(1):        one-loop determinant)
                     + O(lambda)      (N2LO, O(lambda):   two loops).

   WHY A SCALE BAND AT ONE LOOP.  The all-orders rate is mu-independent, but a
   TRUNCATED result is mu-independent only THROUGH the computed order: the
   leftover is the first neglected order.  The Baratella cancellation
   d/dln mu [S_R + Sigma^r/2] = 0 is an EXACTLY-one-loop identity (the explicit
   I_3 ln mu of the determinant cancels the one-loop running of lambda in S_R),
   so "exact at one loop" means the residual is N2LO.  With 1-loop running the
   band is tiny -- and, because mu = mhat sits near the stationary point of the
   rate (a principle-of-minimal-sensitivity that emerges here), quadratically
   small (~0.001-0.02%).  Turning on 2-loop running of lambda in S0 without a
   matching 2-loop determinant deliberately UN-balances the cancellation and
   exposes the N2LO-size sensitivity (~0.27%) -- the honest higher-order
   estimate.  (Not included, same N2LO order as the unknown two-loop finite
   determinant: gamma_phi field rescaling, the 2-loop drift of Delta.)

   Output: RGImprovedRate_rate.{pdf,png}  -- the NLO rate -ln(Gamma/V) Delta^3 vs Delta
           RGImprovedRate_band.{pdf,png}  -- fractional mu-band %, 1-loop vs 2-loop running
   Run:    wolframscript -file RGImprovedRate.wl     (from any directory;
           needs the PolyDet paclet and FindBounce installed; RGBeta optional)
   ============================================================================ *)

(* PolyDet: use the installed paclet if present, else the source tree
   (this file lives at <repo>/mathematica/PolyDet/Examples/) *)
If[! MemberQ[$Packages, "PolyDet`"],
   Quiet@Check[Needs["PolyDet`"],
     PacletDirectoryLoad[FileNameJoin[Drop[FileNameSplit[$InputFileName], -3]]];
     Needs["PolyDet`"]]];
Needs["FindBounce`"];

(* ---- two-loop singlet betas: from RGBeta when available ------------------- *)
(* Paper normalisation V superset lambda phi^4/8: the RGBeta quartic (term
   lamRGB phi^4) is lamRGB = lambda/8; loop factor h = 1/(16 Pi^2) explicit.
   RGBeta's init.m runs Check[Get, ...; Abort[]] and can Abort on a benign
   message even though the package loads (BetaTerm ends up defined); CheckAbort
   at TOP LEVEL survives it, and we detect success from $Packages, not the
   return value. *)
CheckAbort[Quiet@Needs["RGBeta`"], Null];
(* the RGBeta extraction is deterministic; it only emits a benign
   FrontEndObject::notavail during Simplify in a headless kernel, which would trip
   Check -- so Quiet it and VALIDATE the numbers instead of trapping messages *)
betaSource = If[MemberQ[$Packages, "RGBeta`"] && Length[DownValues[RGBeta`BetaTerm]] > 0,
   Module[{b1, b2, bm1, bm2},
     ResetModel[];
     AddScalar[phi, SelfConjugate -> True];
     AddQuartic[lam, {phi, phi, phi, phi}];
     AddScalarMass[m2c, {phi, phi}];
     {b1, b2, bm1, bm2} = Quiet[{BetaTerm[lam, 1], BetaTerm[lam, 2], BetaTerm[m2c, 1], BetaTerm[m2c, 2]}];
     (* convert lam -> lambda/8 and read off the numeric coefficients *)
     {c1, c2} = Quiet@Simplify[{8 (b1 /. lam -> \[Lambda]/8)/\[Lambda]^2, 8 (b2 /. lam -> \[Lambda]/8)/\[Lambda]^3}];
     {d1, d2} = Quiet@Simplify[{(bm1 /. {lam -> \[Lambda]/8, m2c -> 1})/\[Lambda], (bm2 /. {lam -> \[Lambda]/8, m2c -> 1})/\[Lambda]^2}];
     If[VectorQ[{c1, c2, d1, d2}, NumericQ] && {c1, c2, d1, d2} == {9, -51, 3, -15/2},
        "RGBeta " <> RGBeta`$RGBetaVersion, $Failed]],
   $Failed];
If[betaSource === $Failed,
   (* RGBeta-verified coefficients (real singlet, MSbar): in the 4!-normalised
      lambdabar = 3 lambda these are the textbook 3 lb^2 - (17/3) lb^3,
      lb - (5/6) lb^2, gamma = lb^2/12 *)
   {c1, c2, d1, d2} = {9, -51, 3, -15/2};
   betaSource = "hardcoded (RGBeta-verified)"];
Print["beta source: ", betaSource, "   beta_lambda = ", c1, " h l^2 + (", c2, ") h^2 l^3,",
      "   beta_m2/m2 = ", d1, " h l + (", d2, ") h^2 l^2,   h = 1/(16 Pi^2)"];

h = 1/(16 Pi^2) // N;
flow[order_] := Module[{sol},                                (* {lambda[u], m2[u]}, u = ln mu *)
   sol = NDSolveValue[{
      la'[u] == c1 h la[u]^2 + If[order == 2, c2 h^2 la[u]^3, 0],
      mm'[u] == d1 h la[u] mm[u] + If[order == 2, d2 h^2 la[u]^2 mm[u], 0],
      la[0] == 1, mm[0] == -1/4},                            (* lambda0 = v0 = 1: m2 = -lambda v^2/4 *)
      {la, mm}, {u, -1., 1.}];
   sol];
{lam1L, m21L} = flow[1];  {lam2L, m22L} = flow[2];

(* ---- the thin-wall model at lambda = v = 1 -------------------------------- *)
muB = Exp[EulerGamma/2]/(2 Sqrt[Pi]);                        (* "Scale" -> muB  <=>  mu = mhat *)
Vpp[p_] := (3 p^2 - 1)/2;
VV[p_, Dl_] := 1/8 (p^2 - 1)^2 + Dl (p - 1);
prepare[Dl_] := Module[{sp, tT, tF, bf, m2s, rr, lndet0, I3},
   sp = Sort[Re[x /. NSolve[x^3 - x + 2 Rationalize[Dl, 0] == 0, x]]];
   {tT, tF} = {sp[[1]], sp[[3]]};
   bf = FindBounce[VV[p, Dl], p, {N@tF, N@tT}, FieldPoints -> 31, Dimension -> 4];
   If[Head[bf] =!= BounceFunction, Return[$Failed, Module]];
   {m2s, rr} = PolyDetInputs[Vpp, bf];
   lndet0 = PolyDetRen[m2s, rr, "Scale" -> muB, "Target" -> 10^-3];
   I3 = PolyDetIp[m2s, rr, 3][[3]];
   <|"S0hat" -> bf["Action"], "lndet0" -> lndet0, "I3" -> I3|>];

(* ln(Gamma/V) in units of mhat0^4 = 1, at scale mu (in mhat0 units) *)
rate[dat_, lamF_, m2F_, mu_] := Module[{u = Log[mu], laU, m2U, vU, mhU, lndet, S0},
   laU = lamF[u];  m2U = m2F[u];
   vU = Sqrt[-4 m2U/laU];  mhU = Sqrt[laU] vU;
   lndet = dat["lndet0"] + dat["I3"] Log[mu/mhU];            (* d Sigma^r/d ln mu = I_3 *)
   S0 = dat["S0hat"]/laU;
   PolyDetRate[lndet, S0] + 4 Log[mhU]];

(* ---- scan ------------------------------------------------------------------ *)
grid = Range[0.02, 0.10, 0.008];
mus = {1/2, 1/Sqrt[2], 1, Sqrt[2], 2};          (* sample to catch the extremum near mu = mhat *)
Print["scanning ", Length[grid], " Delta points ..."];
scan = Table[Module[{dat = prepare[Dl]},
    If[dat === $Failed, Nothing,
      Module[{cen, r1, r2, chk},
        cen = rate[dat, lam1L, m21L, 1];        (* mu = mhat: 1L = 2L (common boundary) *)
        r1 = rate[dat, lam1L, m21L, #] & /@ mus; (* 1-loop-run rate over the mu window *)
        r2 = rate[dat, lam2L, m22L, #] & /@ mus; (* 2-loop-run rate over the mu window *)
        (* one-loop Baratella cancellation check: (1/2) I_3 vs S0 beta_lambda/lambda *)
        chk = dat["I3"]/2/(dat["S0hat"] c1 h);
        <|"Dl" -> Dl, "cen" -> cen,
          "b1" -> {Min[r1], Max[r1]}, "b2" -> {Min[r2], Max[r2]}, "chk" -> chk|>]]],
   {Dl, grid}];
Print["Delta | ln(G/V) | 1L band %% | 2L band %% | I3/2/(S0 b1) (=1 if Baratella exact)"];
Scan[Print[#["Dl"], "  ", NumberForm[#["cen"], 8],
     "  ", NumberForm[100 (#["b1"][[2]] - #["b1"][[1]])/Abs[#["cen"]], 3],
     "  ", NumberForm[100 (#["b2"][[2]] - #["b2"][[1]])/Abs[#["cen"]], 3],
     "  ", NumberForm[#["chk"], 5]] &, scan];

(* ---- figures ---------------------------------------------------------------- *)
colC = RGBColor[0.84, 0.37, 0.10]; col1 = RGBColor[0.55, 0.70, 0.42]; col2 = RGBColor[0.13, 0.36, 0.62];
sty = Sequence[Frame -> True, FrameStyle -> Black, Axes -> False,
   GridLines -> Automatic, GridLinesStyle -> GrayLevel[0.85],
   LabelStyle -> Directive[FontFamily -> "Latin Modern Roman", 13, Black],
   ImageSize -> 480, PlotRangePadding -> Scaled[0.03]];

(* Fig 1: the NLO decay rate vs Delta at mu = mhat (the mu-band is < 0.15% on
   this absolute scale -- invisible here, shown zoomed in Fig 2) *)
rateVal[e_] := -e["cen"] e["Dl"]^3;             (* -ln(Gamma/V) Delta^3 (positive, O(1)) *)
central = {#["Dl"], rateVal[#]} & /@ scan;
figRate = ListLinePlot[central,
   PlotStyle -> Directive[colC, AbsoluteThickness[1.8]],
   PlotMarkers -> {Graphics[{colC, Disk[]}], 0.022}, sty, PlotRange -> All,
   FrameLabel -> {"\[CapitalDelta]", "-ln(\[CapitalGamma]/V)\[ThinSpace]\[CapitalDelta]\!\(\*SuperscriptBox[\(\), \(3\)]\)"},
   Epilog -> {Text[Style["\[Mu] = \!\(\*OverscriptBox[\(m\), \(^\)]\)  (NLO)", 12, colC], Scaled[{0.80, 0.87}]]}];

(* Fig 2 (the money plot): fractional mu-variation band in percent vs Delta.
   1-loop running -> the NLO mu-residual (Baratella cancellation, tiny);
   2-loop running -> the N2LO-size estimate (the uncancelled S0 running exposes
   the missing two-loop determinant) *)
pct[e_, k_] := 100 (e[k] - e["cen"])/Abs[e["cen"]];   (* {lo%, hi%} band edges *)
b1lo = {#["Dl"], pct[#, "b1"][[1]]} & /@ scan;  b1hi = {#["Dl"], pct[#, "b1"][[2]]} & /@ scan;
b2lo = {#["Dl"], pct[#, "b2"][[1]]} & /@ scan;  b2hi = {#["Dl"], pct[#, "b2"][[2]]} & /@ scan;
figBand = Show[
   ListLinePlot[{b2lo, b2hi}, PlotStyle -> {{col2, AbsoluteThickness[1.3]}, {col2, AbsoluteThickness[1.3]}},
      Filling -> {1 -> {2}}, FillingStyle -> Directive[col2, Opacity[0.35]], InterpolationOrder -> 2],
   ListLinePlot[{b1lo, b1hi}, PlotStyle -> {{col1, AbsoluteThickness[1.3]}, {col1, AbsoluteThickness[1.3]}},
      Filling -> {1 -> {2}}, FillingStyle -> Directive[col1, Opacity[0.6]], InterpolationOrder -> 2],
   Graphics[{GrayLevel[0.4], AbsoluteThickness[1], Dashing[{0.02, 0.015}],
      Line[{{First[grid], 0}, {Last[grid], 0}}]}],
   sty, PlotRange -> {{First[grid], Last[grid]}, {-0.17, 0.17}},
   FrameLabel -> {"\[CapitalDelta]",
     "\[Delta] ln(\[CapitalGamma]/V) / |ln(\[CapitalGamma]/V)|   [%]"},
   Epilog -> {
     Text[Style["2-loop running", 12, col2], Scaled[{0.68, 0.90}]],
     Text[Style["\!\(\*SuperscriptBox[\(N\), \(2\)]\)LO-size estimate", 10, col2], Scaled[{0.68, 0.83}]],
     Text[Style["1-loop running (Baratella)", 10.5, RGBColor[0.30, 0.42, 0.20]], Scaled[{0.60, 0.585}]]}];

outDir = If[$InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
Export[FileNameJoin[{outDir, "RGImprovedRate_rate.pdf"}], figRate];
Export[FileNameJoin[{outDir, "RGImprovedRate_rate.png"}], figRate, ImageResolution -> 150];
Export[FileNameJoin[{outDir, "RGImprovedRate_band.pdf"}], figBand];
Export[FileNameJoin[{outDir, "RGImprovedRate_band.png"}], figBand, ImageResolution -> 150];
Print["wrote RGImprovedRate_{rate,band}.{pdf,png} in ", outDir];
