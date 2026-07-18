(* PolyDet` test suite (.wlt) -- VerificationTests across all sectors.
   Run via Tests/run.wls, or:  TestReport[".../Tests/PolyDet.wlt"]["AllTestsSucceeded"].
   Self-contained: registers + loads the paclet from the enclosing mathematica/ dir.
   Mostly cross-sector consistency checks (no magic numbers) + a few numeric anchors. *)

PacletDirectoryLoad[ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
Needs["PolyDet`"];

close[a_, b_, tol_: 1.*^-9] := TrueQ[Abs[a - b] <= tol Max[1, Abs[b]]];
m2 = {1.0, -0.5, 2.0, 1.3}; rr = {3.0, 5.0, 7.0};   (* generic, one tachyonic shell *)

(* ---- load / context ---- *)
VerificationTest[Context[PolyDetLnR], "PolyDet`", TestID -> "context"];
VerificationTest[Length@Select[
    {PolyDetLnR, PolyDetSigma, PolyDetRen, PolyDetIp, PolyDetInputs, PolyDetXMax,
     PolyDetLnRMat, PolyDetSigmaMat, PolyDetRenMat, PolyDetR1Prime, PolyDetRate,
     PolyDetLnRPsi, PolyDetSigmaPsi, PolyDetRenPsi, PolyDetIpPsi,
     PolyDetLnRGaugeT, PolyDetLnRGaugeDa, PolyDetIpGauge, PolyDetSigmaGauge, PolyDetInputsGauge,
     PolyDetSigmaGaugeNonAbelian},
    (DownValues[#] === {} && Options[#] === {}) &], 0, TestID -> "all-symbols-defined"];

(* ---- single scalar ---- *)
VerificationTest[close[PolyDetLnR[3, m2, rr], -1.37853663222931439], True, TestID -> "scalar-lnR"];
VerificationTest[PolyDetIp[m2, rr, 4][[2]], 0, TestID -> "scalar-I2-vanishes"];
VerificationTest[close[PolyDetSigma[m2, rr, "Dimension" -> 4, "Nsum" -> 80, WorkingPrecision -> 50],
   300.54897022247314, 1.*^-6], True, TestID -> "scalar-Sigma4"];
VerificationTest[close[Abs[Det[PolyDet`Private`Ts[4, 3.0, Sqrt[1.0], Sqrt[2.0]]] - 1], 0, 1.*^-12],
   True, TestID -> "scalar-detTs-eq-1"];

(* ---- multiscalar: diagonal mass matrix == sum of the two scalar channels ---- *)
With[{a = {1.0, 0.8, 1.4, 1.1}, b = {0.6, 1.2, 0.9, 1.3}},
  mm = MapThread[DiagonalMatrix[{#1, #2}] &, {a, b}];
  VerificationTest[close[Re@PolyDetLnRMat[3, mm, rr], Re[PolyDetLnR[3, a, rr] + PolyDetLnR[3, b, rr]], 1.*^-10],
     True, TestID -> "matrix-diag-eq-scalar-sum"];
  VerificationTest[close[Total@Abs[PolyDetIpMat[mm, rr, 4] - (PolyDetIp[a, rr, 4] + PolyDetIp[b, rr, 4])], 0, 1.*^-10],
     True, TestID -> "matrix-Ip-diag-eq-trace"]];

(* ---- Dirac fermion ---- *)
MM = {-1.0, 0.4, 0.9, 1.1};                          (* signed shell masses *)
VerificationTest[close[PolyDetLnRPsi[2, MM, rr], -3.6559263223125296], True, TestID -> "fermion-lnRpsi"];
VerificationTest[close[PolyDetIpPsi[MM, rr, 1][[1]],
   Sum[((MM[[s]] rr[[s]])^2 - (MM[[s + 1]] rr[[s]])^2)/2, {s, 3}], 1.*^-10], True, TestID -> "fermion-I1psi-closed"];
VerificationTest[Abs[PolyDetIpPsi[MM, rr, 2][[2]]] > 0.01, True, TestID -> "fermion-I2psi-nonzero"];

(* ---- Abelian gauge boson ---- *)
mA2 = {0.5, 1.0, 0.8, 1.1}; ma2 = mA2 + 0.3; dmA = {0.1, -0.05, 0.07};
VerificationTest[PolyDetLnRGaugeT[3, mA2, rr], PolyDetLnR[4, mA2, rr], TestID -> "gauge-transverse-eq-scalar"];
VerificationTest[close[
   PolyDet`Private`gaugeLnRDaCore[3, mA2, ma2, dmA, rr, 90, False],
   PolyDetLnR[3, mA2, rr, 90] + PolyDetLnR[5, mA2, rr, 90] + PolyDetLnR[4, ma2, rr, 90], 1.*^-10],
   True, TestID -> "gauge-decoupling-eq-3-scalars"];
VerificationTest[close[Total@Abs[PolyDetIpGauge[mA2, ma2, rr, 5] -
   (2 PolyDetIp[mA2, rr, 5] + PolyDetIp[ma2, rr, 5])], 0, 1.*^-12], True, TestID -> "gauge-Ieff-trace"];
VerificationTest[close[PolyDetLnRGaugeDa[3, mA2, ma2, dmA, rr], -1.87056187091516764, 1.*^-9],
   True, TestID -> "gauge-Da-coupled"];
(* general-xi chi-reduction: Xi->1 reduces EXACTLY to the default diagonal block *)
VerificationTest[close[PolyDetLnRGaugeDa[3, mA2, ma2, dmA, rr, "Xi" -> 1],
   PolyDetLnRGaugeDa[3, mA2, ma2, dmA, rr], 1.*^-12], True, TestID -> "gauge-xi1-reduces"];
VerificationTest[close[PolyDet`Private`gaugeLnRDaXi[3, mA2, ma2, dmA, rr, 1, 90],
   PolyDet`Private`gaugeLnRDaCore[3, mA2, ma2, dmA, rr, 90, True], 1.*^-12],
   True, TestID -> "gauge-chi-reduced-xi1"];
(* DIRECT general-xi diagonal block (supersedes the O(1)-wrong chi-reduced route for
   xi != 1; validated vs python/verify_gauge_xi_direct.py, 7/7).  dmA = Differences[Sqrt[mA2]]
   is the physical g Delta phi jump the direct chain assumes; ma2 = mA2 + Vaa is the xi=1
   Goldstone mass, internally shifted to xi mA2 + Vaa. *)
mA2x = {1/2, 1, 4/5, 11/10}; ma2x = mA2x + 3/10; dmAx = Differences[Sqrt[mA2x]];
VerificationTest[close[PolyDetLnRGaugeDa[2, mA2x, ma2x, dmAx, rr, "Xi" -> 2, WorkingPrecision -> 60],
   -3.20744623571243379, 1.*^-9], True, TestID -> "gauge-xi-direct-j2"];
VerificationTest[close[PolyDetLnRGaugeDa[0, mA2x, ma2x, dmAx, rr, "Xi" -> 2, WorkingPrecision -> 60],
   -3.07349489277497169, 1.*^-9], True, TestID -> "gauge-xi-direct-j0"];
(* the direct chain's xi=1 limit == the validated gaugeLnRDaCore (independent machinery) *)
VerificationTest[close[PolyDet`Private`gaugeLnRDaXiDirect[3, mA2x, ma2x, dmAx, rr, 1, 90],
   PolyDet`Private`gaugeLnRDaCore[3, mA2x, ma2x, dmAx, rr, 90, True], 1.*^-10],
   True, TestID -> "gauge-xi-direct-xi1-limit"];
(* general-xi ASSEMBLY (item 4): the running coefficients from PolyDetSigmaGauge Xi != 1
   (direct-chain summand + effective-trace subtraction; validated in
   python/verify_gauge_xi_assembly.py + notebooks/verify_gauge_xi_assembly_paclet.wls).
   I_1^eff is xi-INDEPENDENT (gauge-invariant quadratic divergence); I_3^eff = the analytic
   effective trace 3 I_3(mA^2) - I_3(xi mA^2) + I_3(xi mA^2 + Vaa).  NB the bare Sigma(xi != 1)
   is NOT ns-convergent (the (3-xi)/2 moment gradient leaves the c_0 drift) -- present the
   running coefficients, not the bare Sigma. *)
VerificationTest[close[
   PolyDetSigmaGauge[mA2x, ma2x, dmAx, rr, "Xi" -> 2, "Nsum" -> 30, "Ig" -> 0]["I1eff"],
   PolyDetSigmaGauge[mA2x, ma2x, dmAx, rr, "Xi" -> 1, "Nsum" -> 30, "Ig" -> 0]["I1eff"], 1.*^-10],
   True, TestID -> "gauge-xi-assembly-I1eff-xi-indep"];
VerificationTest[close[
   PolyDetSigmaGauge[mA2x, ma2x, dmAx, rr, "Xi" -> 2, "Nsum" -> 30, "Ig" -> 0]["I3eff"],
   3 PolyDetIp[mA2x, rr, 3][[3]] - PolyDetIp[2 mA2x, rr, 3][[3]] + PolyDetIp[2 mA2x + 3/10, rr, 3][[3]],
   1.*^-9], True, TestID -> "gauge-xi-assembly-I3eff-formula"];
(* general-xi gradient: Ig -> Automatic uses the closed-form (1+xi)/2 gaugeGradLead pair
   sum (GradientMethod -> Closed), making Sigma_gauge(xi != 1) ns-convergent (the c_0(xi) =
   (1+xi)/2 c_0(1) self-energy is removed; the old moment left a drift). *)
VerificationTest[
   PolyDetSigmaGauge[mA2x, ma2x, dmAx, rr, "Xi" -> 2, "Nsum" -> 30]["GradientMethod"],
   "Closed", TestID -> "gauge-xi-gradient-closed"];
VerificationTest[With[
   {s = Table[PolyDetSigmaGauge[mA2x, ma2x, dmAx, rr, "Xi" -> 2, "Nsum" -> ns]["Sigma"], {ns, {30, 45, 60}}]},
   Abs[s[[3]] - s[[2]]] < Abs[s[[2]] - s[[1]]] && Abs[s[[3]] - s[[2]]] < 3],
   True, TestID -> "gauge-xi-sigma-convergent"];

(* ---- a genuine, self-consistent U(1) Abelian-Higgs two-field example ----
   A smooth double well in sigma = phi^2/2 -> the Goldstone curvature V_aa = V'/phi vanishes at the
   vacua and the gauge m_A^2 = g^2 phi^2 is massive in both.  Anchored to the structural two-field
   requirements + the closed-form gauge-coupling coefficient I_3^eff (PolyDetIpGauge is closed-form,
   no multipole sum, so this stays fast). *)
Needs["FindBounce`"];
UuU[s_] := 6 (s - 2/25)^2 (s - 18/25)^2 - (3/25) (s - 2/25);
VuU[ph_] := UuU[ph^2/2];
{fvU, tvU} = MinMax[ph /. Solve[{VuU'[ph] == 0, ph > 0}, ph, Reals]];
bfU = FindBounce[VuU[ph], ph, {fvU, tvU}];
{mA2U, ma2U, dmAU, rrU} = PolyDetInputsGauge[bfU, 7/10, VuU'[#]/# &];
VerificationTest[Min[mA2U[[{1, -1}]]] > 0, True, TestID -> "gauge-u1-massive-vacua"];
VerificationTest[close[Max[Abs[(ma2U - mA2U)[[{1, -1}]]]], 0, 1.*^-8], True, TestID -> "gauge-u1-goldstone-flat-vacua"];
VerificationTest[Min[ma2U - mA2U] < 0, True, TestID -> "gauge-u1-goldstone-tachyonic-wall"];
VerificationTest[close[PolyDetIpGauge[mA2U, ma2U, rrU, 3][[3]], -90.0457, 1.*^-3], True, TestID -> "gauge-u1-I3eff"];

(* ---- generic-D ("Dimension" option): D=3 fermion + gauge ---- *)
(* the D=3 chain runs at HALF-INTEGER orders; anchors = independently derived reference values *)
VerificationTest[close[PolyDetLnRPsi[5/2, {13/10, 7/10, 11/10}, {1, 2}, 50],
   -0.26629332652599305, 1.*^-10], True, TestID -> "fermion-D3-halfinteger-lnRpsi"];
mA2d = {0.5, 1.0, 0.8, 1.1}; ma2d = mA2d + 0.3; dmAd = {0.1, -0.05, 0.07};
VerificationTest[close[
   PolyDetLnRGaugeDa[2, mA2d, ma2d, 0 dmAd, {1., 2., 3.}, "Dimension" -> 3, WorkingPrecision -> 90],
   PolyDetLnR[3/2, mA2d, {1., 2., 3.}, 90] + PolyDetLnR[7/2, mA2d, {1., 2., 3.}, 90] +
   PolyDetLnR[5/2, ma2d, {1., 2., 3.}, 90], 1.*^-10],
   True, TestID -> "gauge-D3-decoupling-halfinteger"];
VerificationTest[
   PolyDetSigmaGauge[mA2d, ma2d, dmAd, {1., 2., 3.}, "Dimension" -> 3, "Nsum" -> 20]["GradientMethod"],
   "ClosedPairSum", TestID -> "gauge-D3-closed-pairsum-scheme"];

(* ---- massless false vacuum (symmetric FV, Coulomb phase) ---- *)
(* single matching: closed form R_nu = Gamma(nu) (2/z)^(nu-1) I_{nu-1}(z), z = m_1 R_1 *)
VerificationTest[With[{z = Sqrt[2.] 5/2},
   close[PolyDetLnR[3, {2., 0}, {5/2}], Log[Gamma[3] (2/z)^2 BesselI[2, z]], 1.*^-12]],
   True, TestID -> "scalar-masslessFV-closedform"];
(* anchors: independently computed (mpmath) massless-FV values on the generic chain *)
VerificationTest[close[PolyDetLnR[3, {1., -0.5, 2., 0}, rr], 2.36945726848231406, 1.*^-10],
   True, TestID -> "scalar-masslessFV-anchor"];
(* == the m_FV -> 0 limit of the massive chain *)
VerificationTest[close[PolyDetLnR[2, {1., -0.5, 2., 0}, rr, 40],
   PolyDetLnR[2, {1., -0.5, 2., 1.*^-16}, rr, 40], 1.*^-10],
   True, TestID -> "scalar-masslessFV-limit"];
(* vectorised band == single-nu chain on a massless-FV input *)
VerificationTest[close[PolyDet`Private`lnRall[{1., -0.5, 2., 0}, rr, 4, 6, 40][[4]],
   PolyDetLnR[4, {1., -0.5, 2., 0}, rr, 40], 1.*^-12],
   True, TestID -> "scalar-masslessFV-lnRall"];
(* Laurent coefficients stay finite at a zero FV mass (the x^0 -> 1 guard) *)
VerificationTest[FreeQ[PolyDetIp[{1., -0.5, 2., 0}, rr, 6], Indeterminate | ComplexInfinity],
   True, TestID -> "scalar-masslessFV-Ip"];
(* gauge Da block with a massless FV: mpmath anchors (j >= 1 and the j = 0 two-channel) *)
mA2m = {0.5, 1.0, 0.8, 0}; ma2m = {0.8, 1.3, 1.1, 0.3}; dmAm = Differences[Sqrt[mA2m]];
VerificationTest[close[PolyDetLnRGaugeDa[3, mA2m, ma2m, dmAm, rr], 1.14985821563077864, 1.*^-9],
   True, TestID -> "gauge-masslessFV-Da"];
VerificationTest[close[PolyDetLnRGaugeDa[0, mA2m, ma2m, dmAm, rr], 5.54037037682660226, 1.*^-9],
   True, TestID -> "gauge-masslessFV-Da-j0"];
VerificationTest[close[
   PolyDet`Private`gaugeLnRDaCore[2, mA2m, ma2m, 0 dmAm, rr, 90, False],
   PolyDetLnR[2, mA2m, rr, 90] + PolyDetLnR[4, mA2m, rr, 90] + PolyDetLnR[3, ma2m, rr, 90], 1.*^-10],
   True, TestID -> "gauge-masslessFV-decoupling"];
VerificationTest[NumericQ[PolyDetSigmaGauge[mA2m, ma2m, dmAm, rr, "Nsum" -> 12]["Sigma"]],
   True, TestID -> "gauge-masslessFV-sigma-numeric"];
(* massless FV + general xi (Coulomb phase, symmetric FV): the direct diagonal block uses
   the power-law FV boundary gaugeEbndXiML (indicial eigenvectors, prefML=0).  Anchored to
   the mpmath massless chain + validated == the massive m_FV->0 limit
   (python/verify_gauge_xi_massless.py).  dmAm = Differences[Sqrt[mA2m]] (physical g Delta phi). *)
VerificationTest[close[PolyDetLnRGaugeDa[2, mA2m, ma2m, dmAm, rr, "Xi" -> 2, WorkingPrecision -> 90],
   9.996518390583735, 1.*^-9], True, TestID -> "gauge-masslessFV-xi-Da-j2"];
VerificationTest[close[PolyDetLnRGaugeDa[4, mA2m, ma2m, dmAm, rr, "Xi" -> 2, WorkingPrecision -> 90],
   4.962394641908676, 1.*^-9], True, TestID -> "gauge-masslessFV-xi-Da-j4"];
(* the massless-FV general-xi block's xi=1 limit == the validated gaugeLnRDaCore massless *)
VerificationTest[close[PolyDetLnRGaugeDa[2, mA2m, ma2m, dmAm, rr, "Xi" -> 1, WorkingPrecision -> 90],
   PolyDet`Private`gaugeLnRDaCore[2, mA2m, ma2m, dmAm, rr, 90, True], 1.*^-10],
   True, TestID -> "gauge-masslessFV-xi1-limit"];
(* the massless-FV general-xi assembly runs, is finite, and I_1^eff stays xi-independent *)
VerificationTest[close[
   PolyDetSigmaGauge[mA2m, ma2m, dmAm, rr, "Xi" -> 2, "Nsum" -> 30]["I1eff"],
   PolyDetSigmaGauge[mA2m, ma2m, dmAm, rr, "Xi" -> 1, "Nsum" -> 30]["I1eff"], 1.*^-10],
   True, TestID -> "gauge-masslessFV-xi-I1eff-indep"];
(* ORIENTATION zero removal ("ZeroRemoval" -> True): the j=0 term is replaced by the
   reduced determinant R_0' (mass-offset derivative, jumps fixed); the shift
   Log[Abs[R_0'/mhat_a^2]] - lndet R^(Da)_0 is anchored to the mpmath value *)
VerificationTest[close[
   PolyDetSigmaGauge[mA2m, ma2m, dmAm, rr, "Nsum" -> 12, "ZeroRemoval" -> True]["Sigma"] -
   PolyDetSigmaGauge[mA2m, ma2m, dmAm, rr, "Nsum" -> 12]["Sigma"],
   1.45497186506350, 1.*^-6], True, TestID -> "gauge-masslessFV-zeroremoval-shift"];
(* massless FV + xi != 1 is now handled by the direct chain's power-law FV boundary
   (gaugeEbndXiML); it returns a finite value (no longer guarded).  Anchored below in the
   gauge-masslessFV-xi-* tests. *)
VerificationTest[NumericQ[PolyDetLnRGaugeDa[2, mA2m, ma2m, dmAm, rr, "Xi" -> 3/2]],
   True, TestID -> "gauge-masslessFV-xi-runs"];

(* ---- FrozenWall thin-wall fast band ---- *)
(* synthetic case-b chain: TV plateau to R=15, three wall shells, FV outside; exact
   rational inputs so the deep-Pmax subtraction bookkeeping is exact.  FrozenWall
   must agree with the exact chain well below the 1e-3 target at x_max ~ 21. *)
With[{m2f = {2, -1/2, 1/4, 3/4, 1}, rf = {15, 16, 17, 18}},
  VerificationTest[close[
     PolyDetSigma[m2f, rf, "Target" -> 10^-3, "Method" -> "FrozenWall"],
     PolyDetSigma[m2f, rf, "Target" -> 10^-3], 5.*^-3 Abs[PolyDetSigma[m2f, rf, "Target" -> 10^-3]]],
     True, TestID -> "scalar-frozenwall-matches-exact"];
  VerificationTest[
     PolyDetSigma[m2f, rf, "Target" -> 10^-3, "Method" -> "FrozenWall", "NSub" -> 2] =!=
     PolyDetSigma[m2f, rf, "Target" -> 10^-3, "Method" -> "FrozenWall", "NSub" -> 16],
     True, TestID -> "scalar-frozenwall-nsub-active"]];

(* ---- "SubtractFrom" offset (Hurwitz) assembly: exact value identity ---- *)
With[{m2f = {2, -1/2, 1/4, 3/4, 1}, rf = {15, 16, 17, 18}},
  VerificationTest[
     Abs[PolyDetSigma[m2f, rf, "Nsum" -> 60, WorkingPrecision -> 50, "SubtractFrom" -> 22] -
         PolyDetSigma[m2f, rf, "Nsum" -> 60, WorkingPrecision -> 50, "SubtractFrom" -> 0]] < 10^-30,
     True, TestID -> "scalar-subtractfrom-identity"]];

(* ---- FrozenWall sector ports: fermion / gauge / matrix ---- *)
With[{Mf = {-13/10, -4/5, -1/5, 2/5, 1}, rf = {15, 16, 17, 18}},
  VerificationTest[
     Abs[PolyDetSigmaPsi[Mf, rf, "Target" -> 10^-2, "Method" -> "FrozenWall"]/
         PolyDetSigmaPsi[Mf, rf, "Target" -> 10^-2] - 1] < 5 10^-3,
     True, TestID -> "fermion-frozenwall-sigma"]];
With[{mA2 = {6/5, 4/5, 1/2, 9/10, 1}, ma2 = {6/5, 1/2, -1/2, 7/10, 1}, dmA = {-1/5, -3/10, 1/5, 1/4},
      rg = {15, 16, 17, 18}},
  VerificationTest[
     Abs[PolyDet`Private`gaugeLnRDaFrozen[6, mA2, ma2, dmA, rg, 8] -
         PolyDetLnRGaugeDa[6, mA2, ma2, dmA, rg]] < 10^-3,
     True, TestID -> "gauge-frozenwall-da"]];
With[{m2c = {{{2, 3/10}, {3/10, 3/2}}, {{6/5, 1/2}, {1/2, 2/5}}, {{-3/10, 2/5}, {2/5, 3/5}},
      {{3/5, 1/5}, {1/5, 9/10}}, {{1, 1/10}, {1/10, 4/5}}}, rm = {15, 16, 17, 18}},
  VerificationTest[                                   (* partial sum, equal Nsum both sides:
     isolates the frozen band; a full exact Sigma reference is minutes-slow (the very
     matrix-WP stress case the frozen port removes) *)
     Abs[PolyDetSigmaMat[m2c, (3/5) rm, "Nsum" -> 16, "Method" -> "FrozenWall"]/
         PolyDetSigmaMat[m2c, (3/5) rm, "Nsum" -> 16] - 1] < 10^-3,
     True, TestID -> "matrix-frozenwall-sigma"]];
VerificationTest[close[(3 + 1)^2 PolyDet`Private`gaugeGradPair[3, {0.1, -0.05, 0.07}, {1., 2., 3.}, 4],
   PolyDet`Private`gaugeGradLead[3, {0.1, -0.05, 0.07}, {1., 2., 3.}], 1.*^-14],
   True, TestID -> "gauge-gradpair-D4-reduction"];

(* ---- fermion gradient moment PolyDetIgPsi (wave-function running) ---- *)
VerificationTest[
  Module[{V, sp, bf, ig},
    V[p_] := (p^2 - 1)^2/8 + (1/20) (p - 1);
    sp = Sort[Re[x /. NSolve[x^3 - x + 2/20 == 0, x]]];
    bf = FindBounce[V[p], p, {N@sp[[3]], N@sp[[1]]}, FieldPoints -> 21, Dimension -> 4];
    ig = PolyDetIgPsi[bf, (6/10) # &];
    NumericQ[ig] && ig < 0 && Abs[ig - (6/10)^2 PolyDetIgPsi[bf, # &]] < 10^-6],  (* I_g^psi(y) = y^2 I_g^psi(1) *)
  True, TestID -> "fermion-Igpsi-scaling"];
