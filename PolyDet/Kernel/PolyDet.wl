(* ::Package:: *)

(* PolyDet` -- one-loop fluctuation determinant around a polygonal bounce.
   Companion to FindBounce: CONSUMES a polygonal bounce (segment masses + matching
   radii) and returns the renormalised determinant.

   This file is the LOADER: it declares the public PolyDet` context + usage
   messages, then Get's the implementation sections from Kernel/ (each defines into
   PolyDet`Private`).  Load order matters: Scalar provides helpers reused by the
   Multiscalar and Gauge sections.  Unit tests in Tests/PolyDet.wlt. *)

BeginPackage["PolyDet`"]

PolyDetLnR::usage =
  "PolyDetLnR[nu, masses2, radii] gives the Gel'fand-Yaglom log-ratio \
ln \!\(\*SubscriptBox[\(R\), \(\[Nu]\)]\)(\[Infinity]) = ln[(\!\(\*SubscriptBox[\(m\), \(N\)]\)/\!\(\*SubscriptBox[\(m\), \(1\)]\))^\[Nu] \!\(\*SubsuperscriptBox[\(T\), \(11\), \(tot\)]\)] for the polygonal chain with segment \
squared-masses masses2 = {\!\(\*SubsuperscriptBox[\(m\), \(1\), \(2\)]\),...,\!\(\*SubsuperscriptBox[\(m\), \(N\), \(2\)]\)} (signed; tachyonic \!\(\*SuperscriptBox[\(m\), \(2\)]\)<0 allowed) and \
matching radii = {\!\(\*SubscriptBox[\(R\), \(1\)]\),...,\!\(\*SubscriptBox[\(R\), \(N - 1\)]\)}. PolyDetLnR[nu, masses2, radii, prec] sets the \
working precision (default 60). A MASSLESS false vacuum (masses2[[-1]] == 0, symmetric FV) is \
supported: the ratio is then normalized against the massless free operator, with the prefactor \
LogGamma[\[Nu]+1] + \[Nu] Log[2/\!\(\*SubscriptBox[\(m\), \(1\)]\)] replacing the massive one.";

PolyDetSigma::usage =
  "PolyDetSigma[masses2, radii] gives the *subtracted* (UV-regulated, finite) \
fluctuation-determinant sum \!\(\*SubscriptBox[\(\[CapitalSigma]\), \(D\)]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(l\)]\) \!\(\*SubscriptBox[\(d\), \(l\)]\) (ln \!\(\*SubscriptBox[\(R\), \(\[Nu]\)]\) - \!\(\*SubscriptBox[\(\[Sum]\), \(p \[LessEqual] Pmax\)]\) \!\(\*SubscriptBox[\(I\), \(p\)]\)/\!\(\*SuperscriptBox[\(\[Nu]\), \(p\)]\)) + \!\(\*SubscriptBox[\(\[Sum]\), \(p \[GreaterEqual] D\)]\) \!\(\*SubscriptBox[\(I\), \(p\)]\) \!\(\*SubscriptBox[\(Z\), \(D\)]\)(p) (via optimal truncation). This is NOT yet the \
renormalised determinant -- use PolyDetRen to add the regulated divergent orders. \
Options: \"Dimension\" (4), WorkingPrecision (60), \"Pmax\" (11), \"Nsum\" \
(Automatic, the number of l-modes), \"Target\" (Automatic). \"Nsum\" overrides; \
else \"Target\"->eps sets the cutoff \!\(\*SubscriptBox[\(\[Nu]\), \(max\)]\) ~ \!\(\*SubscriptBox[\(x\), \(max\)]\) eps^(-1/(Pmax+2-D)) \
(e.g. eps=0.01 -> ~1.6 \!\(\*SubscriptBox[\(x\), \(max\)]\) at Pmax=11, ~4x faster); Automatic -> ~5 \!\(\*SubscriptBox[\(x\), \(max\)]\) (~1e-8). \
\"Method\"->\"FrozenWall\" enables the thin-wall fast band for case-b chains (TV plateau + wall): \
l >= 2 via exact plateau/FV Bessels + frozen-centrifugal wall exponentials on \"NSub\" (8) \
sub-intervals per shell, at machine precision; the error falls as ~1/\!\(\*SubsuperscriptBox[\(x\), \(max\), \(3\)]\) \
(2e-4 at \!\(\*SubscriptBox[\(x\), \(max\)]\)=23, 3e-6 at 87; 6-30x faster). \"Method\"->Automatic \
switches it on when case-b-like, \!\(\*SubscriptBox[\(x\), \(max\)]\)>=30 and a numeric \"Target\">=1e-3. \
\"SubtractFrom\"->\[Nu]* applies the Laurent subtraction only at \[Nu] >= \[Nu]* with a Hurwitz-\[Zeta] \
add-back and explicit harmonic-number compensation of the divergent orders below -- the returned value is \
IDENTICAL for every \[Nu]* (zeta-split identity), but with \[Nu]* ~ \!\(\*SubscriptBox[\(x\), \(max\)]\) no \
intermediate exceeds ~\!\(\*SubsuperscriptBox[\(x\), \(max\), \(D\)]\), making the assembly machine-safe and \
lifting the high-precision-inputs requirement (Automatic: \[Nu]* = \!\(\*SubscriptBox[\(x\), \(max\)]\) when \
\!\(\*SubscriptBox[\(x\), \(max\)]\) >= 30, else 0 = the classic full-range assembly).";

PolyDetRen::usage =
  "PolyDetRen[masses2, radii] gives the renormalised one-loop determinant \
ln \!\(\*SubscriptBox[\(det\), \(ren\)]\) = \!\(\*SubscriptBox[\(\[CapitalSigma]\), \(D\)]\) + \!\(\*SubscriptBox[\(\[Delta]\), \(ct\)]\)(\[Mu]). \
D=4 (MSbar): \!\(\*SubscriptBox[\(\[Delta]\), \(ct\)]\) = \!\(\*SubscriptBox[\(I\), \(3\)]\) (3/4 + EulerGamma + Log[\!\(\*OverscriptBox[\(\[Mu]\), \(~\)]\)/2]) + \!\(\*SubsuperscriptBox[\(I\), \(3\), \(ln\)]\) with \
\!\(\*SuperscriptBox[\(\*OverscriptBox[\(\[Mu]\), \(~\)]\), \(2\)]\) = 4 Pi \!\(\*SuperscriptBox[\(E\), \(-EulerGamma\)]\) \!\(\*SuperscriptBox[\(\[Mu]\), \(2\)]\) (the 1/eps pole, coeff \!\(\*SubscriptBox[\(I\), \(3\)]\), cancels against the \
action and is returned via \"Divergence\"->True). D=3 (super-renormalisable): finite, \
\[Mu]-independent, \!\(\*SubscriptBox[\(\[Delta]\), \(ct\)]\) = -\!\(\*OverscriptBox[\(m\), \(^\)]\) (1/3) \!\(\*SubscriptBox[\(\[Sum]\), \(s\)]\) \!\(\*SubsuperscriptBox[\(R\), \(s\), \(3\)]\) (\!\(\*SubsuperscriptBox[\(m\), \(s\), \(2\)]\)-\!\(\*SubsuperscriptBox[\(m\), \(\(s + 1\)\), \(2\)]\)). By default the \
l=1 translational zero modes are removed (\"ZeroRemoval\"->True): the artifact term \
\!\(\*SubscriptBox[\(d\), \(1\)]\) ln|\!\(\*SubscriptBox[\(R\), \(l = 1\)]\)| is replaced by \!\(\*SubscriptBox[\(d\), \(1\)]\) ln|\!\(\*SubsuperscriptBox[\(R\), \(1\), \(\[Prime]\)]\)/\!\(\*SuperscriptBox[\(\*OverscriptBox[\(m\), \(^\)]\), \(2\)]\)| (see PolyDetR1Prime), \!\(\*SubscriptBox[\(d\), \(1\)]\)=D. \
Options: \"Dimension\" (4), \"Scale\" (\[Mu], =1), \"ZeroRemoval\" (True), plus PolyDetSigma options.";

PolyDetIp::usage =
  "PolyDetIp[masses2, radii, pmax] gives the large-\[Nu] Laurent coefficients \
{\!\(\*SubscriptBox[\(I\), \(1\)]\),...,\!\(\*SubscriptBox[\(I\), \(pmax\)]\)} of ln \!\(\*SubscriptBox[\(R\), \(\[Nu]\)]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(p\)]\) \!\(\*SubscriptBox[\(I\), \(p\)]\)/\!\(\*SuperscriptBox[\(\[Nu]\), \(p\)]\), from the hardcoded universal \
per-segment forms \!\(\*SubscriptBox[\(g\), \(p\)]\) (\!\(\*SubscriptBox[\(I\), \(p\)]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(s\)]\) \!\(\*SubscriptBox[\(g\), \(p\)]\)(\!\(\*SubsuperscriptBox[\(m\), \(s\), \(2\)]\) \!\(\*SuperscriptBox[\(R\), \(2\)]\)/4, \!\(\*SubsuperscriptBox[\(m\), \(\(s + 1\)\), \(2\)]\) \!\(\*SuperscriptBox[\(R\), \(2\)]\)/4)). \!\(\*SubscriptBox[\(I\), \(2\)]\) = 0; \
pmax <= 14.";

PolyDetInputs::usage =
  "PolyDetInputs[Vpp, bf] extracts {masses2, radii} from a FindBounce BounceFunction \
bf: plateau masses Vpp[\[Phi]] at the two end vacua, polygonal (extended-polygonal) wall \
\!\(\*SuperscriptBox[\(masses\), \(2\)]\) = per-segment Simpson average of Vpp over consecutive bf[\"Path\"] nodes, radii \
from bf[\"Radii\"]. Vpp is V''(\[Phi]).";

PolyDetXMax::usage =
  "PolyDetXMax[masses2, radii] gives \!\(\*SubscriptBox[\(x\), \(max\)]\) = \!\(\*SubscriptBox[\(max\), \(s\)]\) |\!\(\*SubscriptBox[\(m\), \(s\)]\) \!\(\*SubscriptBox[\(R\), \(s\)]\)|, the largest Bessel \
argument (turning-point scale) entering the transfer chain.";

PolyDetR1Prime::usage =
  "PolyDetR1Prime[masses2, radii] gives the dimensionless reduced l=1 ratio \!\(\*SubsuperscriptBox[\(R\), \(1\), \(\[Prime]\)]\) \
that replaces the removed translational zero mode (\!\(\*SubscriptBox[\(R\), \(l = 1\)]\)=0) in the determinant, \
with degeneracy \!\(\*SubscriptBox[\(d\), \(1\)]\) = D.  Computed analytically as the \!\(\*SuperscriptBox[\(\[Mu]\), \(2\)]\)-derivative of the \
transfer chain (interior shell masses offset, FV fixed); -> \!\(\*SuperscriptBox[\(e\), \(D - 1\)]\)/12 in the \
thin-wall limit.  Option \"Dimension\" (4), WorkingPrecision (80).";

PolyDet::usage =
  "PolyDet[V, field, bf, opts] gives the renormalised one-loop determinant \
ln \!\(\*SubscriptBox[\(det\), \(ren\)]\) directly from a FindBounce BounceFunction bf for the potential V (in \
variable field): it extracts {masses2, radii} via PolyDetInputs and evaluates \
PolyDetRen, reading the spacetime dimension from bf[\"Dimension\"] (D=3 or 4). \
PolyDetRen/PolyDetSigma options (\"Scale\", \"Pmax\", \"Nsum\", WorkingPrecision) \
pass through. Needs the bounce built with FieldPoints >= 5 (wall masses ddVL).";

PolyDetLnRMat::usage =
  "PolyDetLnRMat[nu, masses2, radii] gives the matrix Gel'fand-Yaglom log-ratio \
ln det \!\(\*SubscriptBox[\(R\), \(\[Nu]\)]\) = (\[Nu]/2)(ln det \!\(\*SubsuperscriptBox[\(m\), \(N\), \(2\)]\) - ln det \!\(\*SubsuperscriptBox[\(m\), \(1\), \(2\)]\)) + ln det \!\(\*SubsuperscriptBox[\(T\), \(11\), \(tot\)]\) for \!\(\*SubscriptBox[\(N\), \(f\)]\) coupled \
scalars, each masses2[[s]] an \!\(\*SubscriptBox[\(N\), \(f\)]\) x \!\(\*SubscriptBox[\(N\), \(f\)]\) mass matrix (radii = {\!\(\*SubscriptBox[\(R\), \(1\)]\),...,\!\(\*SubscriptBox[\(R\), \(N - 1\)]\)}). The \
matrix transfer \!\(\*SuperscriptBox[\(T\), \(s\)]\) is built from Cayley-Hamilton matrix Bessels (eigenvalues only, \
no eigenvectors); \!\(\*SubscriptBox[\(det\), \(2 Nf\)]\) \!\(\*SuperscriptBox[\(T\), \(s\)]\) = 1. PolyDetLnRMat[nu, masses2, radii, prec] sets the \
precision (default 60). Reduces to a sum of scalar PolyDetLnR for diagonal masses.";

PolyDetIpMat::usage =
  "PolyDetIpMat[masses2, radii, pmax] gives the large-\[Nu] Laurent coefficients \
{\!\(\*SubscriptBox[\(I\), \(1\)]\),...,\!\(\*SubscriptBox[\(I\), \(pmax\)]\)} of ln det \!\(\*SubscriptBox[\(R\), \(\[Nu]\)]\) for coupled scalars: the scalar \!\(\*SubscriptBox[\(g\), \(p\)]\) under a flavour \
trace, \!\(\*SubscriptBox[\(I\), \(p\)]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(s\)]\) tr \!\(\*SubscriptBox[\(g\), \(p\)]\)(\!\(\*SubsuperscriptBox[\(m\), \(s\), \(2\)]\) \!\(\*SuperscriptBox[\(R\), \(2\)]\)/4, \!\(\*SubsuperscriptBox[\(m\), \(\(s + 1\)\), \(2\)]\) \!\(\*SuperscriptBox[\(R\), \(2\)]\)/4) (grouped matrix powers). The \
divergent orders \!\(\*SubscriptBox[\(I\), \(1\)]\) ~ tr \!\(\*SuperscriptBox[\(dm\), \(2\)]\), \!\(\*SubscriptBox[\(I\), \(3\)]\) ~ tr \!\(\*SuperscriptBox[\(dm\), \(4\)]\) are single-matrix powers, basis-independent.";

PolyDetSigmaMat::usage =
  "PolyDetSigmaMat[masses2, radii] gives the subtracted multipole sum \!\(\*SubscriptBox[\(\[CapitalSigma]\), \(D\)]\) for \
coupled scalars (matrix \!\(\*SubscriptBox[\(R\), \(\[Nu]\)]\) + trace subtraction / optimal truncation). Same options \
as PolyDetSigma. Options \"Method\" (\"Exact\" | \"FrozenWall\" | Automatic), \"NSub\" and \"SubtractFrom\" select the case-b thin-wall fast band and the offset (Hurwitz) assembly, as for PolyDetSigma.";

PolyDetRenMat::usage =
  "PolyDetRenMat[masses2, radii] gives the renormalised one-loop determinant for \
\!\(\*SubscriptBox[\(N\), \(f\)]\) coupled scalars: \!\(\*SubscriptBox[\(\[CapitalSigma]\), \(D\)]\) (matrix \!\(\*SubscriptBox[\(R\), \(\[Nu]\)]\)) + \!\(\*SubscriptBox[\(\[Delta]\), \(ct\)]\) from the flavour traces \
tr \!\(\*SuperscriptBox[\(dm\), \(2\)]\), tr \!\(\*SuperscriptBox[\(dm\), \(4\)]\) (D=4 MSbar; D=3 finite). By default the l=1 translational zero \
modes are removed (\"ZeroRemoval\"->True) via the matrix \!\(\*SubsuperscriptBox[\(R\), \(1\), \(\[Prime]\)]\) (PolyDetR1PrimeMat). \
Options as PolyDetRen.";

PolyDetR1PrimeMat::usage =
  "PolyDetR1PrimeMat[masses2, radii] gives the matrix translational \!\(\*SubsuperscriptBox[\(R\), \(1\), \(\[Prime]\)]\) for \!\(\*SubscriptBox[\(N\), \(f\)]\) \
coupled scalars: the reduced l=1 (\[Nu]=D/2) Gel'fand-Yaglom determinant that replaces \
the removed translational zero mode (det \!\(\*SubscriptBox[\(R\), \(l = 1\)]\)=0, one zero eigenvalue along the \
bounce tangent). Computed as the \!\(\*SuperscriptBox[\(\[Mu]\), \(2\)]\)-derivative of det \!\(\*SubscriptBox[\(R\), \(l = 1\)]\) via Jacobi/adjugate \
(matrix dTs source), normalised by the FV mass projected on ker(T11) (the zero-mode \
direction). Reduces to the scalar PolyDetR1Prime when fields decouple. Any N_f \
(general adjugate via Faddeev-LeVerrier; null vector from adj(T11) columns). Option \
\"Dimension\" (4), WorkingPrecision (80).";

PolyDetXMaxMat::usage =
  "PolyDetXMaxMat[masses2, radii] gives \!\(\*SubscriptBox[\(x\), \(max\)]\) = \!\(\*SubscriptBox[\(max\), \(s\)]\) \!\(\*SubscriptBox[\(R\), \(s\)]\) Sqrt[max eigenvalue of \!\(\*SubsuperscriptBox[\(m\), \(s\), \(2\)]\)] \
for the coupled (matrix) chain.";

PolyDetRate::usage =
  "PolyDetRate[lndetRen, S0] gives ln(\[CapitalGamma]/Volume), the log bounce-nucleation rate \
prefactor (Callan-Coleman): (D/2) Log[\!\(\*SubscriptBox[\(S\), \(0\)]\)/(2 Pi)] - lndetRen/2 - \!\(\*SubscriptBox[\(S\), \(0\)]\), with lndetRen = \
PolyDetRen/PolyDetRenMat output (zero modes removed = the reduced \!\(\*SuperscriptBox[\(det\), \(\[Prime]\)]\)/\!\(\*SubscriptBox[\(det\), \(FV\)]\)) and \!\(\*SubscriptBox[\(S\), \(0\)]\) the \
bounce action. Options: \"Dimension\" (4); for a broken global symmetry, \"Goldstones\"->nG \
with \"ChargeNorm\"->Q, the FULL D-dimensional zero-mode norm: Q = Integrate over d^D x of \
(T \[Phi])^2 = \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(D - 1\)]\) Integrate[\[Rho]^(D - 1) (T \[Phi])^2, \[Rho]] evaluated on the bounce -- NOT the radial \
integral alone (generator normalisation matching \"OrbitVolume\") -- and \"OrbitVolume\"->Vol(G/H) \
(default 2 Pi) add the orbit measure. This is the COMPLETE \
prefactor: the single negative mode (l=0 radial) sits in |\!\(\*SuperscriptBox[\(det\), \(\[Prime]\)]\)| (= lndetRen via Re[Log]), \
and Coleman's contour 1/2 cancels the 2 of \[CapitalGamma]=2 Im exactly -- no extra factor (assumes \
the standard one-negative-mode bounce).";

PolyDetDecayRate::usage =
  "PolyDetDecayRate[Vpp, bf, extra] gives ln(\[CapitalGamma]/Volume), the COMPLETE one-loop decay \
rate from a FindBounce BounceFunction bf for a model with driving (Higgs) scalar potential second \
derivative Vpp = V''(\[Phi]) plus a list of extra field sectors. It assembles \
ln \!\(\*SubscriptBox[\(det\), \(tot\)]\) = \!\(\*SuperscriptBox[\(\[CapitalSigma]\), \(Higgs\)]\) + \
Sum_scalars \[CapitalSigma] + Sum_gauge \!\(\*SubscriptBox[\(\[CapitalSigma]\), \(gauge\)]\) - Sum_fermions \!\(\*SuperscriptBox[\(\[CapitalSigma]\), \(\[Psi]\)]\) \
(the fermion -1 is the (-1)^F sign), then PolyDetRate. The Higgs carries the translational zero modes \
(removed). All sectors share one MS-bar \"Scale\" (\[Mu]) and the Simpson-averaged shell masses (as PolyDetInputs). Each \
extra sector is an Association: <|\"Type\"->\"Scalar\", \"Mass2\"->fn, \"Mult\"->n|> (m^2(\[Phi])); \
<|\"Type\"->\"Fermion\", \"Mass\"->fn, \"Mult\"->n|> (signed real M(\[Phi]); Mult = color/flavour); \
<|\"Type\"->\"Gauge\", \"Coupling\"->g, \"Vaa\"->fn, \"Mult\"->n|> (mA2 = g^2 \[Phi]^2, ma2 = mA2 + Vaa[\[Phi]], \
dmA = g \[CapitalDelta]\[Phi]). Options: \"Scale\" (\[Mu]), \"ZeroRemoval\" (True), \"Goldstones\"/\"ChargeNorm\"/\"OrbitVolume\" \
(broken global symmetry), \"Ig\" (gauge gradient: Automatic fit or a number), \"Breakdown\"->True \
returns <|Rate, lndetTotal, Higgs, Sectors, S0, Dimension|>.";

PolyDetInputs::noext =
  "BounceFunction carries no extension data (CoefficientsExtension); rebuild with \
FieldPoints >= 5 and Gradient =!= None so the wall masses ddVL are available.";

(* ---- Dirac fermion (Yukawa) sector ---- *)
PolyDetLnRPsi::usage =
  "PolyDetLnRPsi[nu, masses, radii] gives ln \!\(\*SubsuperscriptBox[\(R\), \(\[Nu]\), \(\[Psi]\)]\) = 2 ln|\!\(\*SubscriptBox[\(r\), \(\[Nu]\)]\)| for a single Dirac \
fermion on the polygonal bounce: SIGNED shell masses masses = {\!\(\*SubscriptBox[\(M\), \(1\)]\),...,\!\(\*SubscriptBox[\(M\), \(N\)]\)}, \!\(\*SubscriptBox[\(M\), \(s\)]\) = m + y \!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\) \
(NOT squared; tachyonic-free, sign-changing allowed), radii = {\!\(\*SubscriptBox[\(R\), \(1\)]\),...,\!\(\*SubscriptBox[\(R\), \(N - 1\)]\)}. First-order \
radial Dirac transfer chain with ADJACENT Bessel orders \[Nu], \[Nu]+1; \!\(\*SubscriptBox[\(r\), \(\[Nu]\)]\) = (|\!\(\*SubscriptBox[\(M\), \(N\)]\)|/|\!\(\*SubscriptBox[\(M\), \(1\)]\)|)^\[Nu] \
\!\(\*SubsuperscriptBox[\(T\), \(11\), \(tot\)]\). PolyDetLnRPsi[nu, masses, radii, prec] sets the working precision.";
PolyDetIpPsi::usage =
  "PolyDetIpPsi[masses, radii, pmax] gives {\!\(\*SubsuperscriptBox[\(I\), \(1\), \(\[Psi]\)]\),...,\!\(\*SubsuperscriptBox[\(I\), \(pmax\), \(\[Psi]\)]\)} (pmax<=15), the large-\[Nu] \
Laurent coefficients of ln \!\(\*SubsuperscriptBox[\(R\), \(\[Nu]\), \(\[Psi]\)]\), summed over the chain from hardcoded universal \
per-matching coefficients in SIGNED \!\(\*SubscriptBox[\(x\), \(s\)]\) = \!\(\*SubscriptBox[\(M\), \(s\)]\) \!\(\*SubscriptBox[\(R\), \(s\)]\). \!\(\*SubsuperscriptBox[\(I\), \(2\), \(\[Psi]\)]\) != 0 (the \[Nu](\[Nu]+1) degeneracy \
breaks even/odd).";
PolyDetSigmaPsi::usage =
  "PolyDetSigmaPsi[masses, radii] gives the MS-bar SUBTRACTED (finite) one-loop Dirac \
determinant \!\(\*SubsuperscriptBox[\(\[CapitalSigma]\), \(D\), \(\[Psi]\)]\): degeneracy \!\(\*SubsuperscriptBox[\(d\), \(\[Nu]\), \(\[Psi]\)]\) = \[Nu](\[Nu]+1) at D=4, zeta add-back \
\!\(\*SuperscriptBox[\(Z\), \(\[Psi]\)]\)(p) = \[Zeta](p-2)+\[Zeta](p-1) (convergent p>=4). \"Dimension\" (default 4; also 3): the chain and \
\!\(\*SubsuperscriptBox[\(I\), \(p\), \(\[Psi]\)]\) are D-independent, D enters via the HALF-INTEGER orders \[Mu] = l+D/2-1, the weight \[Mu]+1/2, \
and the Hurwitz add-back \[Zeta](p,1/2) (D=3 is finite and scale-free; the polygon's jump-log \
\!\(\*SubsuperscriptBox[\(I\), \(2\), \(\[Psi]\)]\)+\!\(\*SubsuperscriptBox[\(I\), \(1\), \(\[Psi]\)]\)/2 is the dropped O(1/N) artifact). Options: \"Pmax\" (7), \"Nsum\"/\"Target\", \
WorkingPrecision (as PolyDetSigma). Options \"Method\" (\"Exact\" | \"FrozenWall\" | Automatic), \"NSub\" and \"SubtractFrom\" select the case-b thin-wall fast band and the offset (Hurwitz) assembly, as for PolyDetSigma.";
PolyDetRenPsi::usage =
  "PolyDetRenPsi[masses, radii] gives the renormalised Dirac determinant \
\!\(\*SubsuperscriptBox[\(\[CapitalSigma]\), \(4\), \(\[Psi], ren\)]\)(\[Mu]) = \!\(\*SubsuperscriptBox[\(\[CapitalSigma]\), \(4\), \(\[Psi]\)]\) + \!\(\*SuperscriptBox[\(\[Beta]\), \(\[Psi]\)]\) Log[Scale], \!\(\*SuperscriptBox[\(\[Beta]\), \(\[Psi]\)]\) = \!\(\*SubsuperscriptBox[\(I\), \(2\), \(\[Psi]\)]\) + \!\(\*SubsuperscriptBox[\(I\), \(3\), \(\[Psi]\)]\) \
(both p=2,3 run, from the \!\(\*SuperscriptBox[\(\[Nu]\), \(1\)]\) and \!\(\*SuperscriptBox[\(\[Nu]\), \(2\)]\) parts of \!\(\*SubsuperscriptBox[\(d\), \(\[Nu]\), \(\[Psi]\)]\)). \"Scale\"->\[Mu] (default 1 = the \
MS-bar value). d/d ln \[Mu] = \!\(\*SuperscriptBox[\(\[Beta]\), \(\[Psi]\)]\). \"Dimension\"->3: scale-free (odd-D dim reg, no running; \
\"Scale\" is ignored).";
PolyDetXMaxPsi::usage =
  "PolyDetXMaxPsi[masses, radii] gives \!\(\*SubscriptBox[\(x\), \(max\)]\) = \!\(\*SubscriptBox[\(max\), \(s\)]\) |\!\(\*SubscriptBox[\(M\), \(s\)]\)| \!\(\*SubscriptBox[\(R\), \(s\)]\), the largest Bessel argument \
of the fermion chain (sets \!\(\*SubscriptBox[\(\[Nu]\), \(max\)]\)).";
PolyDetIgPsi::usage =
  "PolyDetIgPsi[bf, M] gives the fermion gradient moment \!\(\*SubsuperscriptBox[\(I\), \(g\), \(\[Psi]\)]\) = \
-\[Integral]\!\(\*SuperscriptBox[\(\[Rho]\), \(D-1\)]\)\!\(\*SuperscriptBox[\((\*FractionBox[\(d\), \(d\[Rho]\)] M(\*OverscriptBox[\(\[Phi]\), \(_\)]))\), \(2\)]\)d\[Rho] on the BounceFunction bf, \
M the signed fermion mass as a pure function of the field (M(\[Phi]) = m + y \[Phi]). This is the \
wave-function (gradient) running the first-order polygonal fermion determinant discards \
(\!\(\*SuperscriptBox[\(\[Beta]\), \(\[Psi]\)]\) is pure potential); the full running is \
d\!\(\*SuperscriptBox[\(\[CapitalSigma]\), \(\[Psi]\)]\)/d ln \[Mu] = \!\(\*SuperscriptBox[\(\[Beta]\), \(\[Psi]\)]\) + \!\(\*SubsuperscriptBox[\(I\), \(g\), \(\[Psi]\)]\), \
the fermion analogue of the gauge \!\(\*SubsuperscriptBox[\(I\), \(3\), \(eff\)]\) + \!\(\*SubscriptBox[\(I\), \(g\)]\). \
Fixed by the bridge \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Phi]\)]\)\[Integral]\!\(\*SuperscriptBox[\((\[PartialD]\*OverscriptBox[\(\[Phi]\), \(_\)])\), \(2\)]\) = -1/4 \!\(\*SubsuperscriptBox[\(I\), \(g\), \(\[Psi]\)]\), \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Phi]\)]\) = 2\!\(\*SuperscriptBox[\(y\), \(2\)]\)/16\!\(\*SuperscriptBox[\(\[Pi]\), \(2\)]\).";

(* ---- Abelian gauge boson (U(1) Higgs) sector ---- *)
PolyDetLnRGaugeT::usage =
  "PolyDetLnRGaugeT[j, mA2, radii] gives the transverse gauge ln \!\(\*SubsuperscriptBox[\(R\), \(j\), \((T)\)]\): the scalar \
Gel'fand-Yaglom ratio of the gauge mass mA2 = \!\(\*SuperscriptBox[\(g\), \(2\)]\) \!\(\*SuperscriptBox[\(\[Phi]\), \(2\)]\) at order \[Nu] = j+1 (orbital l=j); one \
polarization (degeneracy 2 j(j+2)). = PolyDetLnR[j+1, mA2, radii].";
PolyDetLnRGaugeDa::usage =
  "PolyDetLnRGaugeDa[j, mA2, ma2, dmA, radii] gives ln det \!\(\*SubsuperscriptBox[\(R\), \(j\), \((Da)\)]\) of the coupled diagonal- \
Goldstone block (3x3 for j>0, 2x2 for j=0): consecutive Bessel orders {j, j+2, j+1} (l=j-1,j+1,j), \
masses {\!\(\*SubscriptBox[\(m\), \(A\)]\), \!\(\*SubscriptBox[\(m\), \(A\)]\), \!\(\*SubscriptBox[\(m\), \(a\)]\)}, coupled by the gradient jump dmA[[s]] = g (\!\(\*SubscriptBox[\(\[Phi]\), \(s + 1\)]\)-\!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\)) at each matching \
(2nd-order derivative jump [\!\(\*SubsuperscriptBox[\(\[Psi]\), \(a\), \(\[Prime]\)]\)] = -2 dmA \!\(\*SubscriptBox[\(C\), \(ab\)]\) \!\(\*SubscriptBox[\(\[Psi]\), \(b\)]\)). mA2, ma2 are the squared shell masses \
(length nsh); dmA the signed jumps (length nsh-1); radii the matchings. det \!\(\*SuperscriptBox[\(T\), \(s\)]\) = 1. \
Options: \"Xi\" (gauge parameter, default 1 = background Feynman; \[Xi]!=1 (massive FV, D=4) uses the DIRECT \
general-\[Xi] coupled block -- transverse pair T(\!\(\*SubscriptBox[\(Z\), \(j + 1\)]\)(\!\(\*SubscriptBox[\(m\), \(A\)]\)\[Rho])) + longitudinal pair G(\!\(\*SubscriptBox[\(Z\), \(j + 1\)]\)(Sqrt[\[Xi]]\!\(\*SubscriptBox[\(m\), \(A\)]\)\[Rho])) + Goldstone at \[Xi] \!\(\*SubsuperscriptBox[\(m\), \(A\), \(2\)]\)+\!\(\*SubscriptBox[\(V\), \(aa\)]\), \
asymmetric matching jumps -- replacing the earlier \[Chi]-reduction, which was O(1)-wrong per multipole), \
\"Dimension\" (default 4; also 3, \[Xi]=1 only: orders {\[Nu]-1, \[Nu]+1, \[Nu]} with \[Nu] = j+D/2-1 half-integer, D-dimensional \
Clebsches), WorkingPrecision. dmA -> 0 decouples to the sum of three scalar ratios. A massless false \
vacuum (mA2[[-1]] == 0, symmetric FV / Coulomb phase) is supported at \[Xi] = 1: the gauge channels switch \
to the power-law false-vacuum basis while the Goldstone channel stays massive.";
PolyDetIpGauge::usage =
  "PolyDetIpGauge[mA2, ma2, radii, pmax] gives the effective UV-subtraction Laurent {\!\(\*SubsuperscriptBox[\(I\), \(1\), \(eff\)]\),..} = \
2 \!\(\*SubsuperscriptBox[\(I\), \(p\), \(mA2\)]\) + \!\(\*SubsuperscriptBox[\(I\), \(p\), \(ma2\)]\) (the D-2=2 transverse + Goldstone scalars). \!\(\*SubsuperscriptBox[\(I\), \(1\), \(eff\)]\) ~ \!\(\*SuperscriptBox[\(g\), \(2\)]\) (\[Delta] \!\(\*SubsuperscriptBox[\(m\), \(A\), \(2\)]\)), \
\!\(\*SubsuperscriptBox[\(I\), \(3\), \(eff\)]\) ~ \!\(\*SuperscriptBox[\(g\), \(4\)]\) (\[Delta] \!\(\*SubsuperscriptBox[\(m\), \(A\), \(4\)]\), the gauge-coupling running).";
PolyDetGaugeIg::usage =
  "PolyDetGaugeIg[mA2, ma2, dmA, radii] gives the gradient (wave-function) coefficient \!\(\*SubscriptBox[\(I\), \(g\)]\), the \
leading 1/(j+1) of the gauge-Goldstone mixing residual \!\(\*SuperscriptBox[\((j + 1)\), \(2\)]\) (ln det \!\(\*SuperscriptBox[\(R\), \((Da)\)]\) on - off), fitted as the \
intercept of \!\(\*SuperscriptBox[\((j + 1)\), \(3\)]\) \!\(\*SubscriptBox[\(C\), \(j\)]\) over a high-j window (option \"Window\", default {22,28}; needs WorkingPrecision >~ 6.5*28). \
Polygon-native; the continuum-robust value is the moment \!\(\*SubscriptBox[\(I\), \(g\)]\) = -\!\(\*SuperscriptBox[\(g\), \(2\)]\) Int \!\(\*SuperscriptBox[\(\[Rho]\), \(3\)]\) phibar'^2. \
Option \"Xi\" (default 1): \[Xi]!=1 returns the general-\[Xi] gradient (3-\[Xi])/2 \!\(\*SubscriptBox[\(I\), \(g\)]\)(1), the wave-function \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Phi]\)]\)(\[Xi]) = (3-\[Xi]) \!\(\*SuperscriptBox[\(g\), \(2\)]\)/\!\(\*SuperscriptBox[\((4\[Pi])\), \(2\)]\).";
PolyDetSigmaGauge::usage =
  "PolyDetSigmaGauge[mA2, ma2, dmA, radii] gives the assembled, renormalised Abelian gauge \
determinant as an Association <|Sigma, I1eff, I3eff, Ig, Nsum, GradientMethod|>. The ghost collapse (transverse \
2j(j+2) + FP ghost -2\!\(\*SuperscriptBox[\((j + 1)\), \(2\)]\) = -2) gives \!\(\*SubscriptBox[\(\[CapitalSigma]\), \(gauge\)]\) = -2 \!\(\*SubscriptBox[\(\[Sum]\), \(j\)]\) ln \!\(\*SubsuperscriptBox[\(R\), \(j\), \((T)\)]\) + \!\(\*SubscriptBox[\(\[Sum]\), \(j\)]\) \!\(\*SuperscriptBox[\((j + 1)\), \(2\)]\) lndet \
\!\(\*SubsuperscriptBox[\(R\), \(j\), \((Da)\)]\), subtracted by optimal truncation: the potential trace \!\(\*SubsuperscriptBox[\(I\), \(p\), \(eff\)]\) (\[Zeta](p-2) add-back) plus the \
gradient. By DEFAULT (Ig -> Automatic, \[Xi]=1) the gradient is the fit-free CLOSED-FORM pair sum \
\!\(\*SuperscriptBox[\((j + 1)\), \(2\)]\)\!\(\*SubscriptBox[\(C\), \(j\)]\) (the second-order matching-pair sum over the \!\(\*SubscriptBox[OverscriptBox[\(m\), \(.\)], \(A\)]\) jumps), making \[CapitalSigma] ns-convergent; \
GradientMethod -> Closed and Ig is then Missing (the running coefficient is the moment, profile-based). \
Options: \"Pmax\" (default 6), \"Nsum\" (Automatic ~ 2.4 \!\(\*SubscriptBox[\(x\), \(max\)]\)), \
WorkingPrecision, \"Ig\" (Automatic = closed-form gradient; a number = the legacy \!\(\*SubscriptBox[\(I\), \(g\)]\)/(j+1) moment subtraction, GradientMethod -> Moment), \"Xi\" (gauge parameter, \
default 1; the assembled \[Xi]!=1 \[CapitalSigma] still uses the legacy \[Chi]-route diagonal + the (3-\[Xi])/2 gradient and is \
PENDING the corrected assembly with the direct block + Nielsen validation -- treat \[Xi]!=1 \[CapitalSigma] as qualitative; the \
per-multipole PolyDetLnRGaugeDa \[Xi]!=1 is already the corrected direct chain), \
\"Dimension\" (default 4; also 3, \[Xi]=1 only: transverse+ghost collapse to MINUS ONE scalar tower, Hurwitz \
add-backs at half-integer orders, closed-form gradient PAIR SUM (diagonal jump artifact dropped, \
coherent part restored exactly) -- D=3 is finite and scale-free, GradientMethod -> ClosedPairSum). The bare \
\[CapitalSigma] is scheme/\[Mu]-dependent -- present the coefficients \!\(\*SubsuperscriptBox[\(I\), \(1\), \(eff\)]\), \!\(\*SubsuperscriptBox[\(I\), \(3\), \(eff\)]\), \!\(\*SubscriptBox[\(I\), \(g\)]\) and the running rate \
\!\(\*SubsuperscriptBox[\(I\), \(3\), \(eff\)]\) + \!\(\*SubscriptBox[\(I\), \(g\)]\). Options \"Method\" (\"Exact\" | \"FrozenWall\" | Automatic) and \"NSub\" select the case-b thin-wall fast band (\[Xi] = 1, D = 4), as for PolyDetSigma. \
A massless false vacuum (mA2[[-1]] == 0, symmetric FV) is supported at \[Xi] = 1, D = 4: the long-range \
towers are normalized against the massless free operator, the subtraction and add-back unchanged.";
PolyDetXMaxGauge::usage =
  "PolyDetXMaxGauge[mA2, radii] gives \!\(\*SubscriptBox[\(x\), \(max\)]\) = \!\(\*SubscriptBox[\(max\), \(s\)]\) \!\(\*SubscriptBox[\(m\), \(A, s\)]\) \!\(\*SubscriptBox[\(R\), \(s\)]\), the largest gauge Bessel argument.";
PolyDetInputsGauge::usage =
  "PolyDetInputsGauge[bf, g, Vaa] extracts {mA2, ma2, dmA, radii} from a BounceFunction bf for the \
U(1) Higgs gauge determinant: shell fields \!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\) (node plateaus + segment midpoints, as PolyDetInputs), \
gauge mass mA2 = \!\(\*SuperscriptBox[\(g\), \(2\)]\) \!\(\*SubsuperscriptBox[\(\[Phi]\), \(s\), \(2\)]\), would-be-Goldstone ma2 = \!\(\*SuperscriptBox[\(g\), \(2\)]\) \!\(\*SubsuperscriptBox[\(\[Phi]\), \(s\), \(2\)]\) + Vaa[\!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\)] (Vaa a pure function = \
the transverse curvature \!\(\*SubsuperscriptBox[\(\[PartialD]\), \(a\), \(2\)]\) V), gradient jumps dmA = g (\!\(\*SubscriptBox[\(\[Phi]\), \(s + 1\)]\)-\!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\)), radii from bf.";
PolyDetSigmaGaugeNonAbelian::usage =
  "PolyDetSigmaGaugeNonAbelian[bf, {{g1, n1}, {g2, n2}, ...}] gives the non-Abelian gauge determinant as a \
sum of Abelian U(1) blocks. For a bounce along a fixed direction the mass matrix \!\(\*SubsuperscriptBox[\(M\), \(A\), \(2\)]\)(\[Rho]) = C \!\(\*SuperscriptBox[\(phibar\), \(2\)]\) has a \
constant eigenbasis, so ln det factorises: each channel a is a U(1) PolyDetSigmaGauge with effective coupling \
\!\(\*SubscriptBox[\(g\), \(a\)]\) (\!\(\*SubsuperscriptBox[\(m\), \(A\), \(2\)]\) = \!\(\*SubsuperscriptBox[\(g\), \(a\), \(2\)]\) \!\(\*SuperscriptBox[\(phibar\), \(2\)]\)) and multiplicity \!\(\*SubscriptBox[\(n\), \(a\)]\); massless (\!\(\*SubscriptBox[\(g\), \(a\)]\) = 0, e.g. the photon) channels contribute 0. \
SU(2) doublet: {{g/2, 3}}; SM electroweak: {{g/2, 2}, {Sqrt[\!\(\*SuperscriptBox[\(g\), \(2\)]\)+\!\(\*SuperscriptBox[\(g\), \(\[Prime]2\)]\)]/2, 1}, {0, 1}} (W,W,Z,\[Gamma]). Returns <|Sigma, Channels|>; \
options as PolyDetSigmaGauge plus \"Vaa\" (transverse curvature pure function, default 0).";

Begin["`Private`"]

$PolyDetDir = DirectoryName[$InputFileName];   (* Kernel/ dir of this loader *)
Scan[Get[FileNameJoin[{$PolyDetDir, #}]] &,
  {"Scalar.wl", "Multiscalar.wl", "Fermion.wl", "Gauge.wl", "Rate.wl", "Syntax.wl"}];

End[]
EndPackage[]
