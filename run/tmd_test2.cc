//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <LHAPDF/LHAPDF.h>
#include <quarkoniumtmd/tmdbuilderquarkonium.h>
#include <quarkoniumtmd/sigmazero.h>

// Main program
int main()
{
  // Open LHAPDF set
  LHAPDF::PDF* distpdf = LHAPDF::mkPDF("CT10nlo");

  // Heavy-quark thresholds
  std::vector<double> Thresholds;
  for (auto const& v : distpdf->flavors())
    if (v > 0 && v < 7)
      Thresholds.push_back(distpdf->quarkThreshold(v) < 0.5 ? 0 : distpdf->quarkThreshold(v));

  // Running coupling: standalone APFEL++ evolution, matching Predictions.cc
  // (instead of taking alpha_s from the LHAPDF set).
  apfel::AlphaQCD as{0.118, 91.1876, Thresholds, 1};
  const apfel::TabulateObject<double> tabas{[&] (double const& mu) -> double { return as.Evaluate(mu); }, 200, 0.49999, 2001, 3, Thresholds};
  const auto Alphas = [&] (double const& mu) -> double { return tabas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-7, 3}, apfel::SubGrid{100, 5e-5, 3}, apfel::SubGrid{60, 1e-1, 3}, apfel::SubGrid{50, 6e-1, 3}, apfel::SubGrid{50, 8e-1, 3}}};

  // Rotate PDF set into the QCD evolution basis
  const auto RotPDFs = [=] (double const& x, double const& mu) -> std::map<int,double> { return apfel::PhysToQCDEv(distpdf->xfxQ(x,mu)); };

  // Construct set of distributions as a function of the scale to be
  // tabulated
  const auto EvolvedPDFs = [=,&g] (double const& mu) -> apfel::Set<apfel::Distribution>
  {
    return apfel::Set<apfel::Distribution>{apfel::EvolutionBasisQCD{apfel::NF(mu, Thresholds)}, DistributionMap(g, RotPDFs, mu)};
  };

  // Tabulate collinear PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabPDFs{EvolvedPDFs, 200, distpdf->qMin(), distpdf->qMax(), 3, Thresholds};
  const auto CollPDFs = [&] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabPDFs.Evaluate(mu); };

  // Get timer
  apfel::Timer t;

  // Initialize TMD objects
  const auto TmdObj_approx = apfel::InitializeTmdObjectsQuarkonium(g, Thresholds, apfel::QuarkoniumSpecies::Q_approx);
  const auto TmdObj_3P21   = apfel::InitializeTmdObjectsQuarkonium(g, Thresholds, apfel::QuarkoniumSpecies::Q_3P2_1);

  // Kinematics
  const double Vs = 200.;

  // Type alias for evolved TMD PDF functions
  using EvTMDPDFs_t = std::function<apfel::Set<apfel::Distribution>(double const&, double const&, double const&)>;

  // Build evolved TMD PDFs mapped to their quarkonium species
  const std::map<apfel::QuarkoniumSpecies, EvTMDPDFs_t> EvTMDPDFsMap = {
    {apfel::QuarkoniumSpecies::Q_approx, BuildTmdPDFs(TmdObj_approx, CollPDFs, Alphas, apfel::LogAccuracy::NLLp)},
    {apfel::QuarkoniumSpecies::Q_3P2_1,  BuildTmdPDFs(TmdObj_3P21,   CollPDFs, Alphas, apfel::LogAccuracy::NLLp)},
  };

  // Species-specific masses and branching ratios
  const std::map<apfel::QuarkoniumSpecies, double> MassMap = {
    {apfel::QuarkoniumSpecies::Q_approx, apfel::MJpsi},
    {apfel::QuarkoniumSpecies::Q_3P2_1,  apfel::MChic2},
  };
  const std::map<apfel::QuarkoniumSpecies, double> BrMap = {
    {apfel::QuarkoniumSpecies::Q_approx, apfel::Br.at("Jpsi").at("mumu")},
    {apfel::QuarkoniumSpecies::Q_3P2_1,  apfel::Br.at("Chic2").at("Jpsi") * apfel::Br.at("Jpsi").at("mumu") * apfel::BQc.at("3P21")},
  };

  // Hardcoded BLNY non-perturbative Sudakov parameterisation
  // (Brock-Landry-Nadolsky-Yuan), for benchmarking against Predictions.cc.
  // f_NP(x1,x2,b,zeta) = exp( -(g1 + g2*log(Q/(2*Q0)) + g1g3*log(100*x1*x2)) * b^2/2 )
  // Parameter values are the BLNY class defaults (central values, no replica).
  const double BLNY_g1   =  0.03;
  const double BLNY_g2   =  0.87;
  const double BLNY_g1g3 = -0.17;
  const double BLNY_Q0   =  1.6;
  const auto BLNY_NPFunc = [=] (double const& x1, double const& x2, double const& b, double const& zeta) -> double
  {
    const double Q = sqrt(zeta);
    return exp( - ( BLNY_g1 + BLNY_g2 * log(Q / (2 * BLNY_Q0)) + BLNY_g1g3 * log(100 * x1 * x2) ) * b * b / 2 );
  };

  // Double exponential quadrature for the Hankel transform
  apfel::DoubleExponentialQuadrature DEObj{};

  // Rapidity bin (single bin, as in Predictions.cc's ybv = {-0.35, 0.35})
  const double ymin = -0.35;
  const double ymax =  0.35;

  // qT bins [GeV] (same edges as Predictions.cc's qTv, except qTmin here
  // is 0.01 instead of 1e-9 to avoid the qT -> 0 singularity numerically)
  const std::vector<std::pair<double,double>> qTbins = {
    {1e-5, 0.25}, {0.25, 0.50}, {0.50, 0.75}, {0.75, 1.00},
    {1.00, 1.25}, {1.25, 1.50}
  };

  // Literature reference values (as printed/compared against in
  // Predictions.cc), for the same qT bins.
  const std::vector<double> FY_results = {3.70, 3.54, 3.21, 2.76, 2.25, 1.76};

  // Print-time normalization mass, matching Predictions.cc's hardcoded
  // Mpsi = 3.0 (numerically identical to apfel::MJpsi, used here too).
  const double Mpsi_norm = 3.0;

  // LO prefactor * (branching ratio, including BQc["3P21"] for the
  // chic2 channel) for each channel -- computed once.
  std::map<apfel::QuarkoniumSpecies, double> PrefMap;
  for (auto const& [Species, EvTMDPDFs] : EvTMDPDFsMap)
  {
    const double Qs = MassMap.at(Species);
    double pref = sigmazero(Alphas(Qs), Qs, "SYY", Species) * BrMap.at(Species);
    PrefMap[Species] = pref;
  }

  // ---------------------------------------------------------------------
  // Nested double integral ("brute force" method): a genuine adaptive
  // numerical integral over BOTH y and qT of the pointwise differential
  // cross section dsigma/dy/dqT(y,qT) = sum_channels pref * qT *
  // HankelTransform(b-space luminosity at fixed y, qT). This gives the
  // cross section truly INTEGRATED over the (y,qT) bin, which is then
  // divided by the bin area to get a bin-AVERAGED quantity comparable to
  // the literature (FY_results). For a fixed y, the b-space luminosity is
  // built once per channel and reused across all qT nodes of the inner
  // adaptive integral.
  // ---------------------------------------------------------------------

  // b-space luminosity at FIXED y (includes the b-Jacobian: returns b*lumi*NPsupp)
  const auto build_lumib_atY = [&] (apfel::QuarkoniumSpecies const& Species, double const& y) -> std::function<double(double const&)>
  {
    const auto& EvTMDPDFs = EvTMDPDFsMap.at(Species);
    // TMD kinematics (x1,x2,bstar,zeta) evaluated at Qs = Mpsi for BOTH
    // channels, matching Predictions.cc's single shared luminosity scale
    // (only the LO prefactor in PrefMap uses each channel's true mass).
    // The evolution/matching itself (EvTMDPDFs) stays species-specific,
    // as already implemented.
    const double Qs   = apfel::MJpsi;
    const double zeta = Qs * Qs;
    const double bmax = 0.5; // CSS2 fixed bmax, as in Predictions.cc
    const double x1   = Qs * exp(y)  / Vs;
    const double x2   = Qs * exp(-y) / Vs;

    return [=] (double const& bT) -> double
    {
      const double bstar = std::max(bT / sqrt( 1 + pow(bT / bmax, 2) ), 1e-4);
      const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDPDFs(bstar, Qs, zeta).GetObjects());
      const double lumi   = xF.at(0).Evaluate(x1) * xF.at(0).Evaluate(x2);
      const double NPsupp = pow(BLNY_NPFunc(x1, x2, bT, zeta), 2);
      return bT * lumi * NPsupp;
    };
  };

  // Integral over a qT bin at fixed y, summed over channels: builds
  // lumib once per channel for this y and reuses it for every qT node
  // of the inner adaptive integral.
  const auto sigma_qTbin_at_y = [&] (double const& y, double const& qT_lo, double const& qT_hi) -> double
  {
    std::vector<std::pair<std::function<double(double const&)>, double>> lumibAndPref;
    lumibAndPref.reserve(EvTMDPDFsMap.size());
    for (auto const& [Species, EvTMDPDFs] : EvTMDPDFsMap)
      lumibAndPref.emplace_back(build_lumib_atY(Species, y), PrefMap.at(Species));

    const apfel::Integrator IqT{
      [&] (double const& qT) -> double
      {
        double dsigma_dy_dqT = 0;
        for (auto const& [lumib, pref] : lumibAndPref)
          dsigma_dy_dqT += pref * qT * DEObj.transform(lumib, qT);
        return dsigma_dy_dqT;
      }
    };
    return IqT.integrate(qT_lo, qT_hi, 1e-3);
  };

  // --- Compute the nested double integral, bin by bin, and compare with FY ---

  std::cout << std::scientific;
  std::cout << "\n    ymin           ymax          qTmin          qTmax   "
            << "    Resummed        Ratio to FY"
            << std::endl;

  for (std::size_t i = 0; i < qTbins.size(); ++i)
  {
    const double qT_lo  = qTbins[i].first;
    const double qT_hi  = qTbins[i].second;
    const double qT_mid = (qT_lo + qT_hi) / 2.;
    const double binArea = (ymax - ymin) * (qT_hi - qT_lo);

    // True nested double integral over (y, qT), normalized to a
    // bin-average by dividing by the bin area.
    const apfel::Integrator Iy{
      [&] (double const& y) -> double { return sigma_qTbin_at_y(y, qT_lo, qT_hi); }
    };
    const double sigma_bin_integrated = Iy.integrate(ymin, ymax, 1e-3);
    const double resumm      = apfel::ConvFact * 1e-3 * sigma_bin_integrated / binArea;
    const double resumm_norm = resumm / 2. / M_PI / qT_mid / Mpsi_norm;

    std::cout << ymin  << "\t" << ymax  << "\t"
              << qT_lo << "\t" << qT_hi << "\t"
              << resumm_norm << "\t" << resumm_norm / FY_results[i]
              << std::endl;
  }
  t.stop();

  delete distpdf;
  return 0;
}
