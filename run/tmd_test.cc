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
  LHAPDF::PDF* distpdf = LHAPDF::mkPDF("CT18NLO");

  // Heavy-quark thresholds
  std::vector<double> Thresholds;
  for (auto const& v : distpdf->flavors())
    if (v > 0 && v < 7)
      Thresholds.push_back(distpdf->quarkThreshold(v) < 0.5 ? 0 : distpdf->quarkThreshold(v));

  // Running coupling
  const auto Alphas = [&] (double const& mu) -> double{ return distpdf->alphasQ(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

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
  const auto TmdObj = apfel::InitializeTmdObjectsQuarkonium(g, Thresholds, apfel::QuarkoniumSpecies::Q_1S0_8);

  // Kinematics
  const double Vs = 13000;
  const double y  = 0;
  const double Q  = apfel::MJpsi;

  // Compute 'x1' and 'x2'
  const double x1 = Q * exp(y) / Vs;
  const double x2 = Q / exp(y) / Vs;

  // Build evolved TMD PDFs
  const auto EvTMDPDFs = BuildTmdPDFs(TmdObj, CollPDFs, Alphas, apfel::LogAccuracy::NLLp);

  // Get Drell-Yan hard-factor function
  const double hcs = 1;

  // Get Evolved TMD PDFs at bmax. This is subtracted off from the
  // b-dependent luminosity to improve the convergence of the Hankel
  // transform, particularly at small qT.
  const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDPDFs(2 * exp( - apfel::emc), Q, Q * Q).GetObjects());
  const double lumimax = xF.at(0).Evaluate(x1) * xF.at(0).Evaluate(x2);

  // Construct the TMD luminosity in b space to be trasformed in qT
  // space.
  const std::function<double(double const&)> TMDLumib = [=] (double const& b) -> double
  {
    // Compute b*
    const double bstar = std::max(b / sqrt( 1 + pow(b * distpdf->qMin() / 2 / exp( - apfel::emc), 2) ), 1e-4);

    // Get Evolved TMD PDFs and rotate them into the physical
    // basis
    const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDPDFs(bstar, Q, Q * Q).GetObjects());

    // Combine TMDs
    double lumi = xF.at(0).Evaluate(x1) * xF.at(0).Evaluate(x2);

    // Combine all pieces and return
    return b * ( lumi - lumimax );
  };

  // Double exponential quadrature
  apfel::DoubleExponentialQuadrature DEObj{};
  //apfel::OgataQuadrature DEObj{};

  // Compute predictions
  const int nqT = 20;
  const double qTmin = 0.01;
  const double qTmax = 2;
  const double qTstp = ( qTmax - qTmin ) / ( nqT - 1 );
  std::cout << std::scientific;
  std::cout << "\n     Q         "
            << "       y        "
            << "      qT        "
            << " cross section  "
            << std::endl;
  for (double qT = qTmin; qT <= qTmax; qT += qTstp)
    std::cout << Q << "\t" << y << "\t" << qT << "\t"
              << apfel::ConvFact * sigmazero(Alphas(Q), Q, "SV", apfel::QuarkoniumSpecies::Q_1S0_8) * qT * 8 * M_PI * pow(apfel::alphaem, 2) * hcs / pow(Q, 3) / 9 * DEObj.transform(TMDLumib, qT) << "\t"
              << std::endl;
  t.stop();

  return 0;
}
