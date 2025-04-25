#ifndef FILTER_CONFIGS_CONST_H
#define FILTER_CONFIGS_CONST_H

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

double mass_piplus      = 0.13957039;
double mass_piminus     = 0.13957039;
double mass_pi0         = 0.1349768;
double mass_kplus       = 0.493677;
double mass_kminus      = 0.493677;
double mass_phi         = 1.019461;
double mass_rho0        = 0.77526;

double mass_neutron     = 0.93956542052;
double mass_proton      = 0.93827208943;
double mass_2H          = 1.875612859;
double mass_3H          = 2.808921004;
double mass_3He         = 2.809413498;
double mass_4He         = 3.727379238;
double mass_11B         = 10.25510;
double mass_11C         = 10.25708;
double mass_12C         = 11.17793;

double RadToDeg         = 180.0 / 3.14159265;

TLorentzVector boost_lorentz_vector(TLorentzVector p4, TVector3 boost_vector)
{
    TLorentzVector p4_boosted(p4);
    p4_boosted.Boost(boost_vector);
    return p4_boosted;
}

#endif // FILTER_CONFIGS_CONST_H