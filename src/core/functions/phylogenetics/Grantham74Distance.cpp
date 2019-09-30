#include <cmath>

#include "Grantham74Distance.h"

namespace RevBayesCore {

Grantham74Distance::Grantham74Distance() {}

Grantham74Distance::Grantham74Distance(TypedDagNode<double> *alphaC, TypedDagNode<double> *alphaP, TypedDagNode<double> *alphaV) {
    std::vector<double> default_propC { 0, 0.65, 1.33, 1.38, 2.75, 0.89, 0.92, 0.74, 0.58, 0,
                                        0, 0.33, 0, 0, 0.39, 1.42, 0.71, 0.13, 0.2, 0 };
    std::vector<double> default_propP { 8.1, 10.5, 11.6, 13, 5.5, 10.5, 12.3, 9, 10.4, 5.2,
                                        4.9, 11.3, 5.7, 5.2, 8, 9.2, 8.6, 5.4, 6.2, 5.9 };
    std::vector<double> default_propV { 31, 124, 56, 54, 55, 85, 83, 3, 96, 111, 111, 119,
                                        105, 132, 32.5, 32, 61, 170, 136, 84 };

    Grantham74Distance(alphaC, alphaP, alphaV, default_propC, default_propP, default_propV);
}

Grantham74Distance::Grantham74Distance(TypedDagNode<double> *alphaC, TypedDagNode<double> *alphaP, TypedDagNode<double> *alphaV, std::vector<double> aa_propC,
                                       std::vector<double> aa_propP, std::vector<double> aa_propV) :
    PhysicochemicalDistance(),
    alphaC(alphaC), alphaP(alphaP), alphaV(alphaV),
    aa_propC(aa_propC), aa_propP(aa_propP), aa_propV(aa_propV) {}

Grantham74Distance *Grantham74Distance::clone() const {
    return new Grantham74Distance(*this);
}

void Grantham74Distance::recalculateDistances() {
    distances = std::vector<double>(20*20, 0);

    for(unsigned long i = 0; i < 19; i++) {
        for(auto j = i+1; j < 20; j++) {
            double d = std::pow(alphaC->getValue() * std::pow(aa_propC[i] - aa_propC[j], 2) +
                                alphaP->getValue() * std::pow(aa_propP[i] - aa_propP[j], 2) +
                                alphaV->getValue() * std::pow(aa_propV[i] - aa_propV[j], 2), 0.5);
            distances[j * 20 + i] = distances[i * 20 + j] = d;
        }
    }
}

}
