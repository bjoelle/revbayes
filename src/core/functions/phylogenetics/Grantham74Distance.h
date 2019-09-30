#ifndef GRANTHAM74DISTANCE_H
#define GRANTHAM74DISTANCE_H

#include "PhysicochemicalDistance.h"
#include <TypedDagNode.h>


namespace RevBayesCore {

class Grantham74Distance : public PhysicochemicalDistance {
    public:
        Grantham74Distance();
        Grantham74Distance(TypedDagNode<double>* alphaC, TypedDagNode<double>* alphaP, TypedDagNode<double>* alphaV);
        Grantham74Distance(TypedDagNode<double>* alphaC, TypedDagNode<double>* alphaP, TypedDagNode<double>* alphaV,
                           std::vector<double> aa_propC, std::vector<double> aa_propP, std::vector<double> aa_propV);

        virtual Grantham74Distance* clone( void ) const;

        virtual void recalculateDistances();

    private:
        TypedDagNode<double> *alphaC, *alphaP, *alphaV;
        std::vector<double> aa_propC, aa_propP, aa_propV;
};

}

#endif // GRANTHAM74DISTANCE_H
