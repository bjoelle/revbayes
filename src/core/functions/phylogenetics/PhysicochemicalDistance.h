#ifndef PHYSICOCHEMICALDISTANCE_H
#define PHYSICOCHEMICALDISTANCE_H

#include <AminoAcidState.h>
#include <Cloneable.h>


namespace RevBayesCore {

class PhysicochemicalDistance : public Cloneable {

    public:
        PhysicochemicalDistance() = default;
        virtual double getDistance(AminoAcidState a1, AminoAcidState a2);
        virtual void recalculateDistances() = 0;

    protected:
        std::vector<double> distances;
};

double PhysicochemicalDistance::getDistance(AminoAcidState a1, AminoAcidState a2){
    RbBitSet b1 = a1.getState(), b2 = a2.getState();
    double d = 0;

    std::vector<size_t> set_idxs;
    for(size_t i1 = 0; i1 < b1.size(); i1++) {
        if(b1.isSet(i1)) set_idxs.push_back(i1);
    }

    for(size_t i2 = 0; i2 < b2.size(); i2++) {
        if(!b2.isSet(i2)) continue;
        for(size_t i1 : set_idxs) d += distances[i2 * 20 + i1];
    }

    d /= b1.getNumberSetBits()*b2.getNumberSetBits();

    return d;
}

}

#endif // PHYSICOCHEMICALDISTANCE_H
