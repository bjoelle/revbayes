#ifndef PiecewiseConstantFossilizedBirthDeathProcess_H
#define PiecewiseConstantFossilizedBirthDeathProcess_H

#include "RbVector.h"
#include "AbstractBirthDeathProcess.h"
#include "AbstractPiecewiseConstantFossilizedRangeProcess.h"

namespace RevBayesCore {
    
    class Clade;
    class Taxon;
    
    /**
     * @brief Piecewise-constant fossilized birth-death process with serially sampled fossils.
     *
     * The piecewise-constant birth-death process has constant rates for each time interval.
     * At the end of each time interval there may be an abrupt rate-shift (jump) for each
     * of the rates. Additionally, there may be sampling at the end of each interval.
     * Finally, fossils are sampled with rate psi, the others (fossils and extant taxa) are
     * sampled at sampling times (including the present).
     *
     * We assume that the rate vectors have one more element than the rate-change vectors.
     * Thus, one rate-change means always two interval, two rate-changes three interval, and so on.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Walker Pett)
     * @since 2014-03-18, version 1.0
     *
     */
    class PiecewiseConstantFossilizedBirthDeathProcess : public AbstractBirthDeathProcess, public AbstractPiecewiseConstantFossilizedRangeProcess {
        
        using AbstractBirthDeathProcess::taxa;

    public:
        PiecewiseConstantFossilizedBirthDeathProcess (const TypedDagNode<double>* ra,
                                                      const DagNode *speciation,
                                                      const DagNode *extinction,
                                                      const DagNode *psi,
                                                      const DagNode *counts,
                                                      const TypedDagNode<double>* rho,
                                                      const TypedDagNode<RbVector<double> > *times,
                                                      const std::string &condition,
                                                      const std::vector<Taxon> &taxa,
                                                      bool uo,
                                                      bool pa,
                                                      bool ex );  //!< Constructor
        
        // public member functions
        PiecewiseConstantFossilizedBirthDeathProcess*   clone(void) const;                                         //!< Create an independent clone

        void                                            simulateClade(std::vector<TopologyNode *> &n, double age, double present);

    protected:
        void                                            updateStartEndTimes() const;
        int                                             updateStartEndTimes(const TopologyNode & ) const;

        double                                          pSurvival(double start, double end) const;             //!< Compute the probability of survival of the process (without incomplete taxon sampling).

        // Parameter management functions
        double                                          computeLnProbabilityTimes(void) const;                            //!< Compute the log-transformed probability of the current value.
        double                                          computeLnProbabilityDivergenceTimes(void) const;                            //!< Compute the log-transformed probability of the current value.

        double                                          lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const { throw RbException("Cannot compute P(nTaxa)."); }
        double                                          lnProbTreeShape(void) const;

        double                                          simulateDivergenceTime(double origin, double present) const;    //!< Simulate a speciation event.

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

    private:
        
        // helper functions
        double                                          getMaxTaxonAge( const TopologyNode& ) const;

        mutable std::vector<bool>                       I;
        bool                                            extended;
    };
}

#endif
