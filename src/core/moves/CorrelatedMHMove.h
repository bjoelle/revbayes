/*
 * CorrelatedMHMove.h
 *
 *  Created on: Feb 26, 2019
 *      Author: joellebs
 */

#ifndef CORRELATEDMHMOVE_H_
#define CORRELATEDMHMOVE_H_

#include "Proposal.h"
#include "MetropolisHastingsMove.h"

namespace RevBayesCore {

class CorrelatedMHMove: public MetropolisHastingsMove {
public:
    CorrelatedMHMove(Proposal* main, std::vector<Proposal*> dragged, unsigned int ns,
            double w, bool autoTune = false);
    CorrelatedMHMove(const CorrelatedMHMove &m);
    virtual ~CorrelatedMHMove();

    CorrelatedMHMove& operator=(const CorrelatedMHMove &m); //!< Assignment operator
    virtual CorrelatedMHMove* clone(void) const;  //!< Clone (deep copy)

    Proposal& getMainProposal();  //!< Get primary proposal
    std::vector<Proposal*> getDraggedProposals();  //!< Get dragged proposal(s)
    unsigned int getNSteps();  //!< Get number of dragged steps per move
    void tune(); //!< Tune proposals
    
    double getMoveTuningParameter() const;  //!< Get tuning parameter
    void setMoveTuningParameter(double tp);  //!< Set tuning parameter

protected:
    void performMcmcMove(double prHeat, double lHeat, double pHeat);  //!< Perform MCMC move
    void performHillClimbingMove(double lHeat, double pHeat);  //!< Perform hill-climbing move

    double performMove(double lHeat, double pHeat, double prHeat = 1);  //!< Propose correlated move
    double computePosterior(double lHeat, double pHeat, double prHeat = 1);  //!< Compute full posterior of current DAG

    void acceptMainProposal();  //!< Accept proposal and keep all nodes
    void rejectMainProposal();  //!< Reject proposal and restore all nodes
    void switchXNodeValues();  //!< Set X nodes to the saved values and vice-versa
    void saveNodes(); //!< Save all node values
    void restoreNodes(); //!< Restore all node values
    void clearSaved();  //!< Clear vectors of saved nodes

    virtual void swapNodeInternal(DagNode *oldN, DagNode *newN);  //!< Swap the pointers to the variable on which the move works on.    

private:
    std::vector<Proposal*> dragged_proposals;
    unsigned int n_steps;

    std::vector<DagNode*> saved_nodes;  //!< full save of all nodes modified by the move
    std::vector<DagNode*> saved_nodes_x;  //!< save of nodes modified by the main proposal
};

} /* namespace RevBayesCore */

#endif /* CORRELATEDMHMOVE_H_ */
