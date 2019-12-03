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

	Proposal& getMainProposal();
	std::vector<Proposal*> getDraggedProposals();
	unsigned int getNSteps();
	double getMoveTuningParameter() const;
	void setMoveTuningParameter(double tp);

protected:
	void performMcmcMove(double prHeat, double lHeat, double pHeat);
    void performHillClimbingMove(double lHeat, double pHeat);
    double computePosteriorRatio(double lHeat, double pHeat, double prHeat = 1);
    double computePosterior(double lHeat, double pHeat, double prHeat = 1);
    void rejectProposal(Proposal* p);
    void acceptProposal(Proposal* p);
    void restoreNodesFromSaved(bool all = true);

private:
	std::vector<Proposal*> dragged_proposals;
	std::vector<DagNode*> saved_nodes;
	std::vector<DagNode*> saved_nodes_x;
	unsigned int n_steps;
};

} /* namespace RevBayesCore */

#endif /* CORRELATEDMHMOVE_H_ */
