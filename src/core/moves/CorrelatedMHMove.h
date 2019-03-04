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

private:
	std::vector<Proposal*> draggedProposals;
	unsigned int nSteps;
};

} /* namespace RevBayesCore */

#endif /* CORRELATEDMHMOVE_H_ */
