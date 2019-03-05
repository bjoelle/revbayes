/*
 * CorrelatedMHMove.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: joellebs
 */

#include "CorrelatedMHMove.h"

#include <iostream>
#include <vector>

#include "RbException.h"
#include "DagNode.h"
#include "RbOrderedSet.h"
#include "RbMathLogic.h"
#include "RbConstants.h"

namespace RevBayesCore {

CorrelatedMHMove::CorrelatedMHMove(Proposal* main, std::vector<Proposal*> dragged,
		unsigned int ns, double w, bool autoTune) :
																				MetropolisHastingsMove(main,w,autoTune),
																				draggedProposals(dragged),
																				nSteps(ns) {
	if(draggedProposals.size() < 1) throw new RbException("Trying to initialize correlated operator without dragged proposals.");
	for(Proposal* p : draggedProposals) {
		p->setMove(this);
	}
}

CorrelatedMHMove::CorrelatedMHMove(const CorrelatedMHMove& m) :
														MetropolisHastingsMove(m),
														draggedProposals(m.draggedProposals),
														nSteps(m.nSteps) {
	for(Proposal* p : draggedProposals) {
		p->setMove(this);
	}
}

CorrelatedMHMove::~CorrelatedMHMove() {
	for(Proposal* p : draggedProposals) {
		delete p;
	}
}

Proposal& CorrelatedMHMove::getMainProposal() {
	return getProposal();
}

std::vector<Proposal*> CorrelatedMHMove::getDraggedProposals() {
	return draggedProposals;
}

unsigned int CorrelatedMHMove::getNSteps() {
	return nSteps;
}

double CorrelatedMHMove::getMoveTuningParameter() const {
	return RbConstants::Double::nan;
}

void CorrelatedMHMove::setMoveTuningParameter(double tp) {
}

void CorrelatedMHMove::performMcmcMove(double prHeat, double lHeat,
		double pHeat) {

	// copy state from before main proposal
	std::vector<DagNode> saved_nodes;
	for(DagNode* n : getDagNodes()) {
		saved_nodes.push_back(*n->clone());
	}

	getMainProposal().prepareProposal();
	double mainHR = getMainProposal().doProposal();

	// only nodes that have actually been changed by the proposal need to be saved
	for(DagNode* n : getDagNodes()) {
		std::string nm = n->getName();
		auto it = std::find_if(saved_nodes.begin(), saved_nodes.end(), [](DagNode const& obj){
			return obj.getName() == nm;
		} );
		if(it == saved_nodes.end()) {
			std::cerr << "Error, no saved matching node found for node " << nm ;
			getMainProposal().undoProposal();
			return ;
		}
		if(n->getValueAsString() == ((DagNode) *it).getValueAsString()) saved_nodes.erase(it);
	}

	if(!RbMath::isAComputableNumber(mainHR)) {
		if(RbMath::isNan(mainHR)) std::cerr << "Warning: using proposal " << getMainProposal().getProposalName() << " resulted in HastingsRatio = NaN";
	}
	else {
		double fullPosteriorRatio = computePosteriorRatio(lHeat, pHeat, prHeat);

		double loopHR;
		Proposal* p = draggedProposals[0]; //TODO update for scenario w/more than 1

		for(int ii = 0; ii < nSteps; ii++) {
			p->prepareProposal();
			loopHR = p->doProposal();

			if(!RbMath::isAComputableNumber(loopHR)) {
				if(RbMath::isNan(loopHR)) std::cerr << "Warning: using proposal " << p->getProposalName() << " resulted in HastingsRatio = NaN";
			}
			else {

			}
		}
	}

}

void CorrelatedMHMove::performHillClimbingMove(double lHeat, double pHeat) {
}

double CorrelatedMHMove::computePosteriorRatio(double lHeat, double pHeat, double prHeat) {

	double ln_likelihood_ratio = 0, ln_prior_ratio = 0;

	// Identify nodes that proposal touches
	const std::vector<DagNode*> touched_nodes = getDagNodes();
	const RbOrderedSet<DagNode*> &affected_nodes = getAffectedNodes();

	// first we touch all the nodes
	// that will set the flags for recomputation
	for (DagNode* the_node : touched_nodes) {
		// flag for recomputation
		the_node->touch();
	}

	// compute the probability of the current value for each node
	for (DagNode* the_node : touched_nodes) {

		if ( the_node->isClamped() ) ln_likelihood_ratio += the_node->getLnProbabilityRatio();
		else ln_prior_ratio += the_node->getLnProbabilityRatio();

		if ( !RbMath::isAComputableNumber(ln_prior_ratio) || !RbMath::isAComputableNumber(ln_likelihood_ratio) ) break;
	}

	// then we recompute the probability for all the affected nodes
	for (DagNode* the_node : affected_nodes) {

		if ( the_node->isClamped() ) ln_likelihood_ratio += the_node->getLnProbabilityRatio();
		else ln_prior_ratio += the_node->getLnProbabilityRatio();

		if ( !RbMath::isAComputableNumber(ln_prior_ratio) || !RbMath::isAComputableNumber(ln_likelihood_ratio) ) break;

	}

	double ln_posterior_ratio = pHeat * (lHeat * ln_likelihood_ratio + prHeat * ln_prior_ratio);
	return ln_posterior_ratio;
}

} /* namespace RevBayesCore */

