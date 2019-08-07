/*
 * CorrelatedMHMove.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: joellebs
 */

#include "CorrelatedMHMove.h"

#include <iostream>
#include <vector>
#include <unordered_set>

#include "RbException.h"
#include "DagNode.h"
#include "RbOrderedSet.h"
#include "RbMathLogic.h"
#include "RbConstants.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

namespace RevBayesCore {

CorrelatedMHMove::CorrelatedMHMove(Proposal* main, std::vector<Proposal*> dragged, unsigned int ns, double w, bool autoTune) :
    MetropolisHastingsMove(main,w,autoTune),
    dragged_proposals(dragged), n_steps(ns) {

    if(dragged_proposals.size() < 1) throw RbException("Trying to initialize correlated operator without dragged proposals.");

    std::unordered_set<DagNode*> dragged_nodes;
    for(Proposal* p : dragged_proposals) {
        p->setMove(this);
        for (DagNode* dn : p->getNodes()) {
            dragged_nodes.insert(dn);
        }
    }

    std::vector<DagNode*> nodes_x = getDagNodes();
    for (DagNode* dn : dragged_nodes) {
        if(std::find(nodes_x.begin(), nodes_x.end(), dn) != nodes_x.end())
            throw RbException("Correlated operator cannot have nodes in common between main and dragged proposals. "
                              "Problem found for node " + dn->getName());
        addNode(dn);
    }
}

CorrelatedMHMove::CorrelatedMHMove(const CorrelatedMHMove& m) :
    MetropolisHastingsMove(m),
    n_steps(m.n_steps) {

    dragged_proposals.clear();
    for (Proposal* p : m.dragged_proposals) {
        Proposal* np = p->clone();
        np->setMove(this);
        dragged_proposals.push_back(np);
    }
}

CorrelatedMHMove::~CorrelatedMHMove() {
    for (Proposal* p : dragged_proposals) {
        delete p;
    }
}

CorrelatedMHMove &CorrelatedMHMove::operator=(const CorrelatedMHMove &m) {
    if ( this != &m ) {
        // delegate
        MetropolisHastingsMove::operator=( m );

        dragged_proposals.clear();
        for (Proposal* p : m.dragged_proposals) {
            Proposal* np = p->clone();
            np->setMove(this);
            dragged_proposals.push_back(np);
        }
        n_steps = m.n_steps;
    }

    return *this;
}

CorrelatedMHMove *CorrelatedMHMove::clone() const {
    return new CorrelatedMHMove( *this );
}

Proposal& CorrelatedMHMove::getMainProposal() {
    return getProposal();
}

std::vector<Proposal*> CorrelatedMHMove::getDraggedProposals() {
    return dragged_proposals;
}

unsigned int CorrelatedMHMove::getNSteps() {
    return n_steps;
}

double CorrelatedMHMove::getMoveTuningParameter() const {
    return RbConstants::Double::nan;
}

void CorrelatedMHMove::setMoveTuningParameter(double) {
}

void CorrelatedMHMove::performMcmcMove(double prHeat, double lHeat,
                                       double pHeat) {

    double fullAcceptanceRatio = performMove(lHeat, pHeat, prHeat);

    if(RbMath::isAComputableNumber(fullAcceptanceRatio) &&
            (fullAcceptanceRatio >= 0 || GLOBAL_RNG -> uniform01() < exp(fullAcceptanceRatio))) {
        acceptMainProposal();
        setNumberAcceptedTotal(getNumberAcceptedTotal() + 1);
        setNumberAcceptedCurrentPeriod(getNumberAcceptedCurrentPeriod() + 1);
    }
    else {
        rejectMainProposal();
    }
}

void CorrelatedMHMove::performHillClimbingMove(double lHeat, double pHeat) {

    double fullAcceptanceRatio = performMove(lHeat, pHeat);
    if(RbMath::isAComputableNumber(fullAcceptanceRatio) && fullAcceptanceRatio >= 0) {
        acceptMainProposal();
        setNumberAcceptedTotal(getNumberAcceptedTotal() + 1);
        setNumberAcceptedCurrentPeriod(getNumberAcceptedCurrentPeriod() + 1);
    }
    else {
        rejectMainProposal();
    }
}

double CorrelatedMHMove::computePosterior(double lHeat, double pHeat, double prHeat) {

    double ln_likelihood = 0, ln_prior = 0;

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

        if ( the_node->isClamped() ) ln_likelihood += the_node->getLnProbability();
        else ln_prior += the_node->getLnProbability();

        if ( !RbMath::isAComputableNumber(ln_prior) || !RbMath::isAComputableNumber(ln_likelihood) ) break;
    }

    // then we recompute the probability for all the affected nodes
    for (DagNode* the_node : affected_nodes) {

        if ( the_node->isClamped() ) ln_likelihood += the_node->getLnProbability();
        else ln_prior += the_node->getLnProbability();

        if ( !RbMath::isAComputableNumber(ln_prior) || !RbMath::isAComputableNumber(ln_likelihood) ) break;

    }

    double ln_posterior = pHeat * (lHeat * ln_likelihood + prHeat * ln_prior);
    return ln_posterior;
}

void CorrelatedMHMove::rejectMainProposal() {
    getMainProposal().cleanProposal();
    restoreNodes();
    for(DagNode* n : getDagNodes()) {
        n->restore();
    }
    clearSaved();
}

void CorrelatedMHMove::acceptMainProposal() {
    getMainProposal().cleanProposal();
    for(DagNode* n : getDagNodes()) {
        n->keep();
    }
    clearSaved();
}

double CorrelatedMHMove::performMove(double lHeat, double pHeat,
                                     double prHeat) {
    
    saveNodes();
    double Exy = computePosterior(lHeat, pHeat, prHeat);
    
    // update to primary variable (x -> nx)
    getMainProposal().prepareProposal();
    double mainHR = getMainProposal().doProposal();
    double Enxy = computePosterior(lHeat, pHeat, prHeat);
    double fullPosteriorRatio = Enxy - Exy;
    
    if(!RbMath::isAComputableNumber(mainHR) || !RbMath::isAComputableNumber(fullPosteriorRatio)) {
        if(RbMath::isNan(mainHR)) std::cerr << "Warning: using proposal " << getMainProposal().getProposalName() << " resulted in HastingsRatio = NaN";
        if(RbMath::isNan(fullPosteriorRatio)) std::cerr << "Warning: using proposal " << getMainProposal().getProposalName() << " resulted in Posterior = NaN";
        return mainHR + fullPosteriorRatio;
    }
    
    // only nodes that have actually been changed by the proposal need to be saved
    for(DagNode* n : getDagNodes()) {
        std::string nm = n->getName();
        auto it = std::find_if(saved_nodes_x.begin(), saved_nodes_x.end(), [&](DagNode* const obj){
                return (obj -> getName() == nm);
        } );
        if(it == saved_nodes_x.end()) throw RbException("Error, no saved matching node found for node " + nm );
        if(n->getValueAsString() == (*it)->getValueAsString()) saved_nodes_x.erase(it);
    }

    double stepHR, Exny, Enxny;
    
    // updates to secondary variable(s) (y -> ny)
    for(unsigned int ii = 0; ii < n_steps; ii++) {
        
        int i = GLOBAL_RNG->uniformInt(dragged_proposals.size());
        Proposal* p = dragged_proposals[i];
        
        p->prepareProposal();
        stepHR = p->doProposal();
        
        if(!RbMath::isAComputableNumber(stepHR)) {
            if(RbMath::isNan(stepHR)) std::cerr << "Warning: using proposal " << p->getProposalName() << " resulted in HastingsRatio = NaN";
            p->undoProposal();
            p->cleanProposal();
            continue;
        }

        Enxny = computePosterior(lHeat, pHeat, prHeat);
        
        // restore value of x to calculate Exny
        switchXNodeValues();
        Exny = computePosterior(lHeat, pHeat, prHeat);
        
        // restore value of nx
        switchXNodeValues();

        double wt = static_cast<double>(ii + 1)/static_cast<double>(n_steps + 1);
        double acceptanceRatio = (Enxny - Enxy)*wt + (Exny - Exy)*(1.0 - wt) + stepHR;
        if(acceptanceRatio >= 0|| GLOBAL_RNG -> uniform01() < exp(acceptanceRatio)) {
            Exy = Exny;
            Enxy = Enxny;
        }
        else {
            p->undoProposal();
        }
        p->cleanProposal();

        fullPosteriorRatio += Enxy - Exy;
    }
    double fullAcceptanceRatio = fullPosteriorRatio/(n_steps + 1) + mainHR;
    return fullAcceptanceRatio;
}

void CorrelatedMHMove::saveNodes() {
    // copy state from before main proposal
    for(DagNode* n : getDagNodes()) {
        DagNode* copy = (n->clone());
        saved_nodes.push_back(copy);
        saved_nodes_x.push_back(copy);
    }
}

void CorrelatedMHMove::restoreNodes() {
    std::vector<DagNode*> nodes = getDagNodes();
    
    for(DagNode* n : saved_nodes) {
        std::string nm = n->getName();
        auto it = std::find_if(nodes.begin(), nodes.end(), [&](DagNode* const obj){
			return (obj -> getName() == nm);
		} );
        if(it == nodes.end()) throw RbException("Error, no matching node found for saved node " + nm );

        (*it)->copyValueFromNode(n);
    }
}

void CorrelatedMHMove::switchXNodeValues() {
    std::vector<DagNode*> nodes = getDagNodes();

    for(DagNode* n : saved_nodes_x) {
        std::string nm = n->getName();
        auto it = std::find_if(nodes.begin(), nodes.end(), [&](DagNode* const obj){
			return (obj -> getName() == nm);
		} );
        if(it == nodes.end()) throw RbException("Error, no matching node found for saved node " + nm );

        DagNode* tmp = n->clone();
        n -> copyValueFromNode(*it);
        (*it)->copyValueFromNode(tmp);
        delete tmp;
    }
}

void CorrelatedMHMove::clearSaved() {
    for(DagNode* n : saved_nodes) {
        delete n;
    }
    saved_nodes.clear();
    saved_nodes_x.clear();
}

void CorrelatedMHMove::swapNodeInternal(DagNode *oldN, DagNode *newN) {
    std::vector<DagNode*> nodes = getMainProposal().getNodes();
    if(std::find(nodes.begin(), nodes.end(), oldN) != nodes.end()) {
        getMainProposal().swapNode(oldN, newN);
        return;
    }
    for(Proposal* p : dragged_proposals) {
        nodes = p->getNodes();
        if(std::find(nodes.begin(), nodes.end(), oldN) != nodes.end()) {
            p->swapNode(oldN, newN);
        }
    }
}

void CorrelatedMHMove::tune( void ) {
    if ( num_tried_current_period <= 3 ) return;

    double rate = getNumberAcceptedCurrentPeriod() / double(num_tried_current_period);
    getMainProposal().tune( rate );
    for(Proposal* p : dragged_proposals) {
        p->tune(rate);
    }
}


} /* namespace RevBayesCore */

