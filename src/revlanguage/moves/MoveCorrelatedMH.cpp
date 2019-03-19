/*
 * MoveCorrelatedMH.cpp
 *
 *  Created on: Mar 8, 2019
 *      Author: joellebs
 */

#include "MoveCorrelatedMH.h"

#include "CorrelatedMHMove.h"
#include "MetropolisHastingsMove.h"
#include "Proposal.h"
#include "WorkspaceVector.h"
#include "RealPos.h"

namespace RevLanguage {

Move_CorrelatedMH::Move_CorrelatedMH() {}

Move_CorrelatedMH::~Move_CorrelatedMH() {}

Move_CorrelatedMH* Move_CorrelatedMH::clone(void) const {
	return new Move_CorrelatedMH(*this) ;
}

const std::string& Move_CorrelatedMH::getClassType(void) {
	static std::string rev_type = "Move_CorrelatedMH";
	return rev_type;
}

const TypeSpec& Move_CorrelatedMH::getClassTypeSpec(void) {
	static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
	return rev_type_spec;
}

const TypeSpec& Move_CorrelatedMH::getTypeSpec(void) const {
	return getClassTypeSpec();
}

std::string Move_CorrelatedMH::getMoveName(void) const {
	return "CorrelatedMH";
}

void Move_CorrelatedMH::constructInternalObject(void) {
    delete value;

	double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    unsigned int ns = static_cast<const Integer &>( n_steps->getRevObject() ).getValue();

	RevBayesCore::Move* main = &static_cast<const Move &>( main_move->getRevObject() ).getValue();
	RevBayesCore::Proposal* main_proposal;
	if(RevBayesCore::MetropolisHastingsMove* main_mh = dynamic_cast<RevBayesCore::MetropolisHastingsMove*>(main)) {
		main_proposal = (&(main_mh->getProposal()));
	}
	else throw RbException("Trying to initialize a correlated MH move with a non-MH move " + main->getMoveName());

	RevBayesCore::RbVector<Move> dragged = static_cast<const WorkspaceVector<Move> &> (dragged_moves->getRevObject()).getValue();
	std::vector<RevBayesCore::Proposal *> dragged_proposals = std::vector<RevBayesCore::Proposal *>();
	for(auto it = dragged.begin(); it != dragged.end(); ++it) {
		RevBayesCore::Move* core_move = &(it->getValue());
		if(RevBayesCore::MetropolisHastingsMove* core_mh = dynamic_cast<RevBayesCore::MetropolisHastingsMove*>(core_move)) {
            RevBayesCore::Proposal* dp = (&(core_mh->getProposal()));
            dragged_proposals.push_back(dp);
		}
		else throw RbException("Trying to initialize a correlated MH move with a non-MH move " + core_move->getMoveName());
	}

	value = new RevBayesCore::CorrelatedMHMove(main_proposal, dragged_proposals, ns, w, false);
}

const MemberRules& Move_CorrelatedMH::getParameterRules(void) const {
	static MemberRules memberRules;
	static bool rules_set = false;

	if(!rules_set) {

		memberRules.push_back( new ArgumentRule( "main", Move::getClassTypeSpec(), "The main move to use in the proposal.",
				ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
		memberRules.push_back( new ArgumentRule( "dragged", WorkspaceVector<Move>::getClassTypeSpec(),
				"The vector of dragged move(s) to use in the proposal.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
		memberRules.push_back( new ArgumentRule( "nsteps", Integer::getClassTypeSpec(), "The number of secondary steps to use in each proposal.",
						ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

		/* Inherit weight from Move, put it after variable */
		const MemberRules& inheritedRules = Move::getParameterRules();
		memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );

		rules_set = true;
	}
	return memberRules;
}

void Move_CorrelatedMH::printValue(std::ostream& o) const {
	o << "Move_CorrelatedMH(";
	if(main_move == NULL) o << "?";
	else o << main_move->getName();
	o << ", ";
	if(dragged_moves == NULL) o << "?";
	else o << dragged_moves->getName();
	o << ")";
}

void Move_CorrelatedMH::setConstParameter(const std::string& name,
		const RevPtr<const RevVariable>& var) {
	if(name == "nsteps") n_steps = var;
	else if(name == "main") main_move = var;
	else if(name == "dragged") dragged_moves = var;
	else Move::setConstParameter(name, var);
}

} /* namespace RevLanguage */
