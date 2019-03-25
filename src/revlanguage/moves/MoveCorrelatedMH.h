/*
 * MoveCorrelatedMH.h
 *
 *  Created on: Mar 8, 2019
 *      Author: joellebs
 */

#ifndef REVLANGUAGE_MOVES_MOVECORRELATEDMH_H_
#define REVLANGUAGE_MOVES_MOVECORRELATEDMH_H_

#include "RlMove.h"

namespace RevLanguage {

class Move_CorrelatedMH: public Move {

public:

	Move_CorrelatedMH();
	virtual ~Move_CorrelatedMH();

    virtual Move_CorrelatedMH* clone(void) const;  //!< Clone (deep copy)
    static const std::string& getClassType(void);  //!< Get the name of the class
    static const TypeSpec& getClassTypeSpec(void);  //!< Get the names of the class and parents
    virtual const TypeSpec& getTypeSpec(void) const;  //!< Get the names of the class and parents
	std::string getMoveName(void) const;  //!< Get the name used for the constructor function in Rev.

    void constructInternalObject(void);  //!< Construct a new internal CorrelatedMHMove.
	const MemberRules& getParameterRules(void) const;   //!< Get member rules (const)
	virtual void printValue(std::ostream& o) const;  //!< Print value (for user)

protected:

    virtual void setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable

	RevPtr<const RevVariable> main_move;
	RevPtr<const RevVariable> dragged_moves;
	RevPtr<const RevVariable> n_steps;
};

} /* namespace RevLanguage */

#endif /* REVLANGUAGE_MOVES_MOVECORRELATEDMH_H_ */
