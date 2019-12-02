//
// #ifndef RlSubsplit_H
// #define RlSubsplit_H
//
// #include "Clade.h"
// #include "Taxon.h"
// #include "Subsplit.h"
// #include "ModelObject.h"
//
// #include <set>
// #include <string>
// #include <vector>
//
//
// namespace RevLanguage {
//
//     /**
//     * @file
//     * This file contains the declaration of a Subsplit, which is
//     * the class that contains the names of the taxa that build the clade.
//     *
//     * @brief Declaration of Subsplit
//     *
//     * (c) Copyright 2009-
//     * @date Last modified: $Date: $
//     * @author The RevBayes Development Core Team
//     * @license GPL version 3
//     *
//     *
//     */
//     class Subsplit : public ModelObject<RevBayesCore::Subsplit> {
//
//         public:
//                                                     Subsplit(void);                                                                        //!< Constructor requires character type
//                                                     Subsplit(RevBayesCore::Subsplit *v);                                                      //!< Constructor requires character type
//                                                     Subsplit(const RevBayesCore::Subsplit &v);                                                //!< Constructor requires character type
//                                                     Subsplit(RevBayesCore::TypedDagNode<RevBayesCore::Subsplit> *n);                          //!< Constructor requires character type
//
//             typedef RevBayesCore::Subsplit             valueType;
//
//             // Basic utility functions
//             Subsplit*                                  clone(void) const;                                                                  //!< Clone object
//             void                                    constructInternalObject(void);                                                      //!< We construct the a new internal MCMC object.
//             virtual RevPtr<RevVariable>             executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f); //!< Map member methods to internal functions
//             static const std::string&               getClassType(void);                                                                 //!< Get Rev type
//             static const TypeSpec&                  getClassTypeSpec(void);                                                             //!< Get class type spec
//             std::string                             getConstructorFunctionName(void) const;                                             //!< Get the name used for the constructor function in Rev.
//             const MemberRules&                      getParameterRules(void) const;                                                      //!< Get member rules (const)
//             const TypeSpec&                         getTypeSpec(void) const;                                                            //!< Get language type of the object
//             std::string                             getGuiName(void) { return "Subsplit"; }
//             std::string                             getGuiUnicodeSymbol(void) { return "SS"; }
//             std::string                             getGuiInfo(void) { return ""; }
//
//         protected:
//             void                                    initMethods(void);
//             void                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable
//
//             std::vector<std::string>                        getHelpAuthor(void) const;                                                              //!< Get the author(s) of this function
//             std::string                                     getHelpDescription(void) const;                                                         //!< Get the description for this function
//             std::string                                     getHelpDetails(void) const;                                                             //!< Get the more detailed description of the function
//             std::string                                     getHelpExample(void) const;                                                             //!< Get an executable and instructive example
//             std::vector<RevBayesCore::RbHelpReference>      getHelpReferences(void) const;                                                          //!< Get some references/citations for this function
//             std::vector<std::string>                        getHelpSeeAlso(void) const;                                                             //!< Get suggested other functions
//             std::string                                     getHelpTitle(void) const;                                                               //!< Get the title of this help entry
//
//
//             std::vector<RevPtr<const RevVariable> >         clades;
//             RevPtr<const RevVariable>                       taxa;
//     };
//
// }
//
// #endif
