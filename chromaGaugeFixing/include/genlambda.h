// -*- C++ -*-
/*! \file
 *  \brief Generate a Lambda field for Rxi gauge fixing
 */

#ifndef __genlambda_h__
#define __genlambda_h__

namespace Chroma
{

//! Generate a Lambda field for Rxi gauge fixing
/*!
 * \ingroup gfix
 *
 * Generate a Lambda field, whose Lambda^a ~ i.i.d. N(0, xi)
 * and sum_x(Lambda(x))=0. The field will be used in R_xi
 * gauge fixing iterations.
 *
 * \param lambda     Lambda field used in Rxi gauge fixing ( Modify )
 * \param xi         Rxi gauge fixing parameter ( Read )
 */

void genLambda(LatticeColorMatrix &lambda,
               const double xi);

} // end namespace Chroma

#endif
