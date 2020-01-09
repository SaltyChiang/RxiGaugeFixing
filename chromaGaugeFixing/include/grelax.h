// -*- C++ -*-
/*! \file
 *  \brief Perform a single gauge fixing iteration
 */

#ifndef __grelax_h__
#define __grelax_h__

namespace Chroma
{

//! Perform a single gauge fixing iteration
/*!
 * \ingroup gfix
 *
 * Performs one gauge fixing 'iteration', one checkerboard and SU(2)
 * subgroup only, for gauge fixing to Coulomb gauge in slices perpendicular
 * to the direction "j_decay".
 *
 * \param g          Current (global) gauge transformation matrices ( Modify )
 * \param u          original gauge field ( Read )
 * \param j_decay    direction perpendicular to slices to be gauge fixed ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param cb         checkerboard index ( Read )
 * \param ordo       use overrelaxation or not ( Read )
 * \param orpara     overrelaxation parameter ( Read )
 */

void grelax(LatticeColorMatrix &g,
            const multi1d<LatticeColorMatrix> &u,
            int j_decay, int su2_index, int cb, bool ordo,
            const Real &orpara);

//! Perform a single gauge fixing iteration with Lambda
/*!
 * \ingroup gfix
 *
 * Performs one gauge fixing 'iteration', one checkerboard and SU(2)
 * subgroup only, for gauge fixing to Rxi gauge in slices perpendicular
 * to the direction "j_decay".
 *
 * \param g          Current (global) gauge transformation matrices ( Modify )
 * \param u          original gauge field ( Read )
 * \param lambda     Lambda field in Rxi gauge fixing ( Read )
 * \param j_decay    direction perpendicular to slices to be gauge fixed ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param cb         checkerboard index ( Read )
 * \param ordo       use overrelaxation or not ( Read )
 * \param orpara     overrelaxation parameter ( Read )
 */

void grelax(LatticeColorMatrix &g,
            const multi1d<LatticeColorMatrix> &u,
            const LatticeColorMatrix &lambda,
            int j_decay, int su2_index, int cb, bool ordo,
            const Real &orpara);

//! Perform a single gauge fixing iteration without or with Lambda
/*!
 * \ingroup gfix
 *
 * Performs one gauge fixing 'iteration', one checkerboard and SU(2)
 * subgroup only, for gauge fixing to Coulomb or Rxi gauge in slices perpendicular
 * to the direction "j_decay".
 *
 * \param g          Current (global) gauge transformation matrices ( Modify )
 * \param u          original gauge field ( Read )
 * \param lambda     Lambda field in Rxi gauge fixing ( Read )
 * \param rxido      use Lambda field or not ( Read )
 * \param j_decay    direction perpendicular to slices to be gauge fixed ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param cb         checkerboard index ( Read )
 * \param ordo       use overrelaxation or not ( Read )
 * \param orpara     overrelaxation parameter ( Read )
 */

void grelax(LatticeColorMatrix &g,
            const multi1d<LatticeColorMatrix> &u,
            const LatticeColorMatrix &lambda, bool rxido,
            int j_decay, int su2_index, int cb, bool ordo,
            const Real &orpara);

} // end namespace Chroma

#endif
