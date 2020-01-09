// -*- C++ -*-
/*! \file
 *  \brief Rxi gauge fixing
 */

#ifndef __rxigauge_h__
#define __rxigauge_h__

namespace Chroma
{

//! Rxi gauge fixing
/*!
 * \ingroup gfix
 *
 * Driver for gauge fixing to Rxi gauge in slices perpendicular
 * to the direction "j_decay".
 * If j_decay >= Nd: fix to Landau gauge.
 * Note: as written this works only for SU(2) and SU(3)!

 * \param u        (gauge fixed) gauge field ( Modify )
 * \param lambda   Lambda field in Rxi gauge fixing ( Read )
 * \param g        Gauge transformation matrices (Write)
 * \param n_gf     number of gauge fixing iterations ( Write )
 * \param j_decay  direction perpendicular to slices to be gauge fixed ( Read )
 * \param GFAccu   desired accuracy for gauge fixing ( Read )
 * \param GFMax    maximal number of gauge fixing iterations ( Read )
 * \param OrDo     use overrelaxation or not ( Read )
 * \param OrPara   overrelaxation parameter ( Read )
 */

void rxiGauge(multi1d<LatticeColorMatrix> &u,
              LatticeColorMatrix &lambda,
              LatticeColorMatrix &g,
              int &n_gf,
              int j_decay, const Real &GFAccu, int GFMax,
              bool OrDo, const Real &OrPara);

//! Rxi gauge fixing
/*!
 * \ingroup gfix
 *
 * Driver for gauge fixing to Rxi gauge in slices perpendicular
 * to the direction "j_decay".
 * If j_decay >= Nd: fix to Landau gauge.
 * Note: as written this works only for SU(2) and SU(3)!

 * \param u        (gauge fixed) gauge field ( Modify )
 * \param lambda   Lambda field in Rxi gauge fixing ( Read )
 * \param n_gf     number of gauge fixing iterations ( Write )
 * \param j_decay  direction perpendicular to slices to be gauge fixed ( Read )
 * \param GFAccu   desired accuracy for gauge fixing ( Read )
 * \param GFMax    maximal number of gauge fixing iterations ( Read )
 * \param OrDo     use overrelaxation or not ( Read )
 * \param OrPara   overrelaxation parameter ( Read )
 */

void rxiGauge(multi1d<LatticeColorMatrix> &u,
              LatticeColorMatrix &lambda,
              int &n_gf,
              int j_decay, const Real &GFAccu, int GFMax,
              bool OrDo, const Real &OrPara);

}; // namespace Chroma

#endif
