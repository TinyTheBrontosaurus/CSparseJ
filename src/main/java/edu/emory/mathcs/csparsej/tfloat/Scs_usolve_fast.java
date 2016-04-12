/* ***** BEGIN LICENSE BLOCK *****
 * 
 * Fast additions to CSparse: a Concise Sparse matrix package.
 * Copyright (c) 2016, Eric Nelson
 *
 * -------------------------------------------------------------------------
 * 
 * CSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1f of the License, or (at your option) any later version.
 *
 * CSparseJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * ***** ENS LICENSE BLOCK ***** */

package edu.emory.mathcs.csparsej.tfloat;

import edu.emory.mathcs.csparsej.tfloat.Scs_common.Scs;

/**
 * Solve an upper triangular system Ux=b.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class Scs_usolve_fast
{
    /**
     * A fast version of cs_usolve that optimizes around the fact that input/output vector X is sparse.
     * Removes error checking
     * @param U
     * @param x
     * @param nzs A boolean array indicating which values of x are nonzero. False for zero, true otherwise. This is faster than a boolean compare.
     *                 Any value in x that has mNonzeros=false can be anything, and will be assumed to be zero.
     * @return
     */
    public boolean cs_usolve_fast(Scs_common.Scs U, float[] x, Scs_fast.Nonzeros nzs)
    {
        int         p;
        int         j;
        final int[] Up = U.p;
        final int[] Ui = U.i;
        final float[] Ux = U.x;
        float       xj;
        int         Uip;

        final int[] nzs_mNonzeros = nzs.mNonzeros;
        final int   nzs_mNonzerosTrue = nzs.mNonzerosTrue;
        final int[] nzs_mPrevs    = nzs.mPrevs;

        // Step through each column (j=column index)
        j = nzs.mTail;

        while (j >= 0)
        {
            final int Upjp1m1 = Up[j + 1] - 1;

            // Divide the old vector/column value by the first element in the column
            x[j] = xj = x[j] / Ux[Upjp1m1];

            // Step through each value (row) in that column (p=column index index). Li[p] = row index
            // and subtract

            for (p = Up[j]; p < Upjp1m1; p++)
            {
                Uip = Ui[p];

                // Set to nonzero using booleans.
                if (nzs_mNonzeros[Uip] == nzs_mNonzerosTrue)
                {
                    x[Uip] += -Ux[p] * xj;
                }
                else
                {
                    // First nonzero. Save it.
                    x[Uip] = -Ux[p] * xj;
                    nzs.addNonzeroNonfirst(Uip);
                }
            }
            j = nzs_mPrevs[j];
        }

        return (true);
    }
}
