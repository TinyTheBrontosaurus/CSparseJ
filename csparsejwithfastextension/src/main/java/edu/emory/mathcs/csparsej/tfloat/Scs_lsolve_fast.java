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
 * Solve a lower triangular system Lx=b.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class Scs_lsolve_fast
{
    /**
     * A fast version of cs_lsolve that optimizes around the fact that input/output vector X is sparse.
     * Removes error checking
     * @param L
     * @param x
     * @param nzs A boolean array indicating which values of x are nonzero. False for zero, true otherwise. This is faster than a boolean compare.
     *                 Any value in x that has mNonzeros=false can be anything, and will be assumed to be zero.
     * @return The highest index in x that is nonzero (good for input into cs_usolve_fast
     */
    public static void cs_lsolve_fast(Scs_common.Scs L, float[] x, Scs_fast.Nonzeros nzs) {
        int p;
        int j;
        final float Lx[] = L.x;
        final int   n = L.n;
        final int[] Lp = L.p;
        final int[] Li = L.i;
        int Lpj;
        int Lpj1;
        float xj;
        int Lip;

        // quicker access
        final int[] nzs_mNonzeros = nzs.mNonzeros;
        final int   nzs_mNonzerosTrue = nzs.mNonzerosTrue;
        final int[] nzs_mNexts = nzs.mNexts;

        // Step through each column (j=column index)
        j = nzs.mHead;
        while( j < n ) {
            Lpj = Lp[j];
            Lpj1 = Lp[j+1];
            // Divide the old vector/column value by the first element in the column
            x[j] = xj = x[j] /  Lx[Lpj];

            // Step through each value (row) in that column (p=column index index). Li[p] = row index
            // and subtract

            for (p = Lpj + 1; p < Lpj1; p++)
            {
                Lip = Li[p];
                // Set to nonzero using booleans.
                if( nzs_mNonzeros[Lip] == nzs_mNonzerosTrue )
                {
                    // The old value was already set to nonzero in this algorithm
                    x[Lip] += -Lx[p] * xj;
                }
                else
                {
                    // First nonzero. The old value could have been anything.
                    x[Lip] = -Lx[p] * xj;
                    nzs.addNonzeroNonfirst(Lip);
                }
            }
            j = nzs_mNexts[j];
        }
    }
}
