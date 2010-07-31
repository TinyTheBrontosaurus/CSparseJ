/* ***** BEGIN LICENSE BLOCK *****
 *
 * CSparse: a Concise Sparse matrix package.
 * Copyright (c) 2006, Timothy A. Davis.
 * http://www.cise.ufl.edu/research/sparse/CSparse
 *
 * -------------------------------------------------------------------------
 *
 * CSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
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
 * ***** END LICENSE BLOCK ***** */

package edu.emory.mathcs.csparsej.tdcomplex;

import org.apache.commons.math.complex.Complex;

import edu.emory.mathcs.csparsej.tdcomplex.Zcs_common.Zcs;

/**
 * Solve an upper triangular system L'x=b.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class Zcs_ltsolve {
    /**
     * Solves an upper triangular system L'x=b where x and b are dense. x=b on
     * input, solution on output.
     *
     * @param L
     *            column-compressed, lower triangular matrix
     * @param x
     *            size n, right hand side on input, solution on output
     * @return true if successful, false on error
     */
    public static boolean cs_ltsolve(Zcs L, Complex[] x) {
        int p, j, n, Lp[], Li[];
        Complex Lx[];
        if (!Zcs_util.CS_CSC(L) || x == null)
            return (false); /* check inputs */
        n = L.n;
        Lp = L.p;
        Li = L.i;
        Lx = L.x;
        for (j = n - 1; j >= 0; j--) {
            for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                x[j] = x[j].subtract(Lx[p].multiply(x[Li[p]]));
            }
            x[j] = x[j].divide(Lx[Lp[j]]);
        }
        return (true);
    }

}
