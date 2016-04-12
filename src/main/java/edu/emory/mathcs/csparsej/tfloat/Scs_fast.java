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

/**
 * Data structure for sparse inversions
 * 
 */
public class Scs_fast
{
    public class Nonzeros {
        // Public access for inline speed only. When doing so, copy the 'quick' methods.
        public int[] mNonzeros;
        public int mNonzerosTrue = 1;
        public int[] mNexts;
        public int[] mPrevs;
        public int mHead;
        public int mTail;
        public int mLastAdd;
        public int mCount;

        public void allocate(int count) {
            mNonzeros = new int[count];
            mNexts = new int[count];
            mPrevs = new int[count];
            mCount = count;
            mHead = -1;
            mTail = -1;
            mLastAdd = -1;
        }

        public boolean isNonzero(int idx) {
            return mNonzeros[idx] == mNonzerosTrue;
        }

        /**
         * The caller's job to ensure it's not already nonzero first. If it is, infinite loop
         */
        public void addNonzeroNonfirst(int idx) {
            int next;
            int prev;

            // Find the insertion point
            if( idx > mLastAdd ) {
                next = mLastAdd;
                do {
                    // quick
                    next = mNexts[next];
                } while (next < idx);
                // May have hit the head/tail, so do the full check
                prev = getNextNonzeroDown(next);
            }
            else {
                prev = mLastAdd;
                do {
                    // quick
                    prev = mPrevs[prev];
                } while (prev > idx);
                // May have hit the head/tail, so do the full check
                next = getNextNonzeroUp(prev);
            }

            // Point at the in neighbors
            mNexts[idx] = next;
            mPrevs[idx] = prev;
            // And they should point at the new index
            // Establish a new head of the list
            if( prev < 0 ) {
                mHead = idx;
            }
            else {
                mNexts[prev] = idx;
            }
            if( next >= mCount ) {
                mTail = idx;
            }
            else {
                mPrevs[next] = idx;
            }

            mLastAdd = idx;
            // Mark it as nonzero
            mNonzeros[idx] = mNonzerosTrue;
        }

        /**
         * Wipes and adds
         * @param idx
         */
        public void addNonzeroFirstQuick(int idx) {
            // This is the only value in the list, so the next and previous point outside
            mNexts[idx] = mCount;
            mPrevs[idx] = -1;
            mHead = idx;
            mTail = idx;
            mLastAdd = idx;
            // Mark it as nonzero
            mNonzerosTrue++;
            mNonzeros[idx] = mNonzerosTrue;
        }

        /**
         * Caller did the check for out of bounds
         * @param idx
         * @return
         */
        public int getNextNonzeroUpQuick(int idx) {
            return mNexts[idx];
        }
        public int getNextNonzeroUp(int idx) {
            final int idxToCheck;
            // Must be a valid value. If the value at idx is zero, then this is a problem
            if( idx < 0 ) {
                idxToCheck = mHead;
            } else {
                idxToCheck = mNexts[idx];
            }
            return idxToCheck;
        }

        /**
         * Callers responsibilty to check bounds
         * @param idx
         * @return
         */
        public int getNextNonzeroDownQuick(int idx) {
            return mPrevs[idx];
        }

        public int getNextNonzeroDown(int idx) {
            final int idxToCheck;
            // Must be a valid value. If the value at idx is zero, then this is a problem
            if( idx >= mCount ) {
                idxToCheck = mTail;
            }
            else {
                idxToCheck = mPrevs[idx];
            }
            return idxToCheck;
        }
    }
}
