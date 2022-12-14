(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
    typeof define === 'function' && define.amd ? define(factory) :
    (global = typeof globalThis !== 'undefined' ? globalThis : global || self, global.Constrainautor = factory());
})(this, (function () { 'use strict';

    const epsilon = 1.1102230246251565e-16;
    const splitter = 134217729;
    const resulterrbound = (3 + 8 * epsilon) * epsilon;

    // fast_expansion_sum_zeroelim routine from oritinal code
    function sum(elen, e, flen, f, h) {
        let Q, Qnew, hh, bvirt;
        let enow = e[0];
        let fnow = f[0];
        let eindex = 0;
        let findex = 0;
        if ((fnow > enow) === (fnow > -enow)) {
            Q = enow;
            enow = e[++eindex];
        } else {
            Q = fnow;
            fnow = f[++findex];
        }
        let hindex = 0;
        if (eindex < elen && findex < flen) {
            if ((fnow > enow) === (fnow > -enow)) {
                Qnew = enow + Q;
                hh = Q - (Qnew - enow);
                enow = e[++eindex];
            } else {
                Qnew = fnow + Q;
                hh = Q - (Qnew - fnow);
                fnow = f[++findex];
            }
            Q = Qnew;
            if (hh !== 0) {
                h[hindex++] = hh;
            }
            while (eindex < elen && findex < flen) {
                if ((fnow > enow) === (fnow > -enow)) {
                    Qnew = Q + enow;
                    bvirt = Qnew - Q;
                    hh = Q - (Qnew - bvirt) + (enow - bvirt);
                    enow = e[++eindex];
                } else {
                    Qnew = Q + fnow;
                    bvirt = Qnew - Q;
                    hh = Q - (Qnew - bvirt) + (fnow - bvirt);
                    fnow = f[++findex];
                }
                Q = Qnew;
                if (hh !== 0) {
                    h[hindex++] = hh;
                }
            }
        }
        while (eindex < elen) {
            Qnew = Q + enow;
            bvirt = Qnew - Q;
            hh = Q - (Qnew - bvirt) + (enow - bvirt);
            enow = e[++eindex];
            Q = Qnew;
            if (hh !== 0) {
                h[hindex++] = hh;
            }
        }
        while (findex < flen) {
            Qnew = Q + fnow;
            bvirt = Qnew - Q;
            hh = Q - (Qnew - bvirt) + (fnow - bvirt);
            fnow = f[++findex];
            Q = Qnew;
            if (hh !== 0) {
                h[hindex++] = hh;
            }
        }
        if (Q !== 0 || hindex === 0) {
            h[hindex++] = Q;
        }
        return hindex;
    }

    function sum_three(alen, a, blen, b, clen, c, tmp, out) {
        return sum(sum(alen, a, blen, b, tmp), tmp, clen, c, out);
    }

    // scale_expansion_zeroelim routine from oritinal code
    function scale(elen, e, b, h) {
        let Q, sum, hh, product1, product0;
        let bvirt, c, ahi, alo, bhi, blo;

        c = splitter * b;
        bhi = c - (c - b);
        blo = b - bhi;
        let enow = e[0];
        Q = enow * b;
        c = splitter * enow;
        ahi = c - (c - enow);
        alo = enow - ahi;
        hh = alo * blo - (Q - ahi * bhi - alo * bhi - ahi * blo);
        let hindex = 0;
        if (hh !== 0) {
            h[hindex++] = hh;
        }
        for (let i = 1; i < elen; i++) {
            enow = e[i];
            product1 = enow * b;
            c = splitter * enow;
            ahi = c - (c - enow);
            alo = enow - ahi;
            product0 = alo * blo - (product1 - ahi * bhi - alo * bhi - ahi * blo);
            sum = Q + product0;
            bvirt = sum - Q;
            hh = Q - (sum - bvirt) + (product0 - bvirt);
            if (hh !== 0) {
                h[hindex++] = hh;
            }
            Q = product1 + sum;
            hh = sum - (Q - product1);
            if (hh !== 0) {
                h[hindex++] = hh;
            }
        }
        if (Q !== 0 || hindex === 0) {
            h[hindex++] = Q;
        }
        return hindex;
    }

    function estimate(elen, e) {
        let Q = e[0];
        for (let i = 1; i < elen; i++) Q += e[i];
        return Q;
    }

    function vec(n) {
        return new Float64Array(n);
    }

    const ccwerrboundA = (3 + 16 * epsilon) * epsilon;
    const ccwerrboundB = (2 + 12 * epsilon) * epsilon;
    const ccwerrboundC = (9 + 64 * epsilon) * epsilon * epsilon;

    const B = vec(4);
    const C1 = vec(8);
    const C2 = vec(12);
    const D = vec(16);
    const u$1 = vec(4);

    function orient2dadapt(ax, ay, bx, by, cx, cy, detsum) {
        let acxtail, acytail, bcxtail, bcytail;
        let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0, s1, s0, t1, t0, u3;

        const acx = ax - cx;
        const bcx = bx - cx;
        const acy = ay - cy;
        const bcy = by - cy;

        s1 = acx * bcy;
        c = splitter * acx;
        ahi = c - (c - acx);
        alo = acx - ahi;
        c = splitter * bcy;
        bhi = c - (c - bcy);
        blo = bcy - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = acy * bcx;
        c = splitter * acy;
        ahi = c - (c - acy);
        alo = acy - ahi;
        c = splitter * bcx;
        bhi = c - (c - bcx);
        blo = bcx - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        B[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        B[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        B[2] = _j - (u3 - bvirt) + (_i - bvirt);
        B[3] = u3;

        let det = estimate(4, B);
        let errbound = ccwerrboundB * detsum;
        if (det >= errbound || -det >= errbound) {
            return det;
        }

        bvirt = ax - acx;
        acxtail = ax - (acx + bvirt) + (bvirt - cx);
        bvirt = bx - bcx;
        bcxtail = bx - (bcx + bvirt) + (bvirt - cx);
        bvirt = ay - acy;
        acytail = ay - (acy + bvirt) + (bvirt - cy);
        bvirt = by - bcy;
        bcytail = by - (bcy + bvirt) + (bvirt - cy);

        if (acxtail === 0 && acytail === 0 && bcxtail === 0 && bcytail === 0) {
            return det;
        }

        errbound = ccwerrboundC * detsum + resulterrbound * Math.abs(det);
        det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
        if (det >= errbound || -det >= errbound) return det;

        s1 = acxtail * bcy;
        c = splitter * acxtail;
        ahi = c - (c - acxtail);
        alo = acxtail - ahi;
        c = splitter * bcy;
        bhi = c - (c - bcy);
        blo = bcy - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = acytail * bcx;
        c = splitter * acytail;
        ahi = c - (c - acytail);
        alo = acytail - ahi;
        c = splitter * bcx;
        bhi = c - (c - bcx);
        blo = bcx - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        u$1[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        u$1[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        u$1[2] = _j - (u3 - bvirt) + (_i - bvirt);
        u$1[3] = u3;
        const C1len = sum(4, B, 4, u$1, C1);

        s1 = acx * bcytail;
        c = splitter * acx;
        ahi = c - (c - acx);
        alo = acx - ahi;
        c = splitter * bcytail;
        bhi = c - (c - bcytail);
        blo = bcytail - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = acy * bcxtail;
        c = splitter * acy;
        ahi = c - (c - acy);
        alo = acy - ahi;
        c = splitter * bcxtail;
        bhi = c - (c - bcxtail);
        blo = bcxtail - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        u$1[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        u$1[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        u$1[2] = _j - (u3 - bvirt) + (_i - bvirt);
        u$1[3] = u3;
        const C2len = sum(C1len, C1, 4, u$1, C2);

        s1 = acxtail * bcytail;
        c = splitter * acxtail;
        ahi = c - (c - acxtail);
        alo = acxtail - ahi;
        c = splitter * bcytail;
        bhi = c - (c - bcytail);
        blo = bcytail - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = acytail * bcxtail;
        c = splitter * acytail;
        ahi = c - (c - acytail);
        alo = acytail - ahi;
        c = splitter * bcxtail;
        bhi = c - (c - bcxtail);
        blo = bcxtail - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        u$1[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        u$1[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        u$1[2] = _j - (u3 - bvirt) + (_i - bvirt);
        u$1[3] = u3;
        const Dlen = sum(C2len, C2, 4, u$1, D);

        return D[Dlen - 1];
    }

    function orient2d(ax, ay, bx, by, cx, cy) {
        const detleft = (ay - cy) * (bx - cx);
        const detright = (ax - cx) * (by - cy);
        const det = detleft - detright;

        if (detleft === 0 || detright === 0 || (detleft > 0) !== (detright > 0)) return det;

        const detsum = Math.abs(detleft + detright);
        if (Math.abs(det) >= ccwerrboundA * detsum) return det;

        return -orient2dadapt(ax, ay, bx, by, cx, cy, detsum);
    }

    const iccerrboundA = (10 + 96 * epsilon) * epsilon;
    const iccerrboundB = (4 + 48 * epsilon) * epsilon;
    const iccerrboundC = (44 + 576 * epsilon) * epsilon * epsilon;

    const bc = vec(4);
    const ca = vec(4);
    const ab = vec(4);
    const aa = vec(4);
    const bb = vec(4);
    const cc = vec(4);
    const u = vec(4);
    const v = vec(4);
    const axtbc = vec(8);
    const aytbc = vec(8);
    const bxtca = vec(8);
    const bytca = vec(8);
    const cxtab = vec(8);
    const cytab = vec(8);
    const abt = vec(8);
    const bct = vec(8);
    const cat = vec(8);
    const abtt = vec(4);
    const bctt = vec(4);
    const catt = vec(4);

    const _8 = vec(8);
    const _16 = vec(16);
    const _16b = vec(16);
    const _16c = vec(16);
    const _32 = vec(32);
    const _32b = vec(32);
    const _48 = vec(48);
    const _64 = vec(64);

    let fin = vec(1152);
    let fin2 = vec(1152);

    function finadd(finlen, a, alen) {
        finlen = sum(finlen, fin, a, alen, fin2);
        const tmp = fin; fin = fin2; fin2 = tmp;
        return finlen;
    }

    function incircleadapt(ax, ay, bx, by, cx, cy, dx, dy, permanent) {
        let finlen;
        let adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
        let axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
        let abtlen, bctlen, catlen;
        let abttlen, bcttlen, cattlen;
        let n1, n0;

        let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0, s1, s0, t1, t0, u3;

        const adx = ax - dx;
        const bdx = bx - dx;
        const cdx = cx - dx;
        const ady = ay - dy;
        const bdy = by - dy;
        const cdy = cy - dy;

        s1 = bdx * cdy;
        c = splitter * bdx;
        ahi = c - (c - bdx);
        alo = bdx - ahi;
        c = splitter * cdy;
        bhi = c - (c - cdy);
        blo = cdy - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = cdx * bdy;
        c = splitter * cdx;
        ahi = c - (c - cdx);
        alo = cdx - ahi;
        c = splitter * bdy;
        bhi = c - (c - bdy);
        blo = bdy - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        bc[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        bc[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        bc[2] = _j - (u3 - bvirt) + (_i - bvirt);
        bc[3] = u3;
        s1 = cdx * ady;
        c = splitter * cdx;
        ahi = c - (c - cdx);
        alo = cdx - ahi;
        c = splitter * ady;
        bhi = c - (c - ady);
        blo = ady - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = adx * cdy;
        c = splitter * adx;
        ahi = c - (c - adx);
        alo = adx - ahi;
        c = splitter * cdy;
        bhi = c - (c - cdy);
        blo = cdy - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        ca[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        ca[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        ca[2] = _j - (u3 - bvirt) + (_i - bvirt);
        ca[3] = u3;
        s1 = adx * bdy;
        c = splitter * adx;
        ahi = c - (c - adx);
        alo = adx - ahi;
        c = splitter * bdy;
        bhi = c - (c - bdy);
        blo = bdy - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = bdx * ady;
        c = splitter * bdx;
        ahi = c - (c - bdx);
        alo = bdx - ahi;
        c = splitter * ady;
        bhi = c - (c - ady);
        blo = ady - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        ab[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        ab[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        ab[2] = _j - (u3 - bvirt) + (_i - bvirt);
        ab[3] = u3;

        finlen = sum(
            sum(
                sum(
                    scale(scale(4, bc, adx, _8), _8, adx, _16), _16,
                    scale(scale(4, bc, ady, _8), _8, ady, _16b), _16b, _32), _32,
                sum(
                    scale(scale(4, ca, bdx, _8), _8, bdx, _16), _16,
                    scale(scale(4, ca, bdy, _8), _8, bdy, _16b), _16b, _32b), _32b, _64), _64,
            sum(
                scale(scale(4, ab, cdx, _8), _8, cdx, _16), _16,
                scale(scale(4, ab, cdy, _8), _8, cdy, _16b), _16b, _32), _32, fin);

        let det = estimate(finlen, fin);
        let errbound = iccerrboundB * permanent;
        if (det >= errbound || -det >= errbound) {
            return det;
        }

        bvirt = ax - adx;
        adxtail = ax - (adx + bvirt) + (bvirt - dx);
        bvirt = ay - ady;
        adytail = ay - (ady + bvirt) + (bvirt - dy);
        bvirt = bx - bdx;
        bdxtail = bx - (bdx + bvirt) + (bvirt - dx);
        bvirt = by - bdy;
        bdytail = by - (bdy + bvirt) + (bvirt - dy);
        bvirt = cx - cdx;
        cdxtail = cx - (cdx + bvirt) + (bvirt - dx);
        bvirt = cy - cdy;
        cdytail = cy - (cdy + bvirt) + (bvirt - dy);
        if (adxtail === 0 && bdxtail === 0 && cdxtail === 0 && adytail === 0 && bdytail === 0 && cdytail === 0) {
            return det;
        }

        errbound = iccerrboundC * permanent + resulterrbound * Math.abs(det);
        det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) +
            2 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) +
            ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) +
            2 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) +
            ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) +
            2 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));

        if (det >= errbound || -det >= errbound) {
            return det;
        }

        if (bdxtail !== 0 || bdytail !== 0 || cdxtail !== 0 || cdytail !== 0) {
            s1 = adx * adx;
            c = splitter * adx;
            ahi = c - (c - adx);
            alo = adx - ahi;
            s0 = alo * alo - (s1 - ahi * ahi - (ahi + ahi) * alo);
            t1 = ady * ady;
            c = splitter * ady;
            ahi = c - (c - ady);
            alo = ady - ahi;
            t0 = alo * alo - (t1 - ahi * ahi - (ahi + ahi) * alo);
            _i = s0 + t0;
            bvirt = _i - s0;
            aa[0] = s0 - (_i - bvirt) + (t0 - bvirt);
            _j = s1 + _i;
            bvirt = _j - s1;
            _0 = s1 - (_j - bvirt) + (_i - bvirt);
            _i = _0 + t1;
            bvirt = _i - _0;
            aa[1] = _0 - (_i - bvirt) + (t1 - bvirt);
            u3 = _j + _i;
            bvirt = u3 - _j;
            aa[2] = _j - (u3 - bvirt) + (_i - bvirt);
            aa[3] = u3;
        }
        if (cdxtail !== 0 || cdytail !== 0 || adxtail !== 0 || adytail !== 0) {
            s1 = bdx * bdx;
            c = splitter * bdx;
            ahi = c - (c - bdx);
            alo = bdx - ahi;
            s0 = alo * alo - (s1 - ahi * ahi - (ahi + ahi) * alo);
            t1 = bdy * bdy;
            c = splitter * bdy;
            ahi = c - (c - bdy);
            alo = bdy - ahi;
            t0 = alo * alo - (t1 - ahi * ahi - (ahi + ahi) * alo);
            _i = s0 + t0;
            bvirt = _i - s0;
            bb[0] = s0 - (_i - bvirt) + (t0 - bvirt);
            _j = s1 + _i;
            bvirt = _j - s1;
            _0 = s1 - (_j - bvirt) + (_i - bvirt);
            _i = _0 + t1;
            bvirt = _i - _0;
            bb[1] = _0 - (_i - bvirt) + (t1 - bvirt);
            u3 = _j + _i;
            bvirt = u3 - _j;
            bb[2] = _j - (u3 - bvirt) + (_i - bvirt);
            bb[3] = u3;
        }
        if (adxtail !== 0 || adytail !== 0 || bdxtail !== 0 || bdytail !== 0) {
            s1 = cdx * cdx;
            c = splitter * cdx;
            ahi = c - (c - cdx);
            alo = cdx - ahi;
            s0 = alo * alo - (s1 - ahi * ahi - (ahi + ahi) * alo);
            t1 = cdy * cdy;
            c = splitter * cdy;
            ahi = c - (c - cdy);
            alo = cdy - ahi;
            t0 = alo * alo - (t1 - ahi * ahi - (ahi + ahi) * alo);
            _i = s0 + t0;
            bvirt = _i - s0;
            cc[0] = s0 - (_i - bvirt) + (t0 - bvirt);
            _j = s1 + _i;
            bvirt = _j - s1;
            _0 = s1 - (_j - bvirt) + (_i - bvirt);
            _i = _0 + t1;
            bvirt = _i - _0;
            cc[1] = _0 - (_i - bvirt) + (t1 - bvirt);
            u3 = _j + _i;
            bvirt = u3 - _j;
            cc[2] = _j - (u3 - bvirt) + (_i - bvirt);
            cc[3] = u3;
        }

        if (adxtail !== 0) {
            axtbclen = scale(4, bc, adxtail, axtbc);
            finlen = finadd(finlen, sum_three(
                scale(axtbclen, axtbc, 2 * adx, _16), _16,
                scale(scale(4, cc, adxtail, _8), _8, bdy, _16b), _16b,
                scale(scale(4, bb, adxtail, _8), _8, -cdy, _16c), _16c, _32, _48), _48);
        }
        if (adytail !== 0) {
            aytbclen = scale(4, bc, adytail, aytbc);
            finlen = finadd(finlen, sum_three(
                scale(aytbclen, aytbc, 2 * ady, _16), _16,
                scale(scale(4, bb, adytail, _8), _8, cdx, _16b), _16b,
                scale(scale(4, cc, adytail, _8), _8, -bdx, _16c), _16c, _32, _48), _48);
        }
        if (bdxtail !== 0) {
            bxtcalen = scale(4, ca, bdxtail, bxtca);
            finlen = finadd(finlen, sum_three(
                scale(bxtcalen, bxtca, 2 * bdx, _16), _16,
                scale(scale(4, aa, bdxtail, _8), _8, cdy, _16b), _16b,
                scale(scale(4, cc, bdxtail, _8), _8, -ady, _16c), _16c, _32, _48), _48);
        }
        if (bdytail !== 0) {
            bytcalen = scale(4, ca, bdytail, bytca);
            finlen = finadd(finlen, sum_three(
                scale(bytcalen, bytca, 2 * bdy, _16), _16,
                scale(scale(4, cc, bdytail, _8), _8, adx, _16b), _16b,
                scale(scale(4, aa, bdytail, _8), _8, -cdx, _16c), _16c, _32, _48), _48);
        }
        if (cdxtail !== 0) {
            cxtablen = scale(4, ab, cdxtail, cxtab);
            finlen = finadd(finlen, sum_three(
                scale(cxtablen, cxtab, 2 * cdx, _16), _16,
                scale(scale(4, bb, cdxtail, _8), _8, ady, _16b), _16b,
                scale(scale(4, aa, cdxtail, _8), _8, -bdy, _16c), _16c, _32, _48), _48);
        }
        if (cdytail !== 0) {
            cytablen = scale(4, ab, cdytail, cytab);
            finlen = finadd(finlen, sum_three(
                scale(cytablen, cytab, 2 * cdy, _16), _16,
                scale(scale(4, aa, cdytail, _8), _8, bdx, _16b), _16b,
                scale(scale(4, bb, cdytail, _8), _8, -adx, _16c), _16c, _32, _48), _48);
        }

        if (adxtail !== 0 || adytail !== 0) {
            if (bdxtail !== 0 || bdytail !== 0 || cdxtail !== 0 || cdytail !== 0) {
                s1 = bdxtail * cdy;
                c = splitter * bdxtail;
                ahi = c - (c - bdxtail);
                alo = bdxtail - ahi;
                c = splitter * cdy;
                bhi = c - (c - cdy);
                blo = cdy - bhi;
                s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
                t1 = bdx * cdytail;
                c = splitter * bdx;
                ahi = c - (c - bdx);
                alo = bdx - ahi;
                c = splitter * cdytail;
                bhi = c - (c - cdytail);
                blo = cdytail - bhi;
                t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
                _i = s0 + t0;
                bvirt = _i - s0;
                u[0] = s0 - (_i - bvirt) + (t0 - bvirt);
                _j = s1 + _i;
                bvirt = _j - s1;
                _0 = s1 - (_j - bvirt) + (_i - bvirt);
                _i = _0 + t1;
                bvirt = _i - _0;
                u[1] = _0 - (_i - bvirt) + (t1 - bvirt);
                u3 = _j + _i;
                bvirt = u3 - _j;
                u[2] = _j - (u3 - bvirt) + (_i - bvirt);
                u[3] = u3;
                s1 = cdxtail * -bdy;
                c = splitter * cdxtail;
                ahi = c - (c - cdxtail);
                alo = cdxtail - ahi;
                c = splitter * -bdy;
                bhi = c - (c - -bdy);
                blo = -bdy - bhi;
                s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
                t1 = cdx * -bdytail;
                c = splitter * cdx;
                ahi = c - (c - cdx);
                alo = cdx - ahi;
                c = splitter * -bdytail;
                bhi = c - (c - -bdytail);
                blo = -bdytail - bhi;
                t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
                _i = s0 + t0;
                bvirt = _i - s0;
                v[0] = s0 - (_i - bvirt) + (t0 - bvirt);
                _j = s1 + _i;
                bvirt = _j - s1;
                _0 = s1 - (_j - bvirt) + (_i - bvirt);
                _i = _0 + t1;
                bvirt = _i - _0;
                v[1] = _0 - (_i - bvirt) + (t1 - bvirt);
                u3 = _j + _i;
                bvirt = u3 - _j;
                v[2] = _j - (u3 - bvirt) + (_i - bvirt);
                v[3] = u3;
                bctlen = sum(4, u, 4, v, bct);
                s1 = bdxtail * cdytail;
                c = splitter * bdxtail;
                ahi = c - (c - bdxtail);
                alo = bdxtail - ahi;
                c = splitter * cdytail;
                bhi = c - (c - cdytail);
                blo = cdytail - bhi;
                s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
                t1 = cdxtail * bdytail;
                c = splitter * cdxtail;
                ahi = c - (c - cdxtail);
                alo = cdxtail - ahi;
                c = splitter * bdytail;
                bhi = c - (c - bdytail);
                blo = bdytail - bhi;
                t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
                _i = s0 - t0;
                bvirt = s0 - _i;
                bctt[0] = s0 - (_i + bvirt) + (bvirt - t0);
                _j = s1 + _i;
                bvirt = _j - s1;
                _0 = s1 - (_j - bvirt) + (_i - bvirt);
                _i = _0 - t1;
                bvirt = _0 - _i;
                bctt[1] = _0 - (_i + bvirt) + (bvirt - t1);
                u3 = _j + _i;
                bvirt = u3 - _j;
                bctt[2] = _j - (u3 - bvirt) + (_i - bvirt);
                bctt[3] = u3;
                bcttlen = 4;
            } else {
                bct[0] = 0;
                bctlen = 1;
                bctt[0] = 0;
                bcttlen = 1;
            }
            if (adxtail !== 0) {
                const len = scale(bctlen, bct, adxtail, _16c);
                finlen = finadd(finlen, sum(
                    scale(axtbclen, axtbc, adxtail, _16), _16,
                    scale(len, _16c, 2 * adx, _32), _32, _48), _48);

                const len2 = scale(bcttlen, bctt, adxtail, _8);
                finlen = finadd(finlen, sum_three(
                    scale(len2, _8, 2 * adx, _16), _16,
                    scale(len2, _8, adxtail, _16b), _16b,
                    scale(len, _16c, adxtail, _32), _32, _32b, _64), _64);

                if (bdytail !== 0) {
                    finlen = finadd(finlen, scale(scale(4, cc, adxtail, _8), _8, bdytail, _16), _16);
                }
                if (cdytail !== 0) {
                    finlen = finadd(finlen, scale(scale(4, bb, -adxtail, _8), _8, cdytail, _16), _16);
                }
            }
            if (adytail !== 0) {
                const len = scale(bctlen, bct, adytail, _16c);
                finlen = finadd(finlen, sum(
                    scale(aytbclen, aytbc, adytail, _16), _16,
                    scale(len, _16c, 2 * ady, _32), _32, _48), _48);

                const len2 = scale(bcttlen, bctt, adytail, _8);
                finlen = finadd(finlen, sum_three(
                    scale(len2, _8, 2 * ady, _16), _16,
                    scale(len2, _8, adytail, _16b), _16b,
                    scale(len, _16c, adytail, _32), _32, _32b, _64), _64);
            }
        }
        if (bdxtail !== 0 || bdytail !== 0) {
            if (cdxtail !== 0 || cdytail !== 0 || adxtail !== 0 || adytail !== 0) {
                s1 = cdxtail * ady;
                c = splitter * cdxtail;
                ahi = c - (c - cdxtail);
                alo = cdxtail - ahi;
                c = splitter * ady;
                bhi = c - (c - ady);
                blo = ady - bhi;
                s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
                t1 = cdx * adytail;
                c = splitter * cdx;
                ahi = c - (c - cdx);
                alo = cdx - ahi;
                c = splitter * adytail;
                bhi = c - (c - adytail);
                blo = adytail - bhi;
                t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
                _i = s0 + t0;
                bvirt = _i - s0;
                u[0] = s0 - (_i - bvirt) + (t0 - bvirt);
                _j = s1 + _i;
                bvirt = _j - s1;
                _0 = s1 - (_j - bvirt) + (_i - bvirt);
                _i = _0 + t1;
                bvirt = _i - _0;
                u[1] = _0 - (_i - bvirt) + (t1 - bvirt);
                u3 = _j + _i;
                bvirt = u3 - _j;
                u[2] = _j - (u3 - bvirt) + (_i - bvirt);
                u[3] = u3;
                n1 = -cdy;
                n0 = -cdytail;
                s1 = adxtail * n1;
                c = splitter * adxtail;
                ahi = c - (c - adxtail);
                alo = adxtail - ahi;
                c = splitter * n1;
                bhi = c - (c - n1);
                blo = n1 - bhi;
                s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
                t1 = adx * n0;
                c = splitter * adx;
                ahi = c - (c - adx);
                alo = adx - ahi;
                c = splitter * n0;
                bhi = c - (c - n0);
                blo = n0 - bhi;
                t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
                _i = s0 + t0;
                bvirt = _i - s0;
                v[0] = s0 - (_i - bvirt) + (t0 - bvirt);
                _j = s1 + _i;
                bvirt = _j - s1;
                _0 = s1 - (_j - bvirt) + (_i - bvirt);
                _i = _0 + t1;
                bvirt = _i - _0;
                v[1] = _0 - (_i - bvirt) + (t1 - bvirt);
                u3 = _j + _i;
                bvirt = u3 - _j;
                v[2] = _j - (u3 - bvirt) + (_i - bvirt);
                v[3] = u3;
                catlen = sum(4, u, 4, v, cat);
                s1 = cdxtail * adytail;
                c = splitter * cdxtail;
                ahi = c - (c - cdxtail);
                alo = cdxtail - ahi;
                c = splitter * adytail;
                bhi = c - (c - adytail);
                blo = adytail - bhi;
                s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
                t1 = adxtail * cdytail;
                c = splitter * adxtail;
                ahi = c - (c - adxtail);
                alo = adxtail - ahi;
                c = splitter * cdytail;
                bhi = c - (c - cdytail);
                blo = cdytail - bhi;
                t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
                _i = s0 - t0;
                bvirt = s0 - _i;
                catt[0] = s0 - (_i + bvirt) + (bvirt - t0);
                _j = s1 + _i;
                bvirt = _j - s1;
                _0 = s1 - (_j - bvirt) + (_i - bvirt);
                _i = _0 - t1;
                bvirt = _0 - _i;
                catt[1] = _0 - (_i + bvirt) + (bvirt - t1);
                u3 = _j + _i;
                bvirt = u3 - _j;
                catt[2] = _j - (u3 - bvirt) + (_i - bvirt);
                catt[3] = u3;
                cattlen = 4;
            } else {
                cat[0] = 0;
                catlen = 1;
                catt[0] = 0;
                cattlen = 1;
            }
            if (bdxtail !== 0) {
                const len = scale(catlen, cat, bdxtail, _16c);
                finlen = finadd(finlen, sum(
                    scale(bxtcalen, bxtca, bdxtail, _16), _16,
                    scale(len, _16c, 2 * bdx, _32), _32, _48), _48);

                const len2 = scale(cattlen, catt, bdxtail, _8);
                finlen = finadd(finlen, sum_three(
                    scale(len2, _8, 2 * bdx, _16), _16,
                    scale(len2, _8, bdxtail, _16b), _16b,
                    scale(len, _16c, bdxtail, _32), _32, _32b, _64), _64);

                if (cdytail !== 0) {
                    finlen = finadd(finlen, scale(scale(4, aa, bdxtail, _8), _8, cdytail, _16), _16);
                }
                if (adytail !== 0) {
                    finlen = finadd(finlen, scale(scale(4, cc, -bdxtail, _8), _8, adytail, _16), _16);
                }
            }
            if (bdytail !== 0) {
                const len = scale(catlen, cat, bdytail, _16c);
                finlen = finadd(finlen, sum(
                    scale(bytcalen, bytca, bdytail, _16), _16,
                    scale(len, _16c, 2 * bdy, _32), _32, _48), _48);

                const len2 = scale(cattlen, catt, bdytail, _8);
                finlen = finadd(finlen, sum_three(
                    scale(len2, _8, 2 * bdy, _16), _16,
                    scale(len2, _8, bdytail, _16b), _16b,
                    scale(len, _16c, bdytail, _32), _32,  _32b, _64), _64);
            }
        }
        if (cdxtail !== 0 || cdytail !== 0) {
            if (adxtail !== 0 || adytail !== 0 || bdxtail !== 0 || bdytail !== 0) {
                s1 = adxtail * bdy;
                c = splitter * adxtail;
                ahi = c - (c - adxtail);
                alo = adxtail - ahi;
                c = splitter * bdy;
                bhi = c - (c - bdy);
                blo = bdy - bhi;
                s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
                t1 = adx * bdytail;
                c = splitter * adx;
                ahi = c - (c - adx);
                alo = adx - ahi;
                c = splitter * bdytail;
                bhi = c - (c - bdytail);
                blo = bdytail - bhi;
                t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
                _i = s0 + t0;
                bvirt = _i - s0;
                u[0] = s0 - (_i - bvirt) + (t0 - bvirt);
                _j = s1 + _i;
                bvirt = _j - s1;
                _0 = s1 - (_j - bvirt) + (_i - bvirt);
                _i = _0 + t1;
                bvirt = _i - _0;
                u[1] = _0 - (_i - bvirt) + (t1 - bvirt);
                u3 = _j + _i;
                bvirt = u3 - _j;
                u[2] = _j - (u3 - bvirt) + (_i - bvirt);
                u[3] = u3;
                n1 = -ady;
                n0 = -adytail;
                s1 = bdxtail * n1;
                c = splitter * bdxtail;
                ahi = c - (c - bdxtail);
                alo = bdxtail - ahi;
                c = splitter * n1;
                bhi = c - (c - n1);
                blo = n1 - bhi;
                s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
                t1 = bdx * n0;
                c = splitter * bdx;
                ahi = c - (c - bdx);
                alo = bdx - ahi;
                c = splitter * n0;
                bhi = c - (c - n0);
                blo = n0 - bhi;
                t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
                _i = s0 + t0;
                bvirt = _i - s0;
                v[0] = s0 - (_i - bvirt) + (t0 - bvirt);
                _j = s1 + _i;
                bvirt = _j - s1;
                _0 = s1 - (_j - bvirt) + (_i - bvirt);
                _i = _0 + t1;
                bvirt = _i - _0;
                v[1] = _0 - (_i - bvirt) + (t1 - bvirt);
                u3 = _j + _i;
                bvirt = u3 - _j;
                v[2] = _j - (u3 - bvirt) + (_i - bvirt);
                v[3] = u3;
                abtlen = sum(4, u, 4, v, abt);
                s1 = adxtail * bdytail;
                c = splitter * adxtail;
                ahi = c - (c - adxtail);
                alo = adxtail - ahi;
                c = splitter * bdytail;
                bhi = c - (c - bdytail);
                blo = bdytail - bhi;
                s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
                t1 = bdxtail * adytail;
                c = splitter * bdxtail;
                ahi = c - (c - bdxtail);
                alo = bdxtail - ahi;
                c = splitter * adytail;
                bhi = c - (c - adytail);
                blo = adytail - bhi;
                t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
                _i = s0 - t0;
                bvirt = s0 - _i;
                abtt[0] = s0 - (_i + bvirt) + (bvirt - t0);
                _j = s1 + _i;
                bvirt = _j - s1;
                _0 = s1 - (_j - bvirt) + (_i - bvirt);
                _i = _0 - t1;
                bvirt = _0 - _i;
                abtt[1] = _0 - (_i + bvirt) + (bvirt - t1);
                u3 = _j + _i;
                bvirt = u3 - _j;
                abtt[2] = _j - (u3 - bvirt) + (_i - bvirt);
                abtt[3] = u3;
                abttlen = 4;
            } else {
                abt[0] = 0;
                abtlen = 1;
                abtt[0] = 0;
                abttlen = 1;
            }
            if (cdxtail !== 0) {
                const len = scale(abtlen, abt, cdxtail, _16c);
                finlen = finadd(finlen, sum(
                    scale(cxtablen, cxtab, cdxtail, _16), _16,
                    scale(len, _16c, 2 * cdx, _32), _32, _48), _48);

                const len2 = scale(abttlen, abtt, cdxtail, _8);
                finlen = finadd(finlen, sum_three(
                    scale(len2, _8, 2 * cdx, _16), _16,
                    scale(len2, _8, cdxtail, _16b), _16b,
                    scale(len, _16c, cdxtail, _32), _32, _32b, _64), _64);

                if (adytail !== 0) {
                    finlen = finadd(finlen, scale(scale(4, bb, cdxtail, _8), _8, adytail, _16), _16);
                }
                if (bdytail !== 0) {
                    finlen = finadd(finlen, scale(scale(4, aa, -cdxtail, _8), _8, bdytail, _16), _16);
                }
            }
            if (cdytail !== 0) {
                const len = scale(abtlen, abt, cdytail, _16c);
                finlen = finadd(finlen, sum(
                    scale(cytablen, cytab, cdytail, _16), _16,
                    scale(len, _16c, 2 * cdy, _32), _32, _48), _48);

                const len2 = scale(abttlen, abtt, cdytail, _8);
                finlen = finadd(finlen, sum_three(
                    scale(len2, _8, 2 * cdy, _16), _16,
                    scale(len2, _8, cdytail, _16b), _16b,
                    scale(len, _16c, cdytail, _32), _32, _32b, _64), _64);
            }
        }

        return fin[finlen - 1];
    }

    function incircle(ax, ay, bx, by, cx, cy, dx, dy) {
        const adx = ax - dx;
        const bdx = bx - dx;
        const cdx = cx - dx;
        const ady = ay - dy;
        const bdy = by - dy;
        const cdy = cy - dy;

        const bdxcdy = bdx * cdy;
        const cdxbdy = cdx * bdy;
        const alift = adx * adx + ady * ady;

        const cdxady = cdx * ady;
        const adxcdy = adx * cdy;
        const blift = bdx * bdx + bdy * bdy;

        const adxbdy = adx * bdy;
        const bdxady = bdx * ady;
        const clift = cdx * cdx + cdy * cdy;

        const det =
            alift * (bdxcdy - cdxbdy) +
            blift * (cdxady - adxcdy) +
            clift * (adxbdy - bdxady);

        const permanent =
            (Math.abs(bdxcdy) + Math.abs(cdxbdy)) * alift +
            (Math.abs(cdxady) + Math.abs(adxcdy)) * blift +
            (Math.abs(adxbdy) + Math.abs(bdxady)) * clift;

        const errbound = iccerrboundA * permanent;

        if (det > errbound || -det > errbound) {
            return det;
        }
        return incircleadapt(ax, ay, bx, by, cx, cy, dx, dy, permanent);
    }

    /**
     * A set of numbers, stored as bits in a typed array. The amount of numbers /
     * the maximum number that can be stored is limited by the length, which is
     * fixed at construction time.
     */
    class BitSet {
        constructor(W, bs) {
            this.W = W;
            this.bs = bs;
        }
        /**
         * Add a number to the set.
         *
         * @param idx The number to add. Must be 0 <= idx < len.
         * @return this.
         */
        add(idx) {
            const W = this.W, byte = (idx / W) | 0, bit = idx % W;
            this.bs[byte] |= 1 << bit;
            return this;
        }
        /**
         * Delete a number from the set.
         *
         * @param idx The number to delete. Must be 0 <= idx < len.
         * @return this.
         */
        delete(idx) {
            const W = this.W, byte = (idx / W) | 0, bit = idx % W;
            this.bs[byte] &= ~(1 << bit);
            return this;
        }
        /**
         * Add or delete a number in the set, depending on the second argument.
         *
         * @param idx The number to add or delete. Must be 0 <= idx < len.
         * @param val If true, add the number, otherwise delete.
         * @return val.
         */
        set(idx, val) {
            const W = this.W, byte = (idx / W) | 0, bit = idx % W, m = 1 << bit;
            //this.bs[byte] = set ? this.bs[byte] | m : this.bs[byte] & ~m;
            this.bs[byte] ^= (-val ^ this.bs[byte]) & m; // -set == set * 255
            return val;
        }
        /**
         * Whether the number is in the set.
         *
         * @param idx The number to test. Must be 0 <= idx < len.
         * @return True if the number is in the set.
         */
        has(idx) {
            const W = this.W, byte = (idx / W) | 0, bit = idx % W;
            return !!(this.bs[byte] & (1 << bit));
        }
        /**
         * Iterate over the numbers that are in the set. The callback is invoked
         * with each number that is set. It is allowed to change the BitSet during
         * iteration. If it deletes a number that has not been iterated over, that
         * number will not show up in a later call. If it adds a number during
         * iteration, that number may or may not show up in a later call.
         *
         * @param fn The function to call for each number.
         * @return this.
         */
        forEach(fn) {
            const W = this.W, bs = this.bs, len = bs.length;
            for (let byte = 0; byte < len; byte++) {
                let bit = 0;
                // bs[byte] may change during iteration
                while (bs[byte] && bit < W) {
                    if (bs[byte] & (1 << bit)) {
                        fn(byte * W + bit);
                    }
                    bit++;
                }
            }
            return this;
        }
    }
    /**
     * A bit set using 8 bits per cell.
     */
    class BitSet8 extends BitSet {
        /**
         * Create a bit set.
         *
         * @param len The length of the bit set, limiting the maximum value that
         *        can be stored in it to len - 1.
         */
        constructor(len) {
            const W = 8, bs = new Uint8Array(Math.ceil(len / W)).fill(0);
            super(W, bs);
        }
    }

    function nextEdge(e) { return (e % 3 === 2) ? e - 2 : e + 1; }
    function prevEdge(e) { return (e % 3 === 0) ? e + 2 : e - 1; }
    /**
     * Constrain a triangulation from Delaunator, using (parts of) the algorithm
     * in "A fast algorithm for generating constrained Delaunay triangulations" by
     * S. W. Sloan.
     */
    class Constrainautor {
        /**
         * Make a Constrainautor.
         *
         * @param del The triangulation output from Delaunator.
         * @param edges If provided, constrain these edges as by constrainAll.
         */
        constructor(del, edges) {
            if (!del || typeof del !== 'object' || !del.triangles || !del.halfedges || !del.coords) {
                throw new Error("Expected an object with Delaunator output");
            }
            if (del.triangles.length % 3 || del.halfedges.length !== del.triangles.length || del.coords.length % 2) {
                throw new Error("Delaunator output appears inconsistent");
            }
            if (del.triangles.length < 3) {
                throw new Error("No edges in triangulation");
            }
            this.del = del;
            const U32NIL = 2 ** 32 - 1, // Max value of a Uint32Array: use as a sentinel for not yet defined 
            numPoints = del.coords.length >> 1, numEdges = del.triangles.length;
            // Map every vertex id to the right-most edge that points to that vertex.
            this.vertMap = new Uint32Array(numPoints).fill(U32NIL);
            // Keep track of edges flipped while constraining
            this.flips = new BitSet8(numEdges);
            // Keep track of constrained edges
            this.consd = new BitSet8(numEdges);
            for (let e = 0; e < numEdges; e++) {
                const v = del.triangles[e];
                if (this.vertMap[v] === U32NIL) {
                    this.updateVert(e);
                }
            }
            if (edges) {
                this.constrainAll(edges);
            }
        }
        /**
         * Constrain the triangulation such that there is an edge between p1 and p2.
         *
         * @param segP1 The index of one segment end-point in the coords array.
         * @param segP2 The index of the other segment end-point in the coords array.
         * @return The id of the edge that points from p1 to p2. If the
         *         constrained edge lies on the hull and points in the opposite
         *         direction (p2 to p1), the negative of its id is returned.
         */
        constrainOne(segP1, segP2) {
            const { triangles, halfedges } = this.del, vm = this.vertMap, consd = this.consd, start = vm[segP1];
            // Loop over the edges touching segP1
            let edg = start;
            do {
                // edg points toward segP1, so its start-point is opposite it
                const p4 = triangles[edg], nxt = nextEdge(edg);
                // already constrained, but in reverse order
                if (p4 === segP2) {
                    return this.protect(edg);
                }
                // The edge opposite segP1
                const opp = prevEdge(edg), p3 = triangles[opp];
                // already constrained
                if (p3 === segP2) {
                    this.protect(nxt);
                    return nxt;
                }
                // edge opposite segP1 intersects constraint
                if (this.intersectSegments(segP1, segP2, p3, p4)) {
                    edg = opp;
                    break;
                }
                const adj = halfedges[nxt];
                // The next edge pointing to segP1
                edg = adj;
            } while (edg !== -1 && edg !== start);
            let conEdge = edg;
            // Walk through the triangulation looking for further intersecting
            // edges and flip them. If an intersecting edge cannot be flipped,
            // assign its id to `rescan` and restart from there, until there are
            // no more intersects.
            let rescan = -1;
            while (edg !== -1) {
                // edg is the intersecting half-edge in the triangle we came from
                // adj is now the opposite half-edge in the adjacent triangle, which
                // is away from segP1.
                const adj = halfedges[edg], 
                // cross diagonal
                bot = prevEdge(edg), top = prevEdge(adj), rgt = nextEdge(adj);
                if (adj === -1) {
                    throw new Error("Constraining edge exited the hull");
                }
                if (consd.has(edg)) { // || consd.has(adj) // assume consd is consistent
                    throw new Error("Edge intersects already constrained edge");
                }
                if (this.isCollinear(segP1, segP2, triangles[edg]) ||
                    this.isCollinear(segP1, segP2, triangles[adj])) {
                    throw new Error("Constraining edge intersects point");
                }
                const convex = this.intersectSegments(triangles[edg], triangles[adj], triangles[bot], triangles[top]);
                // The quadrilateral formed by the two triangles adjoing edg is not
                // convex, so the edge can't be flipped. Continue looking for the
                // next intersecting edge and restart at this one later.
                if (!convex) {
                    if (rescan === -1) {
                        rescan = edg;
                    }
                    if (triangles[top] === segP2) {
                        if (edg === rescan) {
                            throw new Error("Infinite loop: non-convex quadrilateral");
                        }
                        edg = rescan;
                        rescan = -1;
                        continue;
                    }
                    // Look for the next intersect
                    if (this.intersectSegments(segP1, segP2, triangles[top], triangles[adj])) {
                        edg = top;
                    }
                    else if (this.intersectSegments(segP1, segP2, triangles[rgt], triangles[top])) {
                        edg = rgt;
                    }
                    else if (rescan === edg) {
                        throw new Error("Infinite loop: no further intersect after non-convex");
                    }
                    continue;
                }
                this.flipDiagonal(edg);
                // The new edge might still intersect, which will be fixed in the
                // next rescan.
                if (this.intersectSegments(segP1, segP2, triangles[bot], triangles[top])) {
                    if (rescan === -1) {
                        rescan = bot;
                    }
                    if (rescan === bot) {
                        throw new Error("Infinite loop: flipped diagonal still intersects");
                    }
                }
                // Reached the other segment end-point? Start the rescan.
                if (triangles[top] === segP2) {
                    conEdge = top;
                    edg = rescan;
                    rescan = -1;
                    // Otherwise, for the next edge that intersects. Because we just
                    // flipped, it's either edg again, or rgt.
                }
                else if (this.intersectSegments(segP1, segP2, triangles[rgt], triangles[top])) {
                    edg = rgt;
                }
            }
            const flips = this.flips;
            this.protect(conEdge);
            do {
                // need to use var to scope it outside the loop, but re-initialize
                // to 0 each iteration
                var flipped = 0;
                flips.forEach(edg => {
                    flips.delete(edg);
                    const adj = halfedges[edg];
                    if (adj === -1) {
                        return;
                    }
                    flips.delete(adj);
                    if (!this.isDelaunay(edg)) {
                        this.flipDiagonal(edg);
                        flipped++;
                    }
                });
            } while (flipped > 0);
            return this.findEdge(segP1, segP2);
        }
        /**
         * Fix the Delaunay condition. It is no longer necessary to call this
         * method after constraining (many) edges, since constrainOne will do it
         * after each.
         *
         * @param deep If true, keep checking & flipping edges until all
         *        edges are Delaunay, otherwise only check the edges once.
         * @return The triangulation object.
         */
        delaunify(deep = false) {
            const halfedges = this.del.halfedges, flips = this.flips, consd = this.consd, len = halfedges.length;
            do {
                var flipped = 0;
                for (let edg = 0; edg < len; edg++) {
                    if (consd.has(edg)) {
                        continue;
                    }
                    flips.delete(edg);
                    const adj = halfedges[edg];
                    if (adj === -1) {
                        continue;
                    }
                    flips.delete(adj);
                    if (!this.isDelaunay(edg)) {
                        this.flipDiagonal(edg);
                        flipped++;
                    }
                }
            } while (deep && flipped > 0);
            return this;
        }
        /**
         * Call constrainOne on each edge, and delaunify afterwards.
         *
         * @param edges The edges to constrain: each element is an array with
         *        [p1, p2] which are indices into the points array originally
         *        supplied to Delaunator.
         * @return The triangulation object.
         */
        constrainAll(edges) {
            const len = edges.length;
            for (let i = 0; i < len; i++) {
                const e = edges[i];
                this.constrainOne(e[0], e[1]);
            }
            return this;
        }
        /**
         * Whether an edge is a constrained edge.
         *
         * @param edg The edge id.
         * @return True if the edge is constrained.
         */
        isConstrained(edg) {
            return this.consd.has(edg);
        }
        /**
         * Find the edge that points from p1 -> p2. If there is only an edge from
         * p2 -> p1 (i.e. it is on the hull), returns the negative id of it.
         *
         * @param p1 The index of the first point into the points array.
         * @param p2 The index of the second point into the points array.
         * @return The id of the edge that points from p1 -> p2, or the negative
         *         id of the edge that goes from p2 -> p1, or Infinity if there is
         *         no edge between p1 and p2.
         */
        findEdge(p1, p2) {
            const start1 = this.vertMap[p2], { triangles, halfedges } = this.del;
            let edg = start1, prv = -1;
            // Walk around p2, iterating over the edges pointing to it
            do {
                if (triangles[edg] === p1) {
                    return edg;
                }
                prv = nextEdge(edg);
                edg = halfedges[prv];
            } while (edg !== -1 && edg !== start1);
            // Did not find p1 -> p2, the only option is that it is on the hull on
            // the 'left-hand' side, pointing p2 -> p1 (or there is no edge)
            if (triangles[nextEdge(prv)] === p1) {
                return -prv;
            }
            return Infinity;
        }
        /**
         * Mark an edge as constrained, i.e. should not be touched by `delaunify`.
         *
         * @private
         * @param edg The edge id.
         * @return If edg has an adjacent, returns that, otherwise -edg.
         */
        protect(edg) {
            const adj = this.del.halfedges[edg], flips = this.flips, consd = this.consd;
            flips.delete(edg);
            consd.add(edg);
            if (adj !== -1) {
                flips.delete(adj);
                consd.add(adj);
                return adj;
            }
            return -edg;
        }
        /**
         * Mark an edge as flipped, unless it is already marked as constrained.
         *
         * @private
         * @param edg The edge id.
         * @return True if edg was not constrained.
         */
        markFlip(edg) {
            const halfedges = this.del.halfedges, flips = this.flips, consd = this.consd;
            if (consd.has(edg)) {
                return false;
            }
            const adj = halfedges[edg];
            if (adj !== -1) {
                flips.add(edg);
                flips.add(adj);
            }
            return true;
        }
        /**
         * Flip the edge shared by two triangles.
         *
         * @private
         * @param edg The edge shared by the two triangles, must have an
         *        adjacent half-edge.
         * @return The new diagonal.
         */
        flipDiagonal(edg) {
            // Flip a diagonal
            //                top                     edg
            //          o  <----- o            o <------  o 
            //         | ^ \      ^           |       ^ / ^
            //      lft|  \ \     |        lft|      / /  |
            //         |   \ \adj |           |  bot/ /   |
            //         | edg\ \   |           |    / /top |
            //         |     \ \  |rgt        |   / /     |rgt
            //         v      \ v |           v  / v      |
            //         o ----->  o            o   ------> o 
            //           bot                     adj
            const { triangles, halfedges } = this.del, flips = this.flips, consd = this.consd, adj = halfedges[edg], bot = prevEdge(edg), lft = nextEdge(edg), top = prevEdge(adj), rgt = nextEdge(adj), adjBot = halfedges[bot], adjTop = halfedges[top];
            if (consd.has(edg)) { // || consd.has(adj) // assume consd is consistent
                throw new Error("Trying to flip a constrained edge");
            }
            // move *edg to *top
            triangles[edg] = triangles[top];
            halfedges[edg] = adjTop;
            if (!flips.set(edg, flips.has(top))) {
                consd.set(edg, consd.has(top));
            }
            if (adjTop !== -1) {
                halfedges[adjTop] = edg;
            }
            halfedges[bot] = top;
            // move *adj to *bot
            triangles[adj] = triangles[bot];
            halfedges[adj] = adjBot;
            if (!flips.set(adj, flips.has(bot))) {
                consd.set(adj, consd.has(bot));
            }
            if (adjBot !== -1) {
                halfedges[adjBot] = adj;
            }
            halfedges[top] = bot;
            this.markFlip(edg);
            this.markFlip(lft);
            this.markFlip(adj);
            this.markFlip(rgt);
            // mark flips unconditionally
            flips.add(bot);
            consd.delete(bot);
            flips.add(top);
            consd.delete(top);
            this.updateVert(edg);
            this.updateVert(lft);
            this.updateVert(adj);
            this.updateVert(rgt);
            return bot;
        }
        /**
         * Whether the two triangles sharing edg conform to the Delaunay condition.
         * As a shortcut, if the given edge has no adjacent (is on the hull), it is
         * certainly Delaunay.
         *
         * @private
         * @param edg The edge shared by the triangles to test.
         * @return True if they are Delaunay.
         */
        isDelaunay(edg) {
            const { triangles, halfedges } = this.del, adj = halfedges[edg];
            if (adj === -1) {
                return true;
            }
            const p1 = triangles[prevEdge(edg)], p2 = triangles[edg], p3 = triangles[nextEdge(edg)], px = triangles[prevEdge(adj)];
            return !this.inCircle(p1, p2, p3, px);
        }
        /**
         * Update the vertex -> incoming edge map.
         *
         * @private
         * @param start The id of an *outgoing* edge.
         * @return The id of the right-most incoming edge.
         */
        updateVert(start) {
            const { triangles, halfedges } = this.del, vm = this.vertMap, v = triangles[start];
            // When iterating over incoming edges around a vertex, we do so in
            // clockwise order ('going left'). If the vertex lies on the hull, two
            // of the edges will have no opposite, leaving a gap. If the starting
            // incoming edge is not the right-most, we will miss edges between it
            // and the gap. So walk counter-clockwise until we find an edge on the
            // hull, or get back to where we started.
            let inc = prevEdge(start), adj = halfedges[inc];
            while (adj !== -1 && adj !== start) {
                inc = prevEdge(adj);
                adj = halfedges[inc];
            }
            vm[v] = inc;
            return inc;
        }
        /**
         * Whether the segment between [p1, p2] intersects with [p3, p4]. When the
         * segments share an end-point (e.g. p1 == p3 etc.), they are not considered
         * intersecting.
         *
         * @private
         * @param p1 The index of point 1 into this.del.coords.
         * @param p2 The index of point 2 into this.del.coords.
         * @param p3 The index of point 3 into this.del.coords.
         * @param p4 The index of point 4 into this.del.coords.
         * @return True if the segments intersect.
         */
        intersectSegments(p1, p2, p3, p4) {
            const pts = this.del.coords;
            // If the segments share one of the end-points, they cannot intersect
            // (provided the input is properly segmented, and the triangulation is
            // correct), but intersectSegments will say that they do. We can catch
            // it here already.
            if (p1 === p3 || p1 === p4 || p2 === p3 || p2 === p4) {
                return false;
            }
            return intersectSegments(pts[p1 * 2], pts[p1 * 2 + 1], pts[p2 * 2], pts[p2 * 2 + 1], pts[p3 * 2], pts[p3 * 2 + 1], pts[p4 * 2], pts[p4 * 2 + 1]);
        }
        /**
         * Whether point px is in the circumcircle of the triangle formed by p1, p2,
         * and p3 (which are in counter-clockwise order).
         *
         * @param p1 The index of point 1 into this.del.coords.
         * @param p2 The index of point 2 into this.del.coords.
         * @param p3 The index of point 3 into this.del.coords.
         * @param px The index of point x into this.del.coords.
         * @return True if (px, py) is in the circumcircle.
         */
        inCircle(p1, p2, p3, px) {
            const pts = this.del.coords;
            return incircle(pts[p1 * 2], pts[p1 * 2 + 1], pts[p2 * 2], pts[p2 * 2 + 1], pts[p3 * 2], pts[p3 * 2 + 1], pts[px * 2], pts[px * 2 + 1]) < 0.0;
        }
        /**
         * Whether point p1, p2, and p are collinear.
         *
         * @private
         * @param p1 The index of segment point 1 into this.del.coords.
         * @param p2 The index of segment point 2 into this.del.coords.
         * @param p The index of the point p into this.del.coords.
         * @return True if the points are collinear.
         */
        isCollinear(p1, p2, p) {
            const pts = this.del.coords;
            return orient2d(pts[p1 * 2], pts[p1 * 2 + 1], pts[p2 * 2], pts[p2 * 2 + 1], pts[p * 2], pts[p * 2 + 1]) === 0.0;
        }
    }
    Constrainautor.intersectSegments = intersectSegments;
    /**
     * Compute if two line segments [p1, p2] and [p3, p4] intersect.
     *
     * @name Constrainautor.intersectSegments
     * @source https://github.com/mikolalysenko/robust-segment-intersect
     * @param p1x The x coordinate of point 1 of the first segment.
     * @param p1y The y coordinate of point 1 of the first segment.
     * @param p2x The x coordinate of point 2 of the first segment.
     * @param p2y The y coordinate of point 2 of the first segment.
     * @param p3x The x coordinate of point 1 of the second segment.
     * @param p3y The y coordinate of point 1 of the second segment.
     * @param p4x The x coordinate of point 2 of the second segment.
     * @param p4y The y coordinate of point 2 of the second segment.
     * @return True if the line segments intersect.
     */
    function intersectSegments(p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y) {
        const x0 = orient2d(p1x, p1y, p3x, p3y, p4x, p4y), y0 = orient2d(p2x, p2y, p3x, p3y, p4x, p4y);
        if ((x0 > 0 && y0 > 0) || (x0 < 0 && y0 < 0)) {
            return false;
        }
        const x1 = orient2d(p3x, p3y, p1x, p1y, p2x, p2y), y1 = orient2d(p4x, p4y, p1x, p1y, p2x, p2y);
        if ((x1 > 0 && y1 > 0) || (x1 < 0 && y1 < 0)) {
            return false;
        }
        //Check for degenerate collinear case
        if (x0 === 0 && y0 === 0 && x1 === 0 && y1 === 0) {
            return !(Math.max(p3x, p4x) < Math.min(p1x, p2x) ||
                Math.max(p1x, p2x) < Math.min(p3x, p4x) ||
                Math.max(p3y, p4y) < Math.min(p1y, p2y) ||
                Math.max(p1y, p2y) < Math.min(p3y, p4y));
        }
        return true;
    }

    return Constrainautor;

}));
