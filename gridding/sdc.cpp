/*
 * Copyright (c) 2014, Dignity Health
 * 
 *     The GPI core node library is licensed under
 * either the BSD 3-clause or the LGPL v. 3.
 * 
 *     Under either license, the following additional term applies:
 * 
 *         NO CLINICAL USE.  THE SOFTWARE IS NOT INTENDED FOR COMMERCIAL
 * PURPOSES AND SHOULD BE USED ONLY FOR NON-COMMERCIAL RESEARCH PURPOSES.  THE
 * SOFTWARE MAY NOT IN ANY EVENT BE USED FOR ANY CLINICAL OR DIAGNOSTIC
 * PURPOSES.  YOU ACKNOWLEDGE AND AGREE THAT THE SOFTWARE IS NOT INTENDED FOR
 * USE IN ANY HIGH RISK OR STRICT LIABILITY ACTIVITY, INCLUDING BUT NOT LIMITED
 * TO LIFE SUPPORT OR EMERGENCY MEDICAL OPERATIONS OR USES.  LICENSOR MAKES NO
 * WARRANTY AND HAS NOR LIABILITY ARISING FROM ANY USE OF THE SOFTWARE IN ANY
 * HIGH RISK OR STRICT LIABILITY ACTIVITIES.
 * 
 *     If you elect to license the GPI core node library under the LGPL the
 * following applies:
 * 
 *         This file is part of the GPI core node library.
 * 
 *         The GPI core node library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version. GPI core node library is distributed
 * in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 * 
 *         You should have received a copy of the GNU Lesser General Public
 * License along with the GPI core node library. If not, see
 * <http://www.gnu.org/licenses/>.
 */


/*******************************/
/* 1D SDC */
/*******************************/
int onedsdc(Array<double> &crds_in, Array<double> &wates, Array<double> &sdc,
            Array<double> &cmtx, long numiter, double taper0)
{
  long i;
  double c0,c1,c2,c3,c4,c5;
  double x, x2,x3,x4,x5;
  long ii0,ii1,ii,ki;

  long npts = wates.size(0);
  Array<double> taper(npts);
  Array<double> crds(crds_in); // make a copy

/////////////////
// 1. make kernel
/////////////////

  long kernmax = 1000;
  Array<double> kernel(kernmax+1);

// the kernel, convolved with itself, should be equivalent to FT of a circle convolved with itself
// See MRM paper by Zwart, Pipe
// Doing this offline, can fit the result to a 5th order polynomial k(x), x from 0 to 1
// NOTE THIS IS THE KERNEL for 2D, NEED TO OPTIMIZE FOR 1D (?)
  c0 = 1.;
  c1 = 0.03056504;
  c2 = -3.01961845;
  c3 = 0.6679865;
  c4 = 2.77924058;
  c5 = -1.45923643;
  for (i = 0; i <= kernmax; i++) {
    x = double(i)/double(kernmax);
    x2 = x*x;
    x4 = x2*x2;
    x3 = x*x2;
    x5 = x*x4;
    kernel(i) = c0 + c1*x + c2*x2 + c3*x3 + c4*x4 + c5*x5;
    }

///////////
// 2. Calculate taper mask
//    Reassign coord range from {-0.5,+0.5} to 3, mtxdim-3
//    cmtx size is (2*effmtx)+6
//                  The 2 is for 2x oversampling per Zwart
//                  The 6 is a safety factor (3 on all sides) to avoid checking for array bounds
///////////

  double mtxsize = cmtx.size(0)-6;
  double maxcrd,crdedge,crdmag, crdnorm;
  double xcrd;

  // Taper Logic
  maxcrd = 0;
  for (i = 0; i < npts; i++) {
    if (fabs(crds(0,i)) > maxcrd) {
      maxcrd = fabs(crds(0,i));
    } } // if, i
  crdedge = (1-taper0)*maxcrd;
  if (taper0 > 0.) {
    crdnorm = 1./(maxcrd-crdedge);
    }
  else {
    crdnorm = 0;
    }
  for (i = 0; i < npts; i++) {
    if (fabs(crds(0,i)) > crdedge) {
      crdmag = fabs(crds(0,i));
      taper(i) = (maxcrd-crdmag)*crdnorm;
      }
    else {
      taper(i) = 1.;
      }
    } // i
    
  // now scale coords, making sure to start within +/- 0.5
  for (i = 0; i < npts; i++) {
    xcrd = min(max(-.5,crds(0,i)),0.5);
    crds(0,i) = 3.+((xcrd+0.5)*mtxsize);
    }

/////////////
// 3. Iteratively convolve onto and off of cmtx
/////////////

  double krad = 1.5*1.65; // the width of fitted kernel radius in cmtx coordinates
  double knorm = kernmax/krad;
  double sdcdenom;
  double dx;
  long   count;

  // Initialize weights by input relative weights
  for (i = 0; i < npts; i++) {
    sdc(i) = wates(i);
    }

  for (count=0;count<numiter;count++) {

      // Initialize cmtx
    for (ii = 0; ii < cmtx.size(0); ii++) {
      cmtx(ii) = 0;
      }

    // crds -> cmtx
    for (i = 0; i < npts; i++) {
      ii0 =  ceil(crds(0,i)-krad);
      ii1 = floor(crds(0,i)+krad);
      for (ii = ii0; ii <= ii1; ii++) {
        dx = fabs(crds(0,i)-double(ii));
        ki = floor(knorm*dx);
        cmtx(ii) += sdc(i)*kernel(ki);
        } // ii
      } // i

    // cmtx -> crds
    for (i = 0; i < npts; i++) {
      sdcdenom = 0.;
      ii0 =  ceil(crds(0,i)-krad);
      ii1 = floor(crds(0,i)+krad);
      for (ii = ii0; ii <= ii1; ii++) {
        dx = fabs(crds(0,i)-double(ii));
        ki = floor(knorm*dx);
        sdcdenom += cmtx(ii)*kernel(ki);
        } // ii
      if (sdcdenom > 0)
        sdc(i) *= taper(i)/sdcdenom;
      } // i

    } // count

   return (0);

} // onedsdc()

/*******************************/
/* 2D SDC */
/*******************************/
int twodsdc(Array<double> &crds, Array<double> &wates, Array<double> &sdc,
            Array<double> &cmtx, long numiter, double taper0)
{
  long i;
  double c0,c1,c2,c3,c4,c5;
  double x, x2,x3,x4,x5;
  long ii0,ii1,jj0,jj1,ii,jj,ki;

  long npts = wates.size(0);
  Array<double> taper(npts);

/////////////////
// 1. make kernel
/////////////////
  long kernmax = 1000;
  Array<double> kernel(kernmax+1);

// the kernel, convolved with itself, should be equivalent to FT of a circle convolved with itself
// See MRM paper by Zwart, Pipe
// Doing this offline, can fit the result to a 5th order polynomial k(x), x from 0 to 1
  c0 = 1.;
  c1 = 0.03056504;
  c2 = -3.01961845;
  c3 = 0.6679865;
  c4 = 2.77924058;
  c5 = -1.45923643;
  for (i = 0; i <= kernmax; i++) {
    x = double(i)/double(kernmax);
    x2 = x*x;
    x4 = x2*x2;
    x3 = x*x2;
    x5 = x*x4;
    kernel(i) = c0 + c1*x + c2*x2 + c3*x3 + c4*x4 + c5*x5;
    }
///////////
// 2. Calculate taper mask
//    Reassign coord range from {-0.5,+0.5} to 3, mtxdim-3
//    cmtx size is (2*effmtx)+6
//                  The 2 is for 2x oversampling per Zwart
//                  The 6 is a safety factor (3 on all sides) to avoid checking for array bounds
///////////
  double mtxsize = cmtx.size(0)-6;
  double maxcrd,crdedge,crdedge2,crdmag, crdnorm;
  double xcrd,ycrd;
  // Taper Logic
  maxcrd = 0;
  for (i = 0; i < npts; i++) {
    if (crds(0,i)*crds(0,i) + crds(1,i)*crds(1,i) > maxcrd) {
      maxcrd = crds(0,i)*crds(0,i) + crds(1,i)*crds(1,i);
    } } // if, i
  maxcrd = sqrt(maxcrd);
  crdedge = (1-taper0)*maxcrd;
  crdedge2 = crdedge*crdedge;
  if (taper0 > 0.) {
    crdnorm = 1./(maxcrd-crdedge);
    }
  else {
    crdnorm = 0;
    }
  for (i = 0; i < npts; i++) {
    if (crds(0,i)*crds(0,i) + crds(1,i)*crds(1,i) > crdedge2) {
      crdmag = sqrt(crds(0,i)*crds(0,i) + crds(1,i)*crds(1,i));
      taper(i) = (maxcrd-crdmag)*crdnorm;
      }
    else {
      taper(i) = 1.;
      }
    } // i
  // now scale coords, making sure to start within +/- 0.5
  for (i = 0; i < npts; i++) {
    xcrd = min(max(-.5,crds(0,i)),0.5);
    ycrd = min(max(-.5,crds(1,i)),0.5);
    crds(0,i) = 3.+((xcrd+0.5)*mtxsize);
    crds(1,i) = 3.+((ycrd+0.5)*mtxsize);
    }
/////////////
// 3. Iteratively convolve onto and off of cmtx
/////////////

  double krad = 1.5*1.65; // the width of fitted kernel radius in cmtx coordinates
  double krad2 = krad*krad;
  double knorm = kernmax/krad;
  double sdcdenom;
  double dx,dx2,dy,rad2;
  long   count;

  // Initialize weights by input relative weights
  for (i = 0; i < npts; i++) {
    sdc(i) = wates(i);
    }
  for (count=0;count<numiter;count++) {

      // Initialize cmtx
    for (ii = 0; ii < cmtx.size(0); ii++) {
      for (jj = 0; jj < cmtx.size(1); jj++) {
        cmtx(ii,jj) = 0;
      } } // ii, jj

    // crds -> cmtx
    for (i = 0; i < npts; i++) {
      ii0 =  ceil(crds(0,i)-krad);
      ii1 = floor(crds(0,i)+krad);
      jj0 =  ceil(crds(1,i)-krad);
      jj1 = floor(crds(1,i)+krad);
      for (ii = ii0; ii <= ii1; ii++) {
        dx = crds(0,i)-double(ii);
        dx2 = dx*dx;
        for (jj = jj0; jj <= jj1; jj++) {
          dy = crds(1,i)-double(jj);
          rad2 = dx2+dy*dy;
          if (rad2 < krad2) {
            ki = floor(knorm*sqrt(rad2));
            cmtx(ii,jj) += sdc(i)*kernel(ki);
          } //if rad2 < krad2
        } } // ii,jj
      } // i

    // cmtx -> crds
    for (i = 0; i < npts; i++) {
      sdcdenom = 0.;
      ii0 =  ceil(crds(0,i)-krad);
      ii1 = floor(crds(0,i)+krad);
      jj0 =  ceil(crds(1,i)-krad);
      jj1 = floor(crds(1,i)+krad);
      for (ii = ii0; ii <= ii1; ii++) {
        dx = crds(0,i)-double(ii);
        dx2 = dx*dx;
        for (jj = jj0; jj <= jj1; jj++) {
          dy = crds(1,i)-double(jj);
          rad2 = dx2+dy*dy;
          if (rad2 < krad2) {
            ki = floor(knorm*sqrt(rad2));
            sdcdenom += cmtx(ii,jj)*kernel(ki);
            } //if rad2 < krad2
        } } // ii,jj
      if (sdcdenom > 0)
        sdc(i) *= taper(i)/sdcdenom;
      } // i

    } // count
   return (0);

} // twodsdc()

/*******************************/
/* 3D SDC */
/*******************************/
int threedsdc(Array<double> &crds, Array<double> &wates, Array<double> &sdc,
            Array<double> &cmtx, long numiter, double taper0)
{
  long i;
  double c0,c1,c2,c3,c4,c5;
  double x, x2,x3,x4,x5;
  long ii0,ii1,jj0,jj1,kk0,kk1,ii,jj,kk,ki;

  long npts = wates.size(0);
  Array<double> taper(npts);

/////////////////
// 1. make kernel
/////////////////

  long kernmax = 1000;
  Array<double> kernel(kernmax+1);

// the kernel, convolved with itself, should be equivalent to FT of a sphere convolved with itself
// See MRM paper by Zwart, Pipe
// Doing this offline, can fit the result to a 5th order polynomial k(x), x from 0 to 1
  c0 = 1.;
  c1 = 0.04522831;
  c2 = -3.36020304;
  c3 = 1.12417012;
  c4 = 2.82448025;
  c5 = -1.63447764;
  for (i = 0; i <= kernmax; i++) {
    x = double(i)/double(kernmax);
    x2 = x*x;
    x4 = x2*x2;
    x3 = x*x2;
    x5 = x*x4;
    kernel(i) = c0 + c1*x + c2*x2 + c3*x3 + c4*x4 + c5*x5;
    }

///////////
// 2. Calculate taper mask
//    Reassign coord range from {-0.5,+0.5} to 3, mtxdim-3
//    cmtx size is (2*effmtx)+6
//                  The 2 is for 2x oversampling per Zwart
//                  The 6 is a safety factor (3 on all sides) to avoid checking for array bounds
///////////

  double mtxxsize = cmtx.size(0)-6;
  double mtxzsize = cmtx.size(2)-6;
  double maxcrd,crdedge,crdedge2,crdmag, crdnorm;
  double xcrd,ycrd,zcrd;

  // Taper Logic
  maxcrd = 0;
  for (i = 0; i < npts; i++) {
    if (crds(0,i)*crds(0,i) + crds(1,i)*crds(1,i) + crds(2,i)*crds(2,i) > maxcrd) {
      maxcrd = crds(0,i)*crds(0,i) + crds(1,i)*crds(1,i) + crds(2,i)*crds(2,i);
    } } // if, i
  maxcrd = sqrt(maxcrd);
  crdedge = (1.-taper0)*maxcrd;
  crdedge2 = crdedge*crdedge;
  if (taper0 > 0.) {
    crdnorm = 1./(maxcrd-crdedge);
    }
  else {
    crdnorm = 0;
    }
  for (i = 0; i < npts; i++) {
    if (crds(0,i)*crds(0,i) + crds(1,i)*crds(1,i) + crds(2,i)*crds(2,i) > crdedge2) {
      crdmag = sqrt(crds(0,i)*crds(0,i) + crds(1,i)*crds(1,i) + crds(2,i)*crds(2,i));
      taper(i) = (maxcrd-crdmag)*crdnorm;
      }
    else {
      taper(i) = 1.;
      }
    } // i
    
  // now scale coords, making sure to start within +/- 0.5
  for (i = 0; i < npts; i++) {
    xcrd = min(max(-.5,crds(0,i)),0.5);
    ycrd = min(max(-.5,crds(1,i)),0.5);
    zcrd = min(max(-.5,crds(2,i)),0.5);
    crds(0,i) = 3.+((xcrd+0.5)*mtxxsize);
    crds(1,i) = 3.+((ycrd+0.5)*mtxxsize);
    crds(2,i) = 3.+((zcrd+0.5)*mtxzsize);
    }

/////////////
// 3. Iteratively convolve onto and off of cmtx
/////////////

  double krad = 1.5*1.93; // the width of fitted kernel radius in cmtx coordinates
  double krad2 = krad*krad;
  double knorm = kernmax/krad;
  double sdcdenom;
  double dx,dy,dz,dxy2,rad2;
  long   count;
  double dx2,dy2[20],dz2[20];

  // Initialize weights by input relative weights
  for (i = 0; i < npts; i++) {
    sdc(i) = wates(i);
    }

  for (count=0;count<numiter;count++) {

    // Initialize
    for (ii = 0; ii < cmtx.size(0); ii++) {
      for (jj = 0; jj < cmtx.size(1); jj++) {
        for (kk = 0; kk < cmtx.size(2); kk++) {
          cmtx(ii,jj,kk) = 0;
      } } } // ii, jj, kk

    //////////////////
    // crds -> cmtx //
    //////////////////
    for (i = 0; i < npts; i++) {
      ii0 =  ceil(crds(0,i)-krad);
      ii1 = floor(crds(0,i)+krad);
      jj0 =  ceil(crds(1,i)-krad);
      jj1 = floor(crds(1,i)+krad);
      kk0 =  ceil(crds(2,i)-krad);
      kk1 = floor(crds(2,i)+krad);

      // dx, dy, and dz are the respective distances from point "i" to ii,jj,kk in cmtx
      // We will compare dx^2 + dy^2 + dz^2 to the radius^2 of the kernel to see if we need to evaluate
      // calculate the distances^2 in y and z direction first (faster)
      for (jj = jj0; jj <= jj1; jj++) {
        dy = crds(1,i)-double(jj);
        dy2[jj-jj0] = dy*dy;
        }
      for (kk = kk0; kk <= kk1; kk++) {
        dz = crds(2,i)-double(kk);
        dz2[kk-kk0] = dz*dz;
        }

      // OK Now evaluate every potential point in cmtx and broadcast to ii,jj,kk if appropriate
      for (ii = ii0; ii <= ii1; ii++) {
        dx = crds(0,i)-double(ii);
        dx2 = dx*dx;
        for (jj = jj0; jj <= jj1; jj++) {
          dxy2 = dx2 + dy2[jj-jj0];
          if (dxy2 < krad2) {
            for (kk = kk0; kk <= kk1; kk++) {
              rad2 = dxy2+dz2[kk-kk0];
              if (rad2 < krad2) {
                ki = floor(knorm*sqrt(rad2));
                cmtx(ii,jj,kk) += sdc(i)*kernel(ki);
                } //if rad2 < krad2
            } } // if dxy2, kk
        } } // ii,jj
      } // i

    //////////////////
    // cmtx -> crds //
    //////////////////
    for (i = 0; i < npts; i++) {
      sdcdenom = 0;
      ii0 =  ceil(crds(0,i)-krad);
      ii1 = floor(crds(0,i)+krad);
      jj0 =  ceil(crds(1,i)-krad);
      jj1 = floor(crds(1,i)+krad);
      kk0 =  ceil(crds(2,i)-krad);
      kk1 = floor(crds(2,i)+krad);
      for (jj = jj0; jj <= jj1; jj++) {
        dy = crds(1,i)-double(jj);
        dy2[jj-jj0] = dy*dy;
        }
      for (kk = kk0; kk <= kk1; kk++) {
        dz = crds(2,i)-double(kk);
        dz2[kk-kk0] = dz*dz;
        }
      for (ii = ii0; ii <= ii1; ii++) {
        dx = crds(0,i)-double(ii);
        dx2 = dx*dx;
        for (jj = jj0; jj <= jj1; jj++) {
          dxy2 = dx2 + dy2[jj-jj0];
          if (dxy2 < krad2) {
            for (kk = kk0; kk <= kk1; kk++) {
              rad2 = dxy2+dz2[kk-kk0];
              if (rad2 < krad2) {
                ki = floor(knorm*sqrt(rad2));
                sdcdenom += cmtx(ii,jj,kk)*kernel(ki);
                } //if rad2 < krad2
            } } // if dxy2, kk
        } } // ii,jj
      if (sdcdenom > 0) {
        sdc(i) *= taper(i)/sdcdenom;
        }
      } // i

    } // count

   return (0);

} // threedsdc()
