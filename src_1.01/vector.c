/*----------------------------------------------------------------------
  SerialReax - Reax Force Field Simulator
      
  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu
 
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of 
  the License, or (at your option) any later version.
               
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "vector.h"


inline int Vector_isZero( real* v, int k )
{
  for( --k; k>=0; --k )
    if( fabs( v[k] ) > ALMOST_ZERO )
      return 0;
  
  return 1;
}


inline void Vector_MakeZero( real *v, int k )
{
  for( --k; k>=0; --k )
    v[k] = 0;
}


inline void Vector_Copy( real* dest, real* v, int k )
{
  for( --k; k>=0; --k )
    dest[k] = v[k];
}


inline void Vector_Scale( real* dest, real c, real* v, int k )
{
  for( --k; k>=0; --k )
    dest[k] = c * v[k];
}


inline void Vector_Sum( real* dest, real c, real* v, real d, real* y, int k )
{
  for( --k; k>=0; --k )
    dest[k] = c * v[k] + d * y[k];
}


inline void Vector_Add( real* dest, real c, real* v, int k )
{
  for( --k; k>=0; --k )
    dest[k] += c * v[k];
}


void Vector_Print( FILE *fout, char *vname, real *v, int k )
{
  int i;

  fprintf( fout, "%s:\n", vname );
  for( i = 0; i < k; ++i )
    fprintf( fout, "%24.15e\n", v[i] );
  fprintf( fout, "\n" );
}


inline real Dot( real* v1, real* v2, int k )
{
  real ret = 0;
  
  for( --k; k>=0; --k )
    ret +=  v1[k] * v2[k];

  return ret;
}


inline real Norm( real* v1, int k )
{
  real ret = 0;
  
  for( --k; k>=0; --k )
    ret +=  SQR( v1[k] );

  return SQRT( ret );
}


inline void rvec_Copy( rvec dest, rvec src )
{
  dest[0] = src[0], dest[1] = src[1], dest[2] = src[2];
}

inline void rvec_Scale( rvec ret, real c, rvec v )
{
  ret[0] = c * v[0], ret[1] = c * v[1], ret[2] = c * v[2];
}


inline void rvec_Add( rvec ret, rvec v )
{
  ret[0] += v[0], ret[1] += v[1], ret[2] += v[2];
}


inline void rvec_ScaledAdd( rvec ret, real c, rvec v )
{
  ret[0] += c * v[0], ret[1] += c * v[1], ret[2] += c * v[2];
}


inline void rvec_Sum( rvec ret, rvec v1 ,rvec v2 )
{
  ret[0] = v1[0] + v2[0];
  ret[1] = v1[1] + v2[1];
  ret[2] = v1[2] + v2[2];
}


inline void rvec_ScaledSum( rvec ret, real c1, rvec v1 ,real c2, rvec v2 )
{
  ret[0] = c1 * v1[0] + c2 * v2[0]; 
  ret[1] = c1 * v1[1] + c2 * v2[1];
  ret[2] = c1 * v1[2] + c2 * v2[2];
}


inline real rvec_Dot( rvec v1, rvec v2 )
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}


inline real rvec_ScaledDot( real c1, rvec v1, real c2, rvec v2 )
{
  return (c1*c2) * (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}


inline void rvec_Multiply( rvec r, rvec v1, rvec v2 )
{
  r[0] = v1[0] * v2[0];
  r[1] = v1[1] * v2[1];
  r[2] = v1[2] * v2[2];
}


inline void rvec_iMultiply( rvec r, ivec v1, rvec v2 )
{
  r[0] = v1[0] * v2[0];
  r[1] = v1[1] * v2[1];
  r[2] = v1[2] * v2[2];
}


inline void rvec_Divide( rvec r, rvec v1, rvec v2 )
{
  r[0] = v1[0] / v2[0];
  r[1] = v1[1] / v2[1];
  r[2] = v1[2] / v2[2];
}


inline void rvec_iDivide( rvec r, rvec v1, ivec v2 )
{
  r[0] = v1[0] / v2[0];
  r[1] = v1[1] / v2[1];
  r[2] = v1[2] / v2[2];
}


inline void rvec_Invert( rvec r, rvec v )
{
  r[0] = 1. / v[0];
  r[1] = 1. / v[1];
  r[2] = 1. / v[2];
}


inline void rvec_Cross( rvec ret, rvec v1, rvec v2 )
{
  ret[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ret[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ret[2] = v1[0] * v2[1] - v1[1] * v2[0];
}


inline void rvec_OuterProduct( rtensor r, rvec v1, rvec v2 )
{
  int i, j;

  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
      r[i][j] = v1[i] * v2[j];
}


inline real rvec_Norm_Sqr( rvec v )
{
  return SQR(v[0]) + SQR(v[1]) + SQR(v[2]);
}


inline real rvec_Norm( rvec v )
{
  return SQRT( SQR(v[0]) + SQR(v[1]) + SQR(v[2]) );
}


inline int rvec_isZero( rvec v )
{
  if( fabs(v[0]) > ALMOST_ZERO || 
      fabs(v[1]) > ALMOST_ZERO || 
      fabs(v[2]) > ALMOST_ZERO )
    return 0;
  return 1;
}


inline void rvec_MakeZero( rvec v )
{
  v[0] = v[1] = v[2] = ZERO;
}


inline void rvec_Random( rvec v )
{
  v[0] = Random(2.0)-1.0;
  v[1] = Random(2.0)-1.0;
  v[2] = Random(2.0)-1.0;
}


inline void rtensor_Multiply( rtensor ret, rtensor m1, rtensor m2 )
{
  int i, j, k;
  rtensor temp;

  // check if the result matrix is the same as one of m1, m2.
  // if so, we cannot modify the contents of m1 or m2, so 
  // we have to use a temp matrix.
  if( ret == m1 || ret == m2 )
    {
      for( i = 0; i < 3; ++i )
	for( j = 0; j < 3; ++j )
	  {
	    temp[i][j] = 0;	    
	    for( k = 0; k < 3; ++k )
	      temp[i][j] += m1[i][k] * m2[k][j];
	  }
      
      for( i = 0; i < 3; ++i )
	for( j = 0; j < 3; ++j )
	  ret[i][j] = temp[i][j];	
    }
  else
    {
      for( i = 0; i < 3; ++i )
	for( j = 0; j < 3; ++j )
	  {
	    ret[i][j] = 0;	    
	    for( k = 0; k < 3; ++k )
	      ret[i][j] += m1[i][k] * m2[k][j];
	  }
    }
}


inline void rtensor_MatVec( rvec ret, rtensor m, rvec v )
{
  int i;
  rvec temp;

  // if ret is the same vector as v, we cannot modify the 
  // contents of v until all computation is finished.
  if( ret == v )
    {
      for( i = 0; i < 3; ++i )
	temp[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];

      for( i = 0; i < 3; ++i )
	ret[i] = temp[i];
    }
  else
    {
      for( i = 0; i < 3; ++i )
	ret[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
    }
}


inline void rtensor_Scale( rtensor ret, real c, rtensor m )
{
  int i, j;

  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
      ret[i][j] = c * m[i][j];
}


inline void rtensor_Add( rtensor ret, rtensor t )
{
  int i, j;

  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
      ret[i][j] += t[i][j];
}


inline void rtensor_ScaledAdd( rtensor ret, real c, rtensor t )
{
  int i, j;

  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
      ret[i][j] += c * t[i][j];
}


inline void rtensor_Sum( rtensor ret, rtensor t1, rtensor t2 )
{
  int i, j;

  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
      ret[i][j] = t1[i][j] + t2[i][j];
}


inline void rtensor_ScaledSum( rtensor ret, real c1, rtensor t1, 
			       real c2, rtensor t2 )
{
  int i, j;

  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
      ret[i][j] = c1 * t1[i][j] + c2 * t2[i][j];
}


inline void rtensor_Copy( rtensor ret, rtensor t )
{
  int i, j;

  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
      ret[i][j] = t[i][j];
}


inline void rtensor_Identity( rtensor t )
{
  t[0][0] = t[1][1] = t[2][2] = 1;
  t[0][1] = t[0][2] = t[1][0] = t[1][2] = t[2][0] = t[2][1] = ZERO;
}


inline void rtensor_MakeZero( rtensor t )
{
  t[0][0] = t[0][1] = t[0][2] = ZERO;
  t[1][0] = t[1][1] = t[1][2] = ZERO;
  t[2][0] = t[2][1] = t[2][2] = ZERO;
}


inline void rtensor_Transpose( rtensor ret, rtensor t )
{
  ret[0][0] = t[0][0], ret[1][1] = t[1][1], ret[2][2] = t[2][2];
  ret[0][1] = t[1][0], ret[0][2] = t[2][0];
  ret[1][0] = t[0][1], ret[1][2] = t[2][1];
  ret[2][0] = t[0][2], ret[2][1] = t[1][2];
}


inline real rtensor_Det( rtensor t )
{
  return ( t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1] ) +
           t[0][1] * (t[1][2] * t[2][0] - t[1][0] * t[2][2] ) +
           t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0] ) );
}


inline real rtensor_Trace( rtensor t )
{
  return (t[0][0] + t[1][1] + t[2][2]);
}


void Print_rTensor(FILE* fp, rtensor t)
{
  int i, j;

  for (i=0; i < 3; i++)
    {
      fprintf(fp,"[");
      for (j=0; j < 3; j++)
	fprintf(fp,"%8.3f,\t",t[i][j]);
      fprintf(fp,"]\n");
    }
}


inline void ivec_MakeZero( ivec v )
{
  v[0] = v[1] = v[2] = 0;
}


inline void ivec_Copy( ivec dest, ivec src )
{
  dest[0] = src[0], dest[1] = src[1], dest[2] = src[2];
}


inline void ivec_Scale( ivec dest, real C, ivec src )
{
  dest[0] = C * src[0];
  dest[1] = C * src[1];
  dest[2] = C * src[2];
}


inline void ivec_rScale( ivec dest, real C, rvec src )
{
  dest[0] = (int)(C * src[0]);
  dest[1] = (int)(C * src[1]);
  dest[2] = (int)(C * src[2]);
}


inline int ivec_isZero( ivec v )
{
  if( v[0]==0 && v[1]==0 && v[2]==0 )
    return 1;
  return 0;
}


inline int ivec_isEqual( ivec v1, ivec v2 )
{
  if( v1[0]==v2[0] && v1[1]==v2[1] && v1[2]==v2[2] )
    return 1;

  return 0;
}


inline void ivec_Sum( ivec dest, ivec v1, ivec v2 )
{
  dest[0] = v1[0] + v2[0];
  dest[1] = v1[1] + v2[1];
  dest[2] = v1[2] + v2[2];
}
